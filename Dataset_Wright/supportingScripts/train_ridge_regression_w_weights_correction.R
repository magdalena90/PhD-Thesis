
library(tidyverse)
library(caret) ; library(DMwR) ; library(MLmetrics) ; library(ROCR)
library(mgcv)

#setwd('~/PhD/Thesis/Dataset_Wright/supportingScripts')

# Functions
create_train_val_test_sets = function(dataset, genes_info, p_train, p_val, seed){
  
  # Get SFARI Score of all the samples so our train and test sets are balanced for each score
  sample_scores = data.frame('ID' = rownames(dataset)) %>% 
    left_join(genes_info %>% dplyr::select(ID, gene.score), by = 'ID') %>% 
    mutate(gene.score = ifelse(!(gene.score %in% c('1','2','3')), 'None', gene.score))
  
  set.seed(seed)
  train_idx = createDataPartition(sample_scores$gene.score, p = p_train, list = FALSE)
  train_set = dataset[train_idx,]
  
  val_idx = createDataPartition(sample_scores$gene.score[-train_idx], p = p_val/(1-p_train), list = FALSE)
  val_set = dataset[-train_idx,][val_idx,]
  test_set = dataset[!rownames(dataset) %in% c(rownames(train_set), rownames(val_set)),]
  
  # Modify SFARI label in train set, save gene IDs (bc we lose them with SMOTE) and perform oversampling using SMOTE
  # Note: we can't return the IDs to the training set because of the oversampling, which created repeated IDs
  set.seed(seed)
  train_set = train_set %>% mutate(SFARI = ifelse(SFARI == TRUE, 'SFARI', 'not_SFARI') %>% as.factor,
                                   ID = rownames(.) %>% as.factor) %>% SMOTE(form = SFARI ~ . - ID)
  train_set_IDs = train_set %>% pull(ID)
  
  return(list('train_set' = train_set %>% dplyr::select(-ID), 'val_set' = val_set, 'test_set' = test_set, 
              'train_set_IDs' = train_set_IDs))
}

calc_performance_metrics = function(model, selected_set){
  
  # Modification to keep the predictions with the same mean as the original model
  predictions = model %>% predict(selected_set, type = 'prob') %>% 
    mutate(prob = SFARI - mean(SFARI) + 0.3909484, pred = prob>0.5, 
           SFARI = selected_set$SFARI) %>% dplyr::select(-c(not_SFARI))
  
  if(all(predictions$pred == 0)){
    prec = NA
    F1 = NA
  } else {
    prec = Precision(predictions$SFARI %>% as.numeric, predictions$pred %>% as.numeric, positive = '1')
    F1 = F1_Score(predictions$SFARI %>% as.numeric, predictions$pred %>% as.numeric, positive = '1')
  }
  
  acc = mean(predictions$SFARI == predictions$pred)
  rec = Recall(predictions$SFARI %>% as.numeric, predictions$pred %>% as.numeric, positive = '1')
  pred_ROCR = prediction(predictions$prob, predictions$SFARI)
  AUC = performance(pred_ROCR, measure='auc')@y.values[[1]]
  MLP = performance(pred_ROCR, measure='lift', x.measure='rpp')@y.values[[1]] %>% max(na.rm = TRUE)
  b_acc = mean(c(mean(predictions$SFARI[predictions$SFARI] == predictions$pred[predictions$SFARI]),
                 mean(predictions$SFARI[!predictions$SFARI] == predictions$pred[!predictions$SFARI])))
  
  return(c('acc' = acc, 'prec' = prec, 'rec' = rec, 'F1' = F1, 'AUC' = AUC, 'MLP'=MLP, 'b_acc' = b_acc))  
}

calc_pm_aggregated_predictions = function(predictions){
  
  if(all(predictions$pred == 0)){
    prec = NA
    F1 = NA
  } else {
    prec = Precision(predictions$SFARI %>% as.numeric, predictions$pred %>% as.numeric, positive = '1')
    F1 = F1_Score(predictions$SFARI %>% as.numeric, predictions$pred %>% as.numeric, positive = '1')
  }
  
  acc = mean(predictions$SFARI == predictions$pred)
  rec = Recall(predictions$SFARI %>% as.numeric, predictions$pred %>% as.numeric, positive = '1')
  pred_ROCR = prediction(predictions$prob, predictions$SFARI)
  AUC = performance(pred_ROCR, measure='auc')@y.values[[1]]
  MLP = performance(pred_ROCR, measure='lift', x.measure='rpp')@y.values[[1]] %>% max(na.rm = TRUE)
  b_acc = mean(c(mean(predictions$SFARI[predictions$SFARI] == predictions$pred[predictions$SFARI]),
                 mean(predictions$SFARI[!predictions$SFARI] == predictions$pred[!predictions$SFARI])))
  
  return(c('acc' = acc, 'prec' = prec, 'rec' = rec, 'F1' = F1, 'AUC' = AUC, 'MLP' = MLP, 'b_acc' = b_acc))
}

calibrate_weights = function(dataset, genes_info, seed, Loops){
  
  # Create training, validation and test sets (but we'll only use the training set)
  train_val_test_sets = create_train_val_test_sets(dataset, genes_info, p_train, p_val, seed)
  train_set = train_val_test_sets[['train_set']]
  train_set_IDs = train_val_test_sets[['train_set_IDs']]
  
  # SET INITIAL PARAMETERS
  # General parameters
  lambda_seq = 10^seq(1, -4, by = -.1)
  k_fold = 10
  cv_repeats = 5
  set.seed(seed)
  trControl = trainControl(method = 'repeatedcv', number = k_fold, repeats = cv_repeats, verboseIter = FALSE, 
                           classProbs = TRUE, savePredictions = 'final', summaryFunction = twoClassSummary,
                           seeds = as.list(seq(seed, seed+length(lambda_seq))))
  # Bias correction parameters
  eta = 0.5
  lambda = 0
  w = rep(1, nrow(train_set))
  
  
  # TRAIN INITIAL MODEL
  set.seed(seed)
  h = train(SFARI ~., data = train_set, method = 'glmnet', trControl = trControl, metric = 'ROC',
            tuneGrid = expand.grid(alpha = 0, lambda = lambda_seq))
  
  
  # CORRECT BIAS
  
  # Mean Expression info
  mean_expr = data.frame('ID' = train_set_IDs) %>% 
              left_join(genes_info %>% dplyr::select(ID, meanExpr), by = 'ID') %>%
              mutate('meanExpr_std' = (meanExpr-mean(meanExpr))/sd(meanExpr))
  
  # Track behaviour of algorithm
  bias_vec = c()
  b_acc_vec = c()
  
  for(l in 1:Loops){
    
    # Calculate bias for positive predicted samples
    probs = predict(h, train_set, type = 'prob')
    probs = probs$SFARI - mean(probs$SFARI) + 0.5#0.3909484
    preds = probs %>% sapply(function(x) ifelse(x < 0.5, 'not_SFARI', 'SFARI'))
    bias = mean(mean_expr$meanExpr_std[preds=='SFARI'])
    if(is.na(bias)) bias = 0 # This happens when all the observations are labelled Negative
    
    # Update weights
    lambda = lambda - eta*bias
    w_hat = exp(lambda*mean_expr$meanExpr_std)
    w = 1/(1+w_hat)
    w[train_set$SFARI=='SFARI'] = w[train_set$SFARI=='SFARI']*w_hat[train_set$SFARI=='SFARI']
    
    # Update tracking vars
    bias_vec = c(bias_vec, bias)
    probs = predict(h, train_set, type = 'prob')
    # Modification to keep the predictions with the same mean as the original model:
    probs = probs$SFARI - mean(probs$SFARI) + 0.5#0.3909484
    preds = probs %>% sapply(function(x) ifelse(x < 0.5, 'not_SFARI', 'SFARI'))
    b_acc = mean(c(mean(preds[preds=='SFARI'] == train_set$SFARI[preds=='SFARI']),
                   mean(preds[preds!='SFARI'] == train_set$SFARI[preds!='SFARI'])))
    b_acc_vec = c(b_acc_vec, b_acc)
    
    # Update h
    set.seed(seed)
    h = train(SFARI ~., data = train_set, method = 'glmnet', weights = w, trControl = trControl,
              metric = 'ROC', tuneGrid = expand.grid(alpha = 0, lambda = lambda_seq))
  }
  
  return(list('lambda' = lambda, 'bias_vec' = bias_vec, 'b_acc_vec' = b_acc_vec))
}

train_unbiased_model = function(dataset, genes_info, p_train, p_val, seed, lambda){
  
  # Create training, validation and test sets
  train_val_test_sets = create_train_val_test_sets(dataset, genes_info, p_train, p_val, seed)
  train_set = train_val_test_sets[['train_set']]
  val_set = train_val_test_sets[['val_set']]
  test_set = train_val_test_sets[['test_set']]
  train_set_IDs = train_val_test_sets[['train_set_IDs']]
  
  # Cross-validation paameters
  lambda_seq = 10^seq(1, -4, by = -.1)
  k_fold = 10
  cv_repeats = 5
  set.seed(seed)
  trControl = trainControl(method = 'repeatedcv', number = k_fold, repeats = cv_repeats, verboseIter = FALSE, 
                           classProbs = TRUE, savePredictions = 'final', summaryFunction = twoClassSummary,
                           seeds = as.list(seq(seed*100, seed*100+length(lambda_seq))))
  
  # Bias correcting parameters
  mean_expr = data.frame('ID' = train_set_IDs) %>% 
              left_join(genes_info %>% dplyr::select(ID, meanExpr), by = 'ID') %>%
              mutate('meanExpr_std' = (meanExpr-mean(meanExpr))/sd(meanExpr))
  w_hat = exp(lambda*mean_expr$meanExpr_std)
  w = 1/(1+w_hat)
  w[train_set$SFARI=='SFARI'] = w[train_set$SFARI=='SFARI']*w_hat[train_set$SFARI=='SFARI']
  
  # Train model
  set.seed(seed)
  model = train(SFARI ~., data = train_set, method = 'glmnet', weights = w, trControl = trControl, 
                metric = 'ROC', tuneGrid = expand.grid(alpha = 0, lambda = lambda_seq))
  
  # Calculate performance metrics
  pred_train = calc_performance_metrics(model, train_set %>% mutate(SFARI = SFARI=='SFARI'))
  pred_val = calc_performance_metrics(model, val_set)
  pred_test = calc_performance_metrics(model, test_set)
  
  # Predict labels in test set
  predictions = model %>% predict(test_set, type = 'prob')
  preds = data.frame('ID' = rownames(test_set), 'prob' = predictions$SFARI) %>% 
    mutate(pred = prob > 0.5)
  
  # Extract coefficients from features
  coefs = coef(model$finalModel, model$bestTune$lambda) %>% as.vector
  
  return(list('coefs'=coefs, 'pm_train'=pred_train, 'pm_val'=pred_val, 'pm_test'=pred_test, 'preds'=preds))
}

# Load dataset
load('../Data/preprocessedData/classification_dataset.RData')

# Define parameters
p_train = 0.7
p_val = 0.15
n_iter = 100
initial_seed = 123
seeds = initial_seed:(initial_seed+n_iter-1)
Loops = 30

# Store predictions
predictions = data.frame('ID' = rownames(dataset), 'SFARI' = dataset$SFARI, 'prob' = 0, 'pred' = 0, 'n' = 0)

# Run model
for(seed in seeds){
  
  print(paste0(seed-initial_seed+1,'/',n_iter))
  
  # Run weights technique to learn the lambda parameter
  weights_model_output = calibrate_weights(dataset, genes_info, seed, Loops)
  lambda = weights_model_output[['lambda']] # Optimised bias correction parameter
  
  # Run unbiased model
  unbiased_model_output = train_unbiased_model(dataset, genes_info, p_train, p_val, seed, lambda)
  
  # Update outputs
  if(seed == initial_seed){
    pm_train = unbiased_model_output$pm_train %>% data.frame
    pm_val = unbiased_model_output$pm_val %>% data.frame
    pm_test = unbiased_model_output$pm_test %>% data.frame
    coefs = unbiased_model_output$coefs %>% data.frame
  } else{
    pm_train = cbind(pm_train, unbiased_model_output$pm_train)
    pm_val = cbind(pm_val, unbiased_model_output$pm_val)
    pm_test = cbind(pm_test, unbiased_model_output$pm_test)
    coefs = cbind(coefs, unbiased_model_output$coefs)
  }
  
  # Update predictions
  update_preds = unbiased_model_output$preds %>% dplyr::select(prob, pred) %>% mutate(n=1)
  predictions[predictions$ID %in% unbiased_model_output$preds$ID, c('prob','pred','n')] = 
    predictions[predictions$ID %in% unbiased_model_output$preds$ID, c('prob','pred','n')] + update_preds
  
}


# Summarise results of each iteration
weights_coefs = data.frame('coef' = c('Intercept', colnames(dataset)[-ncol(dataset)]),
                           'mean' = rowMeans(coefs), 'sd' = apply(coefs,1,sd))
weights_pm_train = data.frame('metric'=rownames(pm_train), 'mean'=rowMeans(pm_train), 'sd'=apply(pm_train,1,sd))
weights_pm_val = data.frame('metric'=rownames(pm_val), 'mean'=rowMeans(pm_val), 'sd'=apply(pm_val,1,sd))
weights_pm_test = data.frame('metric'=rownames(pm_test), 'mean'=rowMeans(pm_test), 'sd'=apply(pm_val,1,sd))

# Calculate performance metrics of the final model
weights_predictions = predictions %>% mutate(prob=prob/n, pred_count=pred-mean(pred)+0.5, pred=prob>0.5) %>%
                      left_join(genes_info %>% dplyr::select(ID, hgnc_symbol), by = 'ID') %>% drop_na(pred)
weights_pm_agg_preds = calc_pm_aggregated_predictions(weights_predictions)

save(weights_predictions, weights_pm_agg_preds, weights_pm_train, weights_pm_val, weights_pm_test, weights_coefs,
     weights_model_output, file = '../Data/Results/Ridge_regression_w_weights_correction.RData')
