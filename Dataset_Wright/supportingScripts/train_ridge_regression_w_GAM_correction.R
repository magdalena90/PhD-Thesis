
library(tidyverse)
library(caret) ; library(DMwR) ; library(MLmetrics) ; library(ROCR)
library(mgcv)

#setwd('~/PhD/Thesis/Dataset_Gandal/supportingScripts')

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

calc_performance_metrics = function(model, selected_set, bias_correction, set_IDs){
  
  predictions = model %>% predict(selected_set, type = 'prob') %>% 
                mutate(prob = SFARI, pred = prob > 0.5, SFARI = selected_set$SFARI) %>% 
                dplyr::select(-c(not_SFARI)) %>% mutate(ID = set_IDs)
  
  # Correct bias in model
  predictions = predictions %>% left_join(bias_correction %>% dplyr::select(-prob), by = 'ID') %>%
                mutate(prob = prob + correction + mean_prob) %>% 
                mutate(prob = ifelse(prob < 0, 0, prob), pred = prob > 0.5)
  
  
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

train_model = function(dataset, genes_info, p_train, p_val, seed, bias_correction){
  
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
  
  # Train Ridge regression
  set.seed(seed)
  model = train(SFARI ~., data = train_set, method = 'glmnet', trControl = trControl, metric = 'ROC', 
                tuneGrid = expand.grid(alpha = 0, lambda = lambda_seq))
  
  # Calculate performance metrics
  pred_train = calc_performance_metrics(model, train_set %>% mutate(SFARI = SFARI=='SFARI'), bias_correction,
                                        train_set_IDs)
  pred_val = calc_performance_metrics(model, val_set, bias_correction, rownames(val_set))
  pred_test = calc_performance_metrics(model, test_set, bias_correction, rownames(test_set))
  
  # Predict labels in test set correcting for bias
  predictions = model %>% predict(test_set, type = 'prob')
  preds = data.frame('ID' = rownames(test_set), 'prob' = predictions$SFARI) %>%  
          left_join(bias_correction %>% dplyr::select(-prob), by = 'ID') %>%
          mutate(prob = prob + correction + mean_prob) %>% 
          mutate(prob = ifelse(prob < 0, 0, prob), pred = prob > 0.5) %>%
          mutate(pred = prob > 0.5)
  
  return(list('pm_train'=pred_train, 'pm_val'=pred_val, 'pm_test'=pred_test, 'preds'=preds))
  
}

# Load dataset
load('../Data/preprocessedData/classification_dataset.RData')

################################################################################################################
# BIAS CORRECTION
load('../Data/Results/Ridge_regression.RData')

bias_correction = ridge_predictions %>% dplyr::select(ID, prob) %>% 
                  left_join(genes_info %>% dplyr::select(ID, meanExpr), by = 'ID')

gam_fit = gam(prob ~ s(meanExpr), method = 'REML', data = bias_correction)

bias_correction = bias_correction %>% mutate(residuals = gam_fit$residuals, 
                                             correction = gam_fit$residuals - prob,
                                             mean_prob = mean(ridge_predictions$prob))
################################################################################################################

# Define parameters
p_train = 0.7
p_val = 0.15
n_iter = 100
initial_seed = 123
seeds = initial_seed:(initial_seed+n_iter-1)

# Store predictions
predictions = data.frame('ID' = rownames(dataset), 'SFARI' = dataset$SFARI, 'prob' = 0, 'pred' = 0, 'n' = 0)

# Run model
for(seed in seeds){
  
  print(paste0(seed-initial_seed+1,'/',n_iter))
  
  # Run model
  model_output = train_model(dataset, genes_info, 0.7, 0.15, seed, bias_correction)
  
  # Update outputs
  if(seed == initial_seed){
    pm_train = model_output$pm_train %>% data.frame
    pm_val = model_output$pm_val %>% data.frame
    pm_test = model_output$pm_test %>% data.frame
  } else{
    pm_train = cbind(pm_train, model_output$pm_train)
    pm_val = cbind(pm_val, model_output$pm_val)
    pm_test = cbind(pm_test, model_output$pm_test)
  }
  
  # Update predictions
  update_preds = model_output$preds %>% dplyr::select(prob, pred) %>% mutate(n=1)
  predictions[predictions$ID %in% model_output$preds$ID, c('prob','pred','n')] = 
    predictions[predictions$ID %in% model_output$preds$ID, c('prob','pred','n')] + update_preds
  
}


# Summarise results of each iteration
GAM_pm_train = data.frame('metric' = rownames(pm_train), 'mean' = rowMeans(pm_train), 'sd' = apply(pm_train,1,sd))
GAM_pm_val = data.frame('metric' = rownames(pm_val), 'mean' = rowMeans(pm_val), 'sd' = apply(pm_val,1,sd))
GAM_pm_test = data.frame('metric' = rownames(pm_test), 'mean' = rowMeans(pm_test), 'sd' = apply(pm_val,1,sd))

# Calculate performance metrics of the final model
GAM_predictions = predictions %>% mutate(prob = prob/n, pred_count = pred, pred = prob>0.5) %>%
  left_join(genes_info %>% dplyr::select(ID, hgnc_symbol), by = 'ID') %>% drop_na(pred)

save(GAM_predictions, GAM_pm_train, GAM_pm_val, GAM_pm_test,
     file = '../Data/Results/Ridge_regression_w_GAM_correction.RData')

