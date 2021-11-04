
library(tidyverse)
library(biomaRt)
library(clusterProfiler) ; library(ReactomePA) ; library(DOSE) ; library(org.Hs.eg.db)

################################################################################################################
# LOAD DATA

# SFARI Genes
SFARI_genes = read_csv('./../../SFARI/Data/SFARI_genes_01-03-2020_w_ensembl_IDs.csv')

# Load Gandal dataset
load('./../Data/preprocessedData/preprocessed_data.RData')
datExpr = datExpr %>% data.frame

# WGCNA metrics
WGCNA_metrics = read.csv('./../Data/preprocessedData/WGCNA_metrics.csv')

# Updates genes_info with SFARI information and clusters
genes_info = genes_info %>% left_join(SFARI_genes, by = 'ID') %>% 
  left_join(datGenes %>% mutate(ID = rownames(.)) %>% dplyr::select(ID, hgnc_symbol), by = 'ID') %>%
  dplyr::select(ID, hgnc_symbol, log2FoldChange, shrunken_log2FoldChange, significant, Neuronal) %>%
  left_join(WGCNA_metrics, by = 'ID') %>% dplyr::select(-contains('pval'))

top_modules = c('#44A0FF', '#D177FF', '#F47B5B', '#00BADE', '#64B200', '#DD71FA')


rm(dds, WGCNA_metrics)

################################################################################################################
# Prepare dataset for Enrichment Analysis

EA_dataset = genes_info %>% dplyr::rename('ensembl_gene_id' = ID) %>% filter(Module!='gray')

# ClusterProfile works with Entrez Gene Ids, o we have to assign one to each gene
getinfo = c('ensembl_gene_id','entrezgene')
mart=useMart(biomart='ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl',host='feb2014.archive.ensembl.org')
biomart_output = getBM(attributes=getinfo, filters=c('ensembl_gene_id'), 
                       values=EA_dataset$ensembl_gene_id, mart=mart)

EA_dataset = biomart_output %>% left_join(EA_dataset, by = 'ensembl_gene_id') %>% 
             dplyr::rename('ID' = ensembl_gene_id)

rm(getinfo, mart, biomart_output)

################################################################################################################
# GSEA enrichment

nPerm = 1e5

GSEA_dataset = EA_dataset %>% dplyr::select(ID, entrezgene, contains('MM.'))

GSEA_enrichment = list()

for(module in top_modules){
  
  cat(paste0('\nModule: ', which(top_modules == module), '/', length(top_modules)))
  
  geneList = GSEA_dataset %>% pull(paste0('MM.',substring(module,2)))
  names(geneList) = GSEA_dataset %>% pull(entrezgene) %>% as.character
  geneList = sort(geneList, decreasing = TRUE)
  
  GSEA_GO = gseGO(geneList, OrgDb = org.Hs.eg.db, pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, 
                  nPerm = nPerm, verbose = FALSE, seed = TRUE)
  
  GSEA_DO = gseDO(geneList, pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, 
                  nPerm = nPerm, verbose = FALSE, seed = TRUE)
  
  GSEA_DGN = gseDGN(geneList, pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, 
                    nPerm = nPerm, verbose = FALSE, seed = TRUE)
  
  GSEA_KEGG = gseKEGG(geneList, organism = 'human', pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, 
                      nPerm = nPerm, verbose = FALSE, seed = TRUE)
  
  GSEA_Reactome = gsePathway(geneList, organism = 'human', pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, 
                             nPerm = nPerm, verbose = FALSE, seed = TRUE)
  
  GSEA_enrichment[[module]] = list('GO' = GSEA_GO, 'DO' = GSEA_DO, 'DGN' = GSEA_DGN, 'KEGG' = GSEA_KEGG, 
                                   'Reactome' = GSEA_Reactome)
  
  # Save after each iteration
  save(GSEA_enrichment, file = './../Data/preprocessedData/GSEA_results.RData')
}

rm(GSEA_dataset, nPerm, geneList, GSEA_GO, GSEA_DO, GSEA_DGN, GSEA_KEGG, GSEA_Reactome)

################################################################################################################
# ORA enrichment

# Prepare input
universe = EA_dataset$entrezgene %>% as.character

# Perform Enrichment
ORA_enrichment = list()

for(module in top_modules){
  
  genes_in_module = EA_dataset %>% filter(Module == module) %>% pull(entrezgene)
  
  ORA_GO = enrichGO(gene = genes_in_module, universe = universe, OrgDb = org.Hs.eg.db, ont = 'All', 
                    pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, qvalueCutoff = 1)
  
  ORA_DO = enrichDO(gene = genes_in_module, universe = universe, qvalueCutoff = 1,
                    pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1)
  
  ORA_DGN = enrichDGN(gene = genes_in_module, universe = universe, qvalueCutoff = 1,
                      pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1)
  
  ORA_KEGG = enrichKEGG(gene = genes_in_module, universe = universe, qvalueCutoff = 1,
                        pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1) 
  
  ORA_Reactome = enrichPathway(gene = genes_in_module, universe = universe, qvalueCutoff = 1,
                               pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1)
  
  ORA_enrichment[[module]] = list('GO' = ORA_GO, 'DO' = ORA_DO, 'DGN' = ORA_DGN, 'KEGG' = ORA_KEGG, 
                                  'Reactome' = ORA_Reactome)
  
  # Save after each iteration
  save(ORA_enrichment, file = './../Data/preprocessedData/ORA_results.RData')
}

rm(universe, genes_in_module, module, ORA_GO, ORA_DGN, ORA_DO, ORA_KEGG, ORA_Reactome)

################################################################################################################
# Get shared enrichment for each module

top_modules_enrichment = list()

for(module in top_modules){
  
  module_enrichment = list()
  GSEA_enrichment_for_module = GSEA_enrichment[[module]]
  ORA_enrichment_for_module = ORA_enrichment[[module]]
  
  for(dataset in c('KEGG', 'Reactome', 'GO', 'DO', 'DGN')){
    
    GSEA_enrichment_dataset = GSEA_enrichment_for_module[[dataset]] %>% 
      data.frame %>%
      dplyr::rename('pvalue_GSEA' = pvalue, 
                    'p.adjust_GSEA' = p.adjust, 
                    'qvalues_GSEA' = qvalues)
    
    ORA_enrichment_dataset = ORA_enrichment_for_module[[dataset]] %>% 
      data.frame %>%
      dplyr::rename('pvalue_ORA' = pvalue, 
                    'p.adjust_ORA' = p.adjust, 
                    'qvalue_ORA' = qvalue)
    
    # Get shared enrichments (if any)
    shared_enrichment_dataset = GSEA_enrichment_dataset %>%
      inner_join(ORA_enrichment_dataset, by = 'ID')
    
    module_enrichment[[dataset]] = shared_enrichment_dataset
  }
  
  top_modules_enrichment[[module]] = module_enrichment  
}


save(top_modules_enrichment, file = './../Data/preprocessedData/top_modules_enrichment.RData')

rm(module, module_enrichment, GSEA_enrichment_for_module, ORA_enrichment_for_module, dataset, 
   GSEA_enrichment_dataset, ORA_enrichment_dataset, shared_enrichment_dataset)
