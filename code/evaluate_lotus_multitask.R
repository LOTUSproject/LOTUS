rm(list = ls())

## possible tasks for TSG: all between 1 and 30 except 5, 9 and 23
## possible tasks for OG: all between 1 and 30 except 5, 9 and 14

source('lotus.R')
source('evaluation.R')

## To get reproducible results, we set the seed to the following value.

seed = 123456789
set.seed(seed)

## Experiment parameters:

driver_type <- 'tsg' # 'tsg' or 'og'
task <- 2
version <- 'onetask' # 'lotus' for LOTUS, 'lotus2' for LOTUS2, 'aggregation' for aggregation LOTUS and 'onetask' for onetask LOTUS

## Load driver lists

load(paste(driver_type, '_per_Diseases.RData', sep = '')) 

## Load features

if (version %in% c('lotus', 'lotus2', 'aggregation'))
	x <- read.delim('features_2020.txt', header = T, as.is = T, sep = '\t')
if (version == 'onetask')
	x <- read.delim(paste('features_', names(genes_per_disease)[task], '.txt', sep = ''), header = T, as.is = T, sep = '\t')
x[is.na(x)] <- 0
gene_names <- as.matrix(x$Gene)
n_genes <- length(gene_names)

switch(driver_type,
       'tsg' = {x <- x[, match(c('Frameshift', 'LOF', 'Splice'), colnames(x))]},
       'og' = {x <- x[, match(c('Entropy.Score', 'Missense.Damaging', 'Missense.total'), colnames(x))]})

## Encode drivers in the multitask framework: each cancer type is coded in integers: gene-type pairs for type number 1 will have numbers between 1 and n_genes, for type number 2 between n_genes+1 and 2*n_genes, etc.

ind_positives <- c()
for(i in seq(length(genes_per_disease))){
  ind_positives <- c(ind_positives, n_genes * (i-1) + which(gene_names %in% genes_per_disease[[i]]))
}

## Cross-validation parameters

n_folds <- 2
n_repeats <- 1

## Kernels

files <- 'PPIKernel_2020.RData'

## Weights for the similarity matrix

switch(version,
       'lotus' = {alpha = 0.5
       beta = 0.5},
       'lotus2' = {alpha = 1
       beta = 1},
       'aggregation' = {alpha=0
       beta = 1},
       'onetask' = {alpha = 1
       beta = 0})

files_similarity <- ''
if (version == 'lotus2')
  files_similarity <- 'Kernel_similarity.RData'

weight <- c(1,1)

res <- cross_validation_evaluation(x = x, positives = ind_positives, n_folds = n_folds, files = files, task = task, files_similarity = files_similarity, n_repeats = n_repeats, alpha = alpha, beta = beta, weight = weight)

save(data = res, file = paste('evaluation_', version, '_multitask_', driver_type, '_type=', task, '.RData', sep = ''))


  
  
  
  
