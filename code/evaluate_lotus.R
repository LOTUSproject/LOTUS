
rm(list=ls())

source('lotus.R')
source('evaluation.R')

## To get reproducible results, we set the seed to the following value.

seed = 123456789
set.seed(seed)

## Experiment parameters: Are you looking for TSG or OG? Do you want to use the datasets (mutations and known drivers) from 20/20, from MutSigCV, from CGCv86 or from TUSON? Do you want to include PPI information in addition to the Gram matrix, or not, or do you want to use the PPI kernel without the Gram matrix?

driver_type <- 'og'     # 'tsg' or 'og'
dataset <- 'cosmicv86'        # '2020', 'mutsig', 'cosmicv86' or 'tuson'
ppi <- 'yes'           # 'yes', 'no' or 'only'

## Load features

x <- read.delim(paste('features_', dataset, '.txt', sep = ''), as.is = T, sep = '\t',  header = T)
x[is.na(x)] <- 0
gene_names <- as.matrix(x$Gene)

switch(driver_type,
	'tsg'={x <- x[, match(c('Frameshift', 'LOF', 'Splice'), colnames(x))]},
	'og'={x <- x[, match(c('Entropy.Score', 'Missense.Damaging', 'Missense.total'), colnames(x))]})

## Load driver gene data

positives <- read.table(paste(driver_type, '_', dataset, '.txt', sep = ''))$V1
ind_positives <- which(gene_names %in% positives)

## Cross-validation parameters

n_folds <- 5
n_repeats <- 2

## The following parameter weight sets the weights of the SVM kernel: it equals alpha * Gram_matrix_of_the_features + beta * PPI_Kernel, and weight = [alpha, beta].
## Hence, set weight = c(1, 1) for LOTUS, weight = c(1, 0) for LOTUS with no PPI information, and weight = c(0, 1) for LOTUS with only PPI_information. 

switch(ppi,
	'yes' = {weight <- c(1, 1)},
	'no' = {weight <- c(1, 0)},
	'only' = {weight <- c(0, 1)})

## Kernel

files <- paste('PPIKernel_', dataset, '.RData', sep = '')

## Consistency error

consistency_error <- cross_validation_evaluation(x = x, positives = ind_positives, n_folds = n_folds, files = files, n_repeats = n_repeats, weight = weight)

switch(ppi,
	'yes' = {complement = ''},
	'no' = {complement = '_noPPI'},
	'only'={complement='_noGram'})
write.table(consistency_error, file = paste('evaluate_lotus_', driver_type, '_', dataset, complement, '.txt', sep = ''), row.names = F, col.names = F, quote = F)


