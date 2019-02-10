

cross_validation_evaluation <- function(x, positives, n_folds = 5, files = '', task = 1, files_similarity = '', n_repeats = 1, alpha = .5, beta= . 5, weight = c(1, 1)) {
	
    ## This function computes the performance of LOTUS in terms of consistency error.
    ## The arguments have to satisfy the following conditions:
    ##  - x is a data frame with each row corresponding to a gene,
    ##  - positives is a set of integers corresponding to the indices of known driver genes,
    ##  - n_folds is an integer determining the number of folds used in the cross-validation procedure,
    ##  - files is a list of strings corresponding to the paths of the kernels that will be used in the SVM classification step,
    ##  - task is an integer corresponding to the index of the disease to consider in the multitask approach (it is 1 in the pan-cancer approach),
    ##  - files_similarity is a list of strings corresponding to the paths of the similarity matrices that will be used in the multitask SVM classication step,
    ##  - n_repeats is an integer corresponding to the number of repeats of the cross-validation procedure,
    ##  - alpha and beta are real numbers in [0,1] corresponding to the weights of the diagonal and constant matrix in the similarity matrix,
    ##  - weight is a vector corresponding to the weights applied to each kernel in files and to the Gram matrix.
  
    require('gplots')
    require('ROCR')
    require('kernlab')
	
    n_genes <- dim(x)[1]
    n_task <- ceiling(max(positives) / n_genes)
    	    
    ## Get specific positives
    
    specific_positives <- positives[which((positives > (n_genes * (task - 1))) & (positives <= (n_genes * task)))]
    n_pos <- length(specific_positives)
        
    ## Make folds (each fold is a subset of positive examples that will be hidden at each iteration)
	
    folds <- list()
    for (i in seq(n_repeats)) {
		    folds <- c(folds, split(sample(seq(n_pos)), rep(1:n_folds, length = n_pos)))
    }
    n_iter <- length(folds)

    ## Cross-validation loop

    pred = list()
    true_labels = list()
    for (iter in seq(n_iter)) {
        cat('-iteration', iter, '-')
        # Move some positive examples to candidates
	cvposset <- setdiff(positives, specific_positives[folds[[iter]]])
        # Make predictions
        pred[[iter]] <- lotus(x = x, positives = cvposset, n_folds = n_folds, files = files, task = task, files_similarity = files_similarity, alpha = alpha, beta = beta, weight = weight)
        
        # Compute labels
        lab <- matrix(0, 1, n_genes)
        lab[specific_positives - n_genes * (task - 1)] <- 1
        lab <- lab[-(intersect(cvposset, specific_positives) - n_genes * (task - 1))]
        true_labels[[iter]] <- lab
        
    }

    # Estimate performance
	
    p <- prediction(pred, true_labels)
    auc <- mean(unlist(performance(p, "auc")@y.values))
    meanrank <- (n_genes - n_pos) * (1 - auc)

    return(meanrank)

}

