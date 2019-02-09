
lotus <- function(x,positives,n_folds=2,files='',task=1,files_similarity='',alpha=.5,beta=.5,weight=c(1,1)) {

    ## This function implements LOTUS.
    ## The arguments have to satisfy the following conditions:
    ##  - x is a data frame with each row corresponding to a gene,
    ##  - positives is a set of integers corresponding to the indices of known driver genes,
    ##  - n_folds is an integer determining the number of folds used in the cross-validation procedure,
    ##  - files is a list of strings corresponding to the paths of the kernels that will be used in the SVM classification step,
    ##  - task is an integer corresponding to the index of the disease to consider in the multitask approach (it is 1 in the pan-cancer approach),
    ##  - files_similarity is a list of strings corresponding to the paths of the similarity matrices that will be used in the multitask SVM classication step,
    ##  - alpha and beta are real numbers in [0,1] corresponding to the weights of the diagonal and the constant matrix in the definition of the similarity matrix,
    ##  - weight is a vector corresponding to the weights applied to each kernel in files and to the Gram matrix.
    ## The function outputs a list with two fields: 
    ##  - 'scores' are the scores for every gene in the complementary of positives,
    ##  - 'performance' summarizes the performance of the classication method for every value of parameter C in a geometric grid between sqrt(2)^-5 and sqrt(2)^5 (in row) and for every value of the number of used features between 3 and the total number (in column).
    
    require('gplots')
    require('ROCR')
    require('kernlab')
    
    ## Compute similarity matrix
    
    n_genes <- dim(x)[1]
    n_task <- max(ceiling(positives/n_genes))

    K_similarity <- similarity_matrix(n_task=n_task,files_similarity=files_similarity,alpha=alpha,beta=beta)
	
    ## Classification

    cost_list <- sqrt(2)^seq(-5,5)
    performance <- rep(0,length(cost_list))
    for (i_loop in seq(length(cost_list))) {
        performance[i_loop] <- CV_performance(x=x,positives=positives,cost=cost_list[i_loop],n_folds=n_folds,files=files,task=task,K_similarity=K_similarity)
    }

    optimal_parameter <- which(performance==min(performance),arr.ind=T)
    optimal_cost <- cost_list[optimal_parameter]
    scores <- lotus_svm(x=x,positives=positives,cost=optimal_cost,files=files,task=task,K_similarity=K_similarity,weight=weight)
    return(scores)

}

lotus_svm <- function(x,positives,cost,files,K_similarity,task,weight) {

    ## This function implements the tuned SVM classification method of LOTUS.
    ## The arguments have to satisfy the following conditions:
    ##  - x is a matrix or a data frame with each gene corresponding to a gene,
    ##  - positives is a set of integers corresponding to the indices of known driver genes,
    ##  - cost is a real number used as C-parameter by SVM,
    ##  - files is a list of strings corresponding to the paths of the kernels that will be used in the SVM classification step,
    ##  - K_similarity is a matrix used as similarity matrix in the multitask SVM approach,
    ##  - task is an integer corresponding to the index of the disease to consider in the multitask approach (it is 1 in the pan-cancer approach),
    ##  - weight is a vector of weights corresponding to the weights applied to each kernel in files and to the Gram matrix.
    ## The function outputs a list of scores for all genes in the complementary of positives.
    
    n_genes <- dim(x)[1]
    n_zero <- 40
    n_task <- dim(K_similarity)[1]
	
    ## Add negative examples
    
    x_zero <- x
    x_zero[(n_genes+1):(n_genes+n_zero),] <- matrix(0,ncol=length(colnames(x)),nrow=n_zero)

    negatives <- c()
    for (i in seq(n_task)) {
        negatives <- c(negatives,(i-1)*(n_genes+n_zero)+seq(n_genes+1,n_genes+n_zero))
    }

    ## Take into account the addition of negatives in the indexation of the positives
    
    for (i in seq(length(positives))) {
        task_i <- ceiling(positives[i]/n_genes)
        index_i <- positives[i] - n_genes*(task_i-1)
        positives[i] <- index_i + (n_genes+n_zero)*(task_i-1)
    }
    candidates <- setdiff(seq((n_genes+n_zero)*(task-1)+1,(n_genes+n_zero)*(task-1)+n_genes),positives)
    
    ## Compute kernels
   
    Kernel <- Kernel_matrix(x_zero,files,n_zero,weight=weight)
    
    k <- function(x,y) {
        task_x <- ceiling(x/(n_genes+n_zero))
        task_y <- ceiling(y/(n_genes+n_zero))
        index_x <- x - (n_genes+n_zero)*(task_x-1)
        index_y <- y - (n_genes+n_zero)*(task_y-1)
        return(Kernel[index_x,index_y]*K_similarity[task_x,task_y])
    }
    class(k) <- "kernel"

    ## Compute classifier

    Gene <- seq((n_genes+n_zero)*n_task)
    label <- matrix(0,nrow=((n_genes+n_zero)*n_task),ncol=1)
    x <- data.frame(Gene,label)
    x$label[positives] <- 1
    formula <- as.formula("label ~ Gene")
    m <- ksvm(formula,data=x[c(positives,negatives),],type="C-svc",kernel=k,C=cost,prob.model=FALSE,scaled=FALSE,shrinking=FALSE)

    ## Predict scores
			
    predict(m,x[candidates,],type="decision")

}

Kernel_matrix <- function(x,files,n_zero,weight=c(1,0)) {

    ## This function computes the specific kernel for lotus_svm.
    ## The arguments have to satisfy the following conditions:
    ##  - x is a matrix or a data frame with each gene corresponding to a gene,
    ##  - files is a list of strings corresponding to the paths of the kernels that will added to the linear kernel computed from the feature matrix x, these kernels must be of size dim(x)[1]-n_zero,
    ##  - n_zero is the number of created negative genes in the method (i.e. the weight of the "zero" gene),
    ##  - weights is a vector corresponding to the weights applied to the kernels in files and to the Gram matrix.
    ## The function outputs a square matrix of size dim(x)[1].
    
    n_genes <- dim(x)[1]-n_zero

    ## Compute the linear kernel

    support <- c()
    for (i in seq(dim(x)[2])) {
      if (any(x[, i] != 0)) {
        support <- c(support, i)
      }        
    }
    x <- as.matrix(x[, support])
    
    Kernel <- weight[1]*tcrossprod(as.matrix(scale(x)))/(dim(x)[1]*dim(x)[2])
    
    ## Extract kernels from files.
	
    n_files <- length(files)

    if (files!='') {
        Kernels <- list()
        for (i in seq(n_files)) {
            Kernels[[i]] <- diag(n_genes+n_zero)
            Kernels[[i]][1:n_genes,1:n_genes] <- loadRData(files[i])[1:n_genes,1:n_genes]
        }

        ## Add all kernels

        for (i in seq(n_files)) {
            Kernel <- Kernel + weight[2]*Kernels[[i]]
        }
        Kernel <- Kernel/(n_files+1)
    }

    return(Kernel)

}

similarity_matrix <- function(n_task,files_similarity,alpha=.5,beta=.5) {
    
    ## This function computes the specific similarity matrix for lotus_svm in the multitask approach.
    ## The arguments have to satisfy the following conditions:
    ##  - n_task is an integer corresponding to the number of considered diseases in the multitask approach,
    ##  - files_similarity is a list of strings corresponding to the paths of the similarity matrices that will added to a classical similarity matrix, these matrices must be of size n_task.
    ## The function outputs a square matrix of size n_task.
  
    ## Compute a classical similarity matrix

    K_similarity <- alpha*diag(n_task)+beta*matrix(1,n_task,n_task)
  
    ## Extract similarity matrices from files.
  
    n_files <- length(files_similarity)
  
    if (files_similarity!='') {
        M <- list()
        for (i in seq(n_files)) {
            M[[i]] <- loadRData(files_similarity[i])
        }
    
    ## Add all matrices
    
    for (i in seq(n_files)) {
        K_similarity <- K_similarity + M[[i]]
    }
    K_similarity <- K_similarity/(n_files+1)
    }
  
    return(K_similarity)
    
}

loadRData <- function(fileName){
	
    ## This function loads an .RData file into the desired variable name.

    load(fileName)
    get(ls()[ls() != "fileName"])

}

CV_performance <- function(x,positives,cost,n_folds,files,K_similarity,task) {
    
    ## This function evaluates lotus_svm so as to optimize parameter C.
    ## The arguments have to satisfy the following conditions:
    ##  - x is a matrix or a data frame with each gene corresponding to a gene,
    ##  - positives is a set of integers corresponding to the indices of known driver genes,
    ##  - cost is a real number corresponding to the value of the C-parameter,
    ##  - n_folds is an integer determining the number of folds to be used in the cross-validation procedure,
    ##  - files is a list of strings corresponding to the paths of the kernels that will added to the linear kernel computed from the feature matrix x, these kernels must be of size dim(x)[1],
    ##  - K_similarity is a matrix used as similarity matrix in the multitask SVM approach,
    ##  - task is an integer corresponding to the index of the disease to consider in the multitask approach (it is 1 in the pan-cancer approach).
    ## The function outputs a real number that is the score of the method with the specified value of C.
        
    n_genes <- dim(x)[1]

    ## Positives for tested type
    
    specific_positives <- positives[which((positives>(n_genes*(task-1)))&(positives<=(n_genes*task)))]
    n_specific_positives <- length(specific_positives)
    
    ## Make folds (each fold is a subset of positive examples that will be hidden at each iteration)

    folds <- split(sample(seq(n_specific_positives)),rep(1:n_folds,length=n_specific_positives))
    n_iter <- length(folds)
        
    ## Cross-validation loop

    pred = list()
    true_labels <- list()
    for (iter in seq(n_iter)) {
        cat('.')
        # Move some positive examples to candidates
        cvposset <- setdiff(positives,specific_positives[folds[[iter]]])
        # Make predictions
        pred[[iter]] <- lotus_svm(x=x,positives=cvposset,cost=cost,files=files,K_similarity=K_similarity,task=task,weight=weight)
        # Compute labels
        lab <- matrix(0,1,n_genes)
        lab[specific_positives-n_genes*(task-1)] <- 1
        lab <- lab[-(intersect(cvposset,specific_positives)-n_genes*(task-1))]
        true_labels[[iter]] <- lab
        
    }

    ## Estimate performance

    p <- prediction(pred,true_labels)
    auc <- mean(unlist(performance(p,"auc")@y.values))
    meanrank <- (n_genes-n_specific_positives)*(1-auc)

    return(meanrank)

}







