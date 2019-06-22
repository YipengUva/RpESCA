#' pESCA model selection based on cross validation error
#'
#' This function implements a missing value based CV model selection
#' approach for the pESCA model on mutliple data sets with different
#' data types, such as quantitative and binary.
#' The details can be found in  \url{https://arxiv.org/abs/1902.06241}.
#'
#' @inheritParams pESCA
#' @param lambdas_g a vector cotains a sequence of values of lambda for
#' the loadings related to quantitative data sets
#' @param lambdas_b a vector cotains a sequence of values of lambda for
#' the loadings related to binary data sets
#'
#' @return This function returns a list contains the results of a pESCA mdoel. \itemize{
#' \item opts_opt: selected lambdas;
#' \item opts_opt: the initialization of the selected model;
#' \item cvErrors_b: the CV errors during the tuning of lambda for
#' loadings related to binary data sets;
#' \item cvErrors_g: the CV errors during the tuning of lambda for
#' loadings related to quantitative data sets;
#' }
#'
#' @examples
#' \dontrun{
#' result_CV <- pESCA_CV_twoSteps(dataSets, dataTypes,
#'                                lambdas_g, lambdas_b,
#'                                penalty='L2', fun_concave='gdp', opts=opts)
#' }
#'
#' @export
pESCA_CV_twoSteps <- function(dataSets, dataTypes, lambdas_g, lambdas_b, penalty = "L2", fun_concave = "gdp", 
    opts = list()) {
    # check if the inputs statisfy the requirements
    stopifnot(class(dataSets) == "list")
    stopifnot(class(penalty) == "character")
    stopifnot(class(fun_concave) == "character")
    if (length(dataTypes) == 1) {
        dataTypes <- unlist(strsplit(dataTypes, split = ""))
    }
    if (exists("quiet", where = opts)) {
        quiet <- opts$quiet
    } else {
        quiet <- 0
    }
    
    # number of data sets, size of each data set
    nTries_g <- length(lambdas_g)
    nTries_b <- length(lambdas_b)
    nDataSets <- length(dataSets)  # number of data sets
    n <- rep(0, nDataSets)  # number of samples
    d <- rep(0, nDataSets)  # numbers of variables in different data sets
    for (i in 1:nDataSets) {
        n[i] <- dim(dataSets[[i]])[1]
        d[i] <- dim(dataSets[[i]])[2]
    }
    if (length(unique(as.factor(n))) != 1) 
        stop("multiple data sets have unequal sample size")
    n <- n[1]
    
    # default dispersion parameters alphas
    if (exists("alphas", where = opts)) {
        alphas <- opts$alphas
    } else {
        alphas <- rep(1, nDataSets)
    }
    
    # create zero matrix to hold results
    cvErrors_g <- matrix(data = 0, nrow = nTries_g, ncol = nDataSets)
    cvErrors_b <- matrix(data = 0, nrow = nTries_b, ncol = nDataSets)
    
    # index out the how many Gaussian data sets and how many nonGuassian data sets
    length_g <- sum(dataTypes == rep("G", nDataSets))
    length_b <- nDataSets - length_g
    
    # split data sets into training set and test set
    splitedData <- dataSplit(dataSets = dataSets, dataTypes = dataTypes, ratio_mis = 0.1)
    trainSets <- splitedData$trainSets
    
    # two steps model selection process
    opts_inner <- opts
    
    # first fix lambda_g0, optimize lambda_b
    lambda_g0 <- lambdas_g[1] * rep(1, length_g)
    inits <- list()  # save the parameters during the model selection
    for (j in 1:nTries_b) {
        lambda_b <- lambdas_b[j] * rep(1, length_b)
        lambdas <- c(lambda_g0, lambda_b)
        trainModel <- pESCA(dataSets = trainSets, dataTypes = dataTypes, lambdas = lambdas, penalty = penalty, 
            fun_concave = fun_concave, opts = opts_inner)
        if ((trainModel$iter <= 2) & (quiet == 0)) {
            print("less than 3 iteration is used.")
        }
        
        # warm start
        mu <- trainModel$mu
        A <- trainModel$A
        B <- trainModel$B
        opts_inner$mu0 <- mu
        opts_inner$A0 <- A
        opts_inner$B0 <- B
        
        # save initilization
        inits[[j]] <- opts_inner
        
        # compute the test error
        ThetaHat <- ones(n) %*% t(mu) + A %*% t(B)
        testError_vec <- cvError_comput(splitedData, dataTypes, alphas, ThetaHat, d)
        cvErrors_b[j, ] <- testError_vec
        
    }
    colnames(cvErrors_b) <- paste0(rep("X_"), as.character(1:nDataSets))
    
    # select the optimal value of lambda_b for nonGuassian data sets
    nonGaussian_cvErrors <- cvErrors_b[, (length_g + 1):nDataSets]
    
    if (length_b > 1) {
        sum_nonGaussian_cvErrors <- rowSums(nonGaussian_cvErrors)
    } else {
        sum_nonGaussian_cvErrors <- nonGaussian_cvErrors
    }
    # remove the first two points, because usually they are not convergenced.  Therefore have relative low
    # CV error
    index_b <- (which.min(sum_nonGaussian_cvErrors[3:nTries_b]) + 2)
    lambda_b_opt <- lambdas_b[index_b] * rep(1, length_b)
    opts_inner <- inits[[index_b]]
    
    # model selection for Gaussian data sets fix lambda_b_opt, optimize lambda_g
    inits <- list()  # save the parameters during the model selection
    for (j in 1:nTries_g) {
        lambda_g <- lambdas_g[j] * rep(1, length_g)
        lambdas <- c(lambda_g, lambda_b_opt)
        trainModel <- pESCA(dataSets = trainSets, dataTypes = dataTypes, lambdas = lambdas, penalty = penalty, 
            fun_concave = fun_concave, opts = opts_inner)
        if ((trainModel$iter <= 2) & (quiet == 0)) {
            print("less than 3 iteration is used.")
        }
        
        # warm start
        mu <- trainModel$mu
        A <- trainModel$A
        B <- trainModel$B
        opts_inner$mu0 <- mu
        opts_inner$A0 <- A
        opts_inner$B0 <- B
        
        # save initilization
        inits[[j]] <- opts_inner
        
        # compute the test error
        ThetaHat <- ones(n) %*% t(mu) + A %*% t(B)
        testError_vec <- cvError_comput(splitedData, dataTypes, alphas, ThetaHat, d)
        cvErrors_g[j, ] <- testError_vec
        
    }
    colnames(cvErrors_g) <- paste0(rep("X_"), as.character(1:nDataSets))
    
    # select the optimal value of lambda_g for Guassian data sets
    Gaussian_cvErrors <- cvErrors_g[, 1:length_g]
    if (length_g > 1) {
        sum_Gaussian_cvErrors <- rowSums(Gaussian_cvErrors)
    } else {
        sum_Gaussian_cvErrors <- Gaussian_cvErrors
    }
    index_g <- which.min(sum_Gaussian_cvErrors)
    lambda_g_opt <- lambdas_g[index_g] * rep(1, length_g)
    opts_opt <- inits[[index_g]]
    
    # save the results
    result_CV <- list()
    result_CV$lambdas_opt <- c(lambda_g_opt, lambda_b_opt)
    result_CV$opts_opt <- opts_opt
    result_CV$cvErrors_g <- cvErrors_g
    result_CV$cvErrors_b <- cvErrors_b
    
    return(result_CV)
}

