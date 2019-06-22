#' pESCA model selection based on CV error when full information is avaliable
#'
#' This function implements a missing value based CV model selection
#' approach for the pESCA model on mutliple data sets with the same data type.
#' The details can be found in  \url{https://arxiv.org/abs/1902.06241}.
#'
#' @inheritParams pESCA_CV
#' @param simulatedData the output of function \code{dataSimu_group_sparse}
#'
#' @return This function returns a list contains the results of a pESCA mdoel. \itemize{
#' \item cvErrors_mat: a matrix contains the CV errors for the full data set and each
#' single data set;
#' \item RVs_mat: a matrix contains the RV coefficients in estimating the simulated
#' common and distinct structures;
#' \item ranks_mat: a matrix contains the rank estimation of the common and distinct
#' structures;
#' \item RMSEs_mat: a matrix contains the RMSEs in estimating the simulated parameters;
#' \item inits: a list contains the initilizations of all the constructed models;
#' \item outs: a list contains the outputs of all the constructed models;
#' }
#'
#' @examples
#' \dontrun{
#' result_CV_fullInfo <- pESCA_CV_fullInfo(simulatedData,
#'                             lambdas_CV, penalty='L2', fun_concave='gdp', opts=opts)
#' }
#'
#' @export
pESCA_CV_fullInfo <- function(simulatedData, lambdas_CV = NULL, penalty = "L2", fun_concave = "gdp", opts = list()) {
    
    # check if the inputs statisfy the requirements
    stopifnot(class(penalty) == "character")
    stopifnot(class(fun_concave) == "character")
    if (exists("thr_path", where = opts)) {
        thr_path <- opts$thr_path
    } else {
        thr_path <- 0
    }
    
    # index out the data sets and it data types
    dataSets <- simulatedData$X
    dataTypes <- simulatedData$dataTypes
    
    # number of data sets, size of each data set
    nTries <- length(lambdas_CV)
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
    sumd <- sum(d)  # total number of variables
    
    # default dispersion parameters alphas
    if (exists("alphas", where = opts)) {
        alphas <- opts$alphas
    } else {
        alphas <- rep(1, nDataSets)
    }
    if (exists("quiet", where = opts)) {
        quiet <- opts$quiet
    } else {
        quiet <- 0
    }
    
    # create zero matrix to hold results
    cvErrors_mat <- matrix(data = NA, nTries, 4)
    ranks_mat <- matrix(data = NA, nTries, 7)
    RVs_mat <- matrix(data = NA, nTries, 7)
    RMSEs_mat <- matrix(data = NA, nTries, 5)
    
    # model selection process
    opts_inner <- opts
    
    # split data sets into training set and test set
    splitedData <- dataSplit(dataSets = dataSets, dataTypes = dataTypes, ratio_mis = 0.1)
    trainSets <- splitedData$trainSets
    
    # save the parameters during the model selection
    inits <- as.list(1:nTries)
    if (thr_path == 1) {
        outs <- as.list(1:nTries)
    }
    
    # model selection
    for (j in 1:nTries) {
        lambda <- lambdas_CV[j]
        
        # using the training set to construct a ESCA model
        lambdas <- lambda * rep(1, nDataSets)
        
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
        
        inits[[j]] <- opts_inner
        if (thr_path == 1) {
            outs[[j]] <- trainModel$Sigmas
        }
        
        # model evulation with respect to the simulated parameters
        S <- trainModel$S
        eval_metrics <- eval_metrics_simu_group(mu, A, B, S, d, simulatedData)
        
        RVs_mat[j, ] <- eval_metrics$RVs_structs
        ranks_mat[j, ] <- eval_metrics$Ranks_structs
        RMSEs_mat[j, ] <- eval_metrics$RMSEs_params
        
        # compute the test error
        ThetaHat <- ones(n) %*% t(mu) + A %*% t(B)
        testError_vec <- cvError_comput(splitedData, dataTypes, alphas, ThetaHat, d)
        
        cvErrors_tmp <- c(sum(testError_vec), testError_vec)
        cvErrors_mat[j, ] <- cvErrors_tmp
    }
    
    colnames(cvErrors_mat) <- paste0(rep("X_"), c("full", as.character(1:nDataSets)))
    colnames(RVs_mat) <- names(eval_metrics$RVs_structs)
    colnames(ranks_mat) <- names(eval_metrics$Ranks_structs)
    colnames(RMSEs_mat) <- names(eval_metrics$RMSEs_params)
    
    RVs_mat[is.na(RVs_mat)] <- 0
    ranks_mat[is.na(ranks_mat)] <- 0
    
    result_CV <- list()
    result_CV$cvErrors_mat <- cvErrors_mat
    result_CV$RVs_mat <- RVs_mat
    result_CV$ranks_mat <- ranks_mat
    result_CV$RMSEs_mat <- RMSEs_mat
    result_CV$inits <- inits
    if (thr_path == 1) {
        result_CV$outs <- outs
    }
    
    return(result_CV)
}

