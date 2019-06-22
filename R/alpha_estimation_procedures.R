#' alpha estimation procedure
#'
#' This function will estimate the noise level (alpha is the variation parameter)
#' of a quantitative data set using the PCA model. The rank of the PCA model
#' is selected using a missing value based cross validation procedure. The
#' details of this function can be found in \url{https://arxiv.org/abs/1902.06241}.
#'
#' @param X a quantitative data set
#' @param K K-fold cross validation
#' @param Rs the searching range for the number of components
#' @param opts a list contains the setting for the algorithm. \itemize{
#' \item tol_obj: tolerance for relative change of hist_obj, default: 1E-6;
#' \item maxit: max number of iterations, default: 1000;
#' }
#'
#' @return This function returns a list contains \itemize{
#' \item alphas_mean: the mean of the K estimated alphas
#' \item alphas_std: the std of the K estimated alphas
#' \item R_CV: the number of PCs with minimum cross validation error
#' \item cvErrors: K-fold cross validation errors
#' }
#'
#' @import RSpectra
#' @importFrom stats sd
#'
#' @examples
#' \dontrun{alpha_estimation(X,K=3,Rs = 1:15,opts=list())}
#'
#' @export
alpha_estimation <- function(X, K = 3, Rs = 1:15, opts = list()) {
    # check if the whole row or whole column is missing
    # index out the rows and columns not fully missing
    W <- 1 - is.na(X)
    X <- X[rowSums(W) > 0, colSums(W) > 0]
    
    # first center the data sets
    m <- dim(X)[1]
    n <- dim(X)[2]
    X <- scale(X, center = TRUE, scale = FALSE)  # NAs are omitted
    
    # parameters
    W <- 1 - is.na(X)
    mn_nonNaN <- sum(W)
    
    # model selection
    md_svd_CV <- svd_CV(X, K, Rs, opts)
    cvErrors <- md_svd_CV$cvErrors
    
    # select the number of components
    R_CV <- Rs[apply(cvErrors, 2, which.min)]
    
    # estimate the noise level
    alphas_CV <- rep(NA, length(R_CV))
    
    # first do a SVD on X to accelerate computation
    R_CV_max <- max(R_CV)
    if (all(!is.na(X))) {
        svd_tmp <- RSpectra::svds(X, R_CV_max, nu = R_CV_max, nv = R_CV_max)
        U <- svd_tmp$u
        S <- diag(svd_tmp$d)
        V <- svd_tmp$v
    } else {
        svd_tmp <- svd_mis(X, R_CV_max, opts)
        U <- svd_tmp$U
        S <- svd_tmp$S
        V <- svd_tmp$V
    }
    
    for (i in 1:length(R_CV)) {
        R_CV_tmp <- R_CV[i]
        Z_hat <- U[, 1:R_CV_tmp] %*% S[1:R_CV_tmp, 1:R_CV_tmp] %*% t(V[, 1:R_CV_tmp])
        DF <- mn_nonNaN - (m + n) * R_CV_tmp
        X_nonNaN <- X
        X_nonNaN[is.na(X)] <- 0
        Z_hat[is.na(X)] <- 0
        sigmaSqure <- (1/DF) * norm(X_nonNaN - Z_hat, "F")^2
        alphas_CV[i] <- sigmaSqure
    }
    alphas_mean <- mean(alphas_CV)
    alphas_std <- stats::sd(alphas_CV)
    
    result <- list()
    result$alphas_mean <- alphas_mean
    result$alphas_std <- alphas_std
    result$alphas_CV <- alphas_CV
    result$R_CV <- R_CV
    result$cvErrors <- cvErrors
    return(result)
}

#' Model selection of a svd model using missing value based CV error
#'
#' This function implemented a missing value based CV model selection approach.
#' First, ratio_mis percent elements are randomly selected as missing values. After
#' that a EM-SVD model is constructed to estimate the prediction error.
#' The details of this function can be found in \url{https://arxiv.org/abs/1902.06241}.
#'
#' @inheritParams alpha_estimation
#' @param ratio_mis the propotion of missing values
#'
#' @return This function returns a list contains \itemize{
#' \item Rs: the number of PCs used for model selection
#' \item cvErrors: K-fold cross validation errors
#' }
#'
#' @examples
#' \dontrun{svd_CV(X,K=3,Rs = 1:15,opts=list(),ratio_mis=0.1)}
#'
#' @export
svd_CV <- function(X, K = 3, Rs = 1:15, opts = list(), ratio_mis = 0.1) {
    # structure to hold results
    length_Rs <- length(Rs)
    cvErrors <- matrix(data = NA, nrow = length_Rs, ncol = K)
    
    # number of non-missing elements
    m <- dim(X)[1]
    n <- dim(X)[2]
    mn <- m * n
    W <- 1 - is.na(X)
    mn_nonNaN <- sum(W)
    
    for (k in 1:K) {
        # seperate the Xtest and Xtrain taken into account the potential problem of NaN
        full_ind_vec <- 1:mn
        non_NaN_ind_vec <- full_ind_vec[W > 0]
        index_X_test <- sample(non_NaN_ind_vec, round(ratio_mis * length(non_NaN_ind_vec)))
        X_train <- X
        X_train[index_X_test] <- NA
        X_test <- X[index_X_test]
        
        # use opts_inner to do warm start
        opts_inner <- opts
        
        # for loop
        for (j in length_Rs:1) {
            R <- Rs[j]
            
            # using the remaining data to construct a SVD model
            trainModel <- svd_mis(X_train, R, opts_inner)
            U <- trainModel$U
            S <- trainModel$S
            V <- trainModel$V
            
            if (R == 1) {
                ZHat <- S * (U %*% t(V))
            } else {
                ZHat <- U %*% S %*% t(V)
            }
            
            # warm start
            opts_inner$U0 = U
            opts_inner$S0 = S
            opts_inner$V0 = V
            
            # extract the estimated parameters for the prediction of missing elements
            X_pred <- ZHat[index_X_test]
            
            # compute the prediction error
            cvErrors[j, k] <- 0.5 * norm(X_test - X_pred, "2")^2
        }
    }
    colnames(cvErrors) <- paste(1:K, "th cv")
    
    result <- list()
    result$cvErrors <- cvErrors
    result$Rs <- Rs
    return(result)
}

#' A svd algorithm with the option for missing values
#'
#' This function implemented a MM algorithm to fit a svd on a quantitaive data
#' set with missing values. The details of this function can be found
#' in \url{https://arxiv.org/abs/1902.06241}.
#'
#' @param X a quantitative data set
#' @param R the number of PCs
#' @param opts a list contains the setting for the algorithm. \itemize{
#' \item tol_obj: tolerance for relative change of hist_obj, default:1E-6;
#' \item maxit: max number of iterations, default: 1000;
#' }
#'
#' @return This function returns a list contains \itemize{
#' \item U: the left singular vectors
#' \item S: the diagnal matrix contains the singular values
#' \item V: the right singular vectors
#' \item iter: the number of iterations used
#' \item cvErrors: K-fold cross validation errors
#' \item diagnose: records hist_obj and rel_obj
#' }
#'
#' @import RSpectra
#'
#' @examples
#' \dontrun{svd_mis(X,R = 3,opts=list())}
#'
#' @export
svd_mis <- function(X, R, opts = list()) {
    # default parameters
    if (exists("tol_obj", where = opts)) {
        tol_obj <- opts$tol_obj
    } else {
        tol_obj <- 1e-06
    }
    if (exists("maxit", where = opts)) {
        maxit <- opts$maxit
    } else {
        maxit <- 1000
    }
    if (exists("quiet", where = opts)) {
        quiet <- opts$quiet
    } else {
        quiet <- 0
    }
    
    # form weighting matrix
    Wc <- is.na(X) - 0
    W <- 1 - Wc  # weighting matrix
    X[is.na(X)] <- 0  # remove missing elements
    
    # initialization initial parameters
    if (exists("U0", where = opts)) {
        U0 <- opts$U0[, 1:R]
        S0 <- opts$S0[1:R, 1:R]
        V0 <- opts$V0[, 1:R]
    } else {
        svd_tmp <- RSpectra::svds(X, R, nu = R, nv = R)
        U0 <- svd_tmp$u
        S0 <- diag(svd_tmp$d)
        if (length(svd_tmp$d) == 1) {
            S0 <- svd_tmp$d
        }
        V0 <- svd_tmp$v
    }
    if (R != 1) {
        Z0 <- U0 %*% tcrossprod(S0, V0)
    } else {
        Z0 <- S0 * tcrossprod(U0, V0)
    }
    
    # specify structure to hold the diagnose results
    diagnose <- list()
    diagnose$hist_objs <- vector()
    diagnose$rel_objs <- vector()
    
    # initial value of loss function
    obj0 <- 0.5 * norm(W * (X - Z0), "F")^2
    diagnose$hist_objs[1] <- obj0
    
    # iterations
    for (k in 1:maxit) {
        if (quiet == 0) 
            print(paste(k, "th iteration"))
        
        # form Xtilde
        Xtilde <- W * X + Wc * Z0
        
        # update Z
        svd_tmp <- RSpectra::svds(Xtilde, R, nu = R, nv = R)
        U <- svd_tmp$u
        V <- svd_tmp$v
        
        if (R != 1) {
            S <- diag(svd_tmp$d)
            Z <- U %*% (tcrossprod(S, V))
        } else {
            S <- svd_tmp$d
            Z <- S * (tcrossprod(U, V))
        }
        
        # new objective value
        obj <- 0.5 * norm(W * (X - Z), "F")^2
        
        # reporting
        diagnose$hist_objs[k + 1] <- obj
        diagnose$rel_objs[k] <- (obj0 - obj)/(obj0 + 1)
        
        # stopping checks
        if (diagnose$rel_objs[k] < tol_obj) 
            break
        
        # save previous results
        Z0 <- Z
        obj0 <- obj
    }
    result <- list()
    result$U <- U
    result$S <- S
    result$V <- V
    result$iter <- k
    result$diagnose <- diagnose
    return(result)
}

