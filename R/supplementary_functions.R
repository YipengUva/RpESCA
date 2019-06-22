#' Element-wise logit function
#'
#' This function will do logit transformation on a probability matrix.
#' The formula is \code{logit(x) = log(x/(1-x))}
#'
#' @param x a matrix contains probabilities
#'
#' @return The logit transformation of a probability matrix
#'
#' @examples
#' \dontrun{logit(matrix(rbeta(3*4,1,1),nrow=3,ncol=4))}
#'
#' @export
logit <- function(x) {
    # the domain of logit() should be in (0,1)
    stopifnot(all(x > 0) & all(x < 1))
    
    log(x/(1 - x))
}

#' Element-wise inverse logit function
#'
#' This function will do inveser logit transformation on a matrix
#' contains real numbers. The formula is \code{inver_logit(x) = (1+exp(-x))^(-1)}
#'
#' @param x a matrix contains real numbers
#'
#' @return a probability matrix
#'
#' @examples
#' \dontrun{inver_logit(matrix(rnorm(3*4),nrow=3,ncol=4))}
#'
#' @export
inver_logit <- function(x) {
    1/(1 + exp(-x))
}

#' Index generating function
#'
#' This function will generate a series of indexes used to index out
#' the variables of the ith data set when we only know the number of
#' variables
#'
#' @param i index for the ith data set.
#' @param ds a vector contains the number of variables for all
#' the data sets.
#'
#' @return A series of indexes for the variables in the ith data set
#'
#' @examples
#' \dontrun{index_Xi(2, c(400,200,100))}
index_Xi <- function(i, ds) {
    if (i == 1) {
        columns_Xi <- 1:ds[1]
    } else {
        columns_Xi <- (sum(ds[1:(i - 1)]) + 1):sum(ds[1:i])
    }
    columns_Xi
}

#' Compute the variation expalined raitios when Gaussian data is used
#'
#' This function computes the variation expalined ratios for component
#' model on quantitative data sets. Details can be found
#' in the paper \url{https://arxiv.org/abs/1902.06241}.
#'
#' @param X a \eqn{n*d} quantitative data set
#' @param mu a \eqn{n*1} column offset term
#' @param A a \eqn{n*R} score matrix
#' @param B a \eqn{d*R} loading matrix
#' @param Q a \eqn{n*d} weighting matrix of the same size of X
#'
#' @return This function returns a list contains the variaiton expalined
#' ratios of the whole model (varExp_total) or of each component (varExp_PCs).
#'
#' @examples
#' \dontrun{ out <- varExp_Gaussian(X,mu,A,B,Q)
#'           out$varExp_total
#'           out$varExp_PCs
#' }
varExp_Gaussian <- function(X, mu, A, B, Q) {
    # parameter used
    m = dim(X)[1]
    
    # compute the loglikelihood of mle and null model
    X_centered <- X - ones(m) %*% mu
    QX <- Q * X_centered
    
    # likelihood of the null model
    l_null <- norm(QX, "F")^2  # null model
    
    # likelihood of the full model
    E_hat <- X_centered - tcrossprod(A, B)
    QE_hat <- Q * E_hat
    l_model <- norm(QE_hat, "F")^2  # full model
    
    # compute the least squares of an individual PC
    R <- dim(B)[2]
    l_PCs <- rep(0, R)
    
    for (r in 1:R) {
        Ar <- A[, r]
        Br <- B[, r]
        QPCr <- Q * (tcrossprod(Ar, Br))
        l_PCs[r] <- l_null - 2 * (crossprod((QX %*% Br), Ar)) + crossprod(Ar, (QPCr %*% Br))
    }
    
    # compute variation explained by each PC
    varExp_PCs <- (1 - l_PCs/l_null) * 100
    
    # total variation explained
    varExp_total <- (1 - l_model/l_null) * 100
    
    # return the results
    out <- list()
    out$varExp_total <- varExp_total
    out$varExp_PCs <- varExp_PCs
    return(out)
}

#' ones function
#'
#' This function will generate a column of ones
#'
#' @param n a integer number indicates the dimension of the vector
#'
#' @return a n dimensional column vector contains ones
#'
#' @examples
#' \dontrun{ones(10)}
#'
#' @export
ones <- function(n) {
    stopifnot(n > 0)
    
    as.matrix(rep(1, n))
}

#' Generate log-spaced sequence
#'
#' This function will generate a log-spaced sequence.
#' The inputs are the same as function \code{seq} in base library.
#'
#' @param from The initial value of the sequence
#' @param to The last value of the sequence
#' @param length.out desired length of the sequence
#'
#' @return A \code{length.out} dimensional squence
#'
#' @examples
#' \dontrun{log10_seq(from=1, to=500, length.out=30)}
#'
#' @export
log10_seq <- function(from, to, length.out) {
    # to must larger than from and from must be larger than 0
    stopifnot((to > from) & (from > 0))
    10^(seq(log10(from), log10(to), length.out = length.out))
}

#' Logistic loss function
#'
#' This function will compute the logistic loss when a
#' binary data set and its natural parameter matrix is avaliable.
#'
#' @param X a binary matrix
#' @param Theta a natual parameter (log-odds) matrix
#'
#' @return The logistic loss
#'
#' @examples
#' \dontrun{
#' Theta <- matrix(rnorm(3*4),3,4)
#' X <- matrix(data=0,3,4)
#' X[Theta>0] <- 1
#' obj_logistic(X,Theta)
#' }
obj_logistic <- function(X, Theta) {
    # X: binary data matrix Theta: offset + Z, the log-odds Theta = ones(m,1)*mu + Z; when x=1, the loss
    # is -log(1/(1+exp(-theta))).  when x=0, the loss is -log(1/(1+exp(theta)).  The following is the
    # matrix form
    
    # X must be a binary matrix and contains 1 and 0
    stopifnot(all(unique(as.vector(X)) == c(1, 0)))
    
    # logistic loss
    X <- 2 * X - 1
    tmp <- 1/(1 + exp(-X * Theta))
    out <- -sum(sum(log(tmp)))
    out
}

#' Evluating pESCA model when simulated parameters are avaliable
#'
#' This function will evaluate the the performance of the constructed
#' pESCA model with group penalty when the simulated parameters are
#' avaliable.
#'
#' @param mu estimated offset term in column vector form
#' @param A estimated score matrix
#' @param B estimated loading matrix
#' @param S estimated group sparse pattern on \code{B}
#' @param ds a vector contains the number of variables in multiple data sets
#' @param simulatedData the output of function \code{dataSimu_group_sparse}
#'
#' @return This function returns a list contains \itemize{
#' \item  RVs_structs: a vector contains the RV coefficients in estimating
#' the global common (C123), local common (C12, C13, C23) and distinct structures
#' (D1, D2, D3);
#' \item  Ranks_structs: a vector contains the ranks of estimated
#' C123, C12, C13, C23, D1, D2, D3;
#' \item RMSEs_params: the relative mean squared errors (RMSEs) in estimating
#' the simulated parameters \eqn{\Theta, \Theta_1, \Theta_2, \Theta_3, \mu}
#' }
#'
#' @examples
#' \dontrun{eval_metrics_simu_group(mu,A,B,S,ds,simulatedData)}
#' 
#' @export
eval_metrics_simu_group <- function(mu, A, B, S, ds, simulatedData) {
    Theta_simu <- simulatedData$Theta_simu
    mu_simu <- simulatedData$mu_simu
    U_simu <- simulatedData$U_simu
    D_simu <- simulatedData$D_simu
    V_simu <- simulatedData$V_simu
    
    n <- dim(U_simu)[1]
    
    # simulated parameters Theta1 Theta2 Theta3
    Theta1_simu <- Theta_simu[, index_Xi(1, ds)]
    Theta2_simu <- Theta_simu[, index_Xi(2, ds)]
    Theta3_simu <- Theta_simu[, index_Xi(3, ds)]
    
    # C123
    i <- 1
    index_factors <- (3 * (i - 1) + 1):(3 * i)
    C123_simu <- U_simu[, index_factors] %*% diag(D_simu[index_factors]) %*% t(V_simu[, index_factors])
    
    # C12
    i <- 2
    index_factors <- (3 * (i - 1) + 1):(3 * i)
    C12_simu <- U_simu[, index_factors] %*% diag(D_simu[index_factors]) %*% t(V_simu[, index_factors])
    C12_simu <- C12_simu[, c(index_Xi(1, ds), index_Xi(2, ds))]
    
    # C13
    i <- 3
    index_factors <- (3 * (i - 1) + 1):(3 * i)
    C13_simu <- U_simu[, index_factors] %*% diag(D_simu[index_factors]) %*% t(V_simu[, index_factors])
    C13_simu <- C13_simu[, c(index_Xi(1, ds), index_Xi(3, ds))]
    
    # C23
    i <- 4
    index_factors <- (3 * (i - 1) + 1):(3 * i)
    C23_simu <- U_simu[, index_factors] %*% diag(D_simu[index_factors]) %*% t(V_simu[, index_factors])
    C23_simu <- C23_simu[, c(index_Xi(2, ds), index_Xi(3, ds))]
    
    # D1
    i <- 5
    index_factors <- (3 * (i - 1) + 1):(3 * i)
    D1_simu <- U_simu[, index_factors] %*% diag(D_simu[index_factors]) %*% t(V_simu[, index_factors])
    D1_simu <- D1_simu[, index_Xi(1, ds)]
    
    # D2
    i <- 6
    index_factors <- (3 * (i - 1) + 1):(3 * i)
    D2_simu <- U_simu[, index_factors] %*% diag(D_simu[index_factors]) %*% t(V_simu[, index_factors])
    D2_simu <- D2_simu[, index_Xi(2, ds)]
    
    # D3
    i <- 7
    index_factors <- (3 * (i - 1) + 1):(3 * i)
    D3_simu <- U_simu[, index_factors] %*% diag(D_simu[index_factors]) %*% t(V_simu[, index_factors])
    D3_simu <- D3_simu[, index_Xi(3, ds)]
    
    # model evulation using simulated parameters
    Theta_Hat <- ones(n) %*% t(mu) + A %*% t(B)
    Theta1_Hat <- Theta_Hat[, index_Xi(1, ds)]
    Theta2_Hat <- Theta_Hat[, index_Xi(2, ds)]
    Theta3_Hat <- Theta_Hat[, index_Xi(3, ds)]
    
    C123_index <- (colSums(S) == 3)
    C123_Hat <- A[, C123_index] %*% t(B[, C123_index])
    
    C12_index <- (colSums(S[1:2, ]) == 2) & (!C123_index)
    C12_Hat <- A[, C12_index] %*% t(B[, C12_index])
    C12_Hat <- C12_Hat[, c(index_Xi(1, ds), index_Xi(2, ds))]
    
    C13_index <- (colSums(S[c(1, 3), ]) == 2) & (!C123_index)
    C13_Hat <- A[, C13_index] %*% t(B[, C13_index])
    C13_Hat <- C13_Hat[, c(index_Xi(1, ds), index_Xi(3, ds))]
    
    C23_index <- (colSums(S[c(2, 3), ]) == 2) & (!C123_index)
    C23_Hat <- A[, C23_index] %*% t(B[, C23_index])
    C23_Hat <- C23_Hat[, c(index_Xi(2, ds), index_Xi(3, ds))]
    
    D_index <- (colSums(S) == 1)
    D1_index <- D_index & (S[1, ] == 1)
    D2_index <- D_index & (S[2, ] == 1)
    D3_index <- D_index & (S[3, ] == 1)
    
    D1_Hat <- A[, D1_index] %*% t(B[, D1_index])
    D1_Hat <- D1_Hat[, index_Xi(1, ds)]
    D2_Hat <- A[, D2_index] %*% t(B[, D2_index])
    D2_Hat <- D2_Hat[, index_Xi(2, ds)]
    D3_Hat <- A[, D3_index] %*% t(B[, D3_index])
    D3_Hat <- D3_Hat[, index_Xi(3, ds)]
    
    # RV coefficients of estimated structures
    RV_C123 <- RV_modified(C123_simu, C123_Hat)
    RV_C12 <- RV_modified(C12_simu, C12_Hat)
    RV_C13 <- RV_modified(C13_simu, C13_Hat)
    RV_C23 <- RV_modified(C23_simu, C23_Hat)
    RV_D1 <- RV_modified(D1_Hat, D1_simu)
    RV_D2 <- RV_modified(D2_Hat, D2_simu)
    RV_D3 <- RV_modified(D3_Hat, D3_simu)
    RVs_structures <- c(RV_C123, RV_C12, RV_C13, RV_C23, RV_D1, RV_D2, RV_D3)
    names(RVs_structures) <- c("C123", "C12", "C13", "C23", "D1", "D2", "D3")
    
    # ranks of estimated structures
    Ranks_structures <- c(sum(C123_index), sum(C12_index), sum(C13_index), sum(C23_index), sum(D1_index), 
        sum(D2_index), sum(D3_index))
    names(Ranks_structures) <- c("C123", "C12", "C13", "C23", "D1", "D2", "D3")
    
    # RV coefficients of estimated Theta
    RMSE_Theta1 <- norm(Theta1_Hat - Theta1_simu, "F")^2/norm(Theta1_simu, "F")^2
    RMSE_Theta2 <- norm(Theta2_Hat - Theta2_simu, "F")^2/norm(Theta2_simu, "F")^2
    RMSE_Theta3 <- norm(Theta3_Hat - Theta3_simu, "F")^2/norm(Theta3_simu, "F")^2
    RMSE_Theta <- norm(Theta_Hat - Theta_simu, "F")^2/norm(Theta_simu, "F")^2
    RMSE_mu <- norm(mu - mu_simu, "F")^2/norm(mu_simu, "F")^2
    RMSEs_parameters <- c(RMSE_Theta, RMSE_Theta1, RMSE_Theta2, RMSE_Theta3, RMSE_mu)
    names(RMSEs_parameters) <- c("Theta", "Theta_1", "Theta_2", "Theta_3", "mu")
    
    output <- list()
    output$RVs_structs <- RVs_structures
    output$Ranks_structs <- Ranks_structures
    output$RMSEs_params <- RMSEs_parameters
    return(output)
}


#' Modified RV coefficient of two matrices
#'
#' This function will compute the modified RV coefficient of two
#' matrices. The details of the modified RV coefficient can be found
#' in the paper \url{https://academic.oup.com/bioinformatics/article/25/3/401/244239}.
#'
#' @param X a matrix
#' @param Y another matrix
#'
#' @return This function returns the modified RV coefficient between
#' two matrices
#'
#' @examples
#' \dontrun{RV_modified(X,Y)}
RV_modified <- function(X, Y) {
    # RV modifed coefficient by bda group
    AA <- X %*% t(X)
    BB <- Y %*% t(Y)
    AA0 <- AA - diag(diag(AA))
    BB0 <- BB - diag(diag(BB))
    RV <- sum(diag(AA0 %*% BB0))/norm(AA0, "F")/norm(BB0, "F")
    return(RV)
}

#' Split multiple data sets into training and test sets
#'
#' This function will split multiple data sets into training and test
#' sets. Nonmissing elements are randomly selected as the test sets.
#' Then the selected elements are taken as missing, and regarded as
#' training sets. The details can be found in \url{https://arxiv.org/abs/1902.06241}.
#'
#' @inheritParams pESCA_CV
#' @param ratio_mis how many percent of test set could be? default: 0.1
#'
#' @return This function returns a list contains \itemize{
#' \item trainSets: a list contains the training sets;
#' \item testSets: a list contains the test sets;
#' \item indexSets: a list contains the index sets.
#' }
#'
#' @examples
#' \dontrun{dataSplit(dataSets,dataTypes,ratio_mis=0.1)}
dataSplit <- function(dataSets, dataTypes, ratio_mis = 0.1) {
    # number of data sets, size of each data set
    nDataSets <- length(dataSets)  # number of data sets
    n <- rep(0, nDataSets)  # number of samples
    d <- rep(0, nDataSets)  # numbers of variables in different data sets
    for (i in 1:nDataSets) {
        n[i] <- dim(dataSets[[i]])[1]
        d[i] <- dim(dataSets[[i]])[2]
    }
    n <- n[1]
    
    # split data sets into training set and test set
    trainSets <- as.list(1:nDataSets)  # training set
    testSets <- as.list(1:nDataSets)  # test set
    indexSets <- as.list(1:nDataSets)  # index of the test set
    for (i in 1:nDataSets) {
        # index out the i-th data set
        Xi <- dataSets[[i]]
        dataType_Xi <- dataTypes[i]
        
        # generate the index of the test set
        full_ind_vec <- 1:(n * d[i])
        
        # if it is binary data, using hierachical sampling
        if (dataType_Xi == "B") {
            ones_ind_vec <- full_ind_vec[Xi == 1]
            zeros_ind_vec <- full_ind_vec[Xi == 0]
            index_Xi_ones <- sample(ones_ind_vec, round(ratio_mis * length(ones_ind_vec)))
            index_Xi_zeros <- sample(zeros_ind_vec, round(ratio_mis * length(zeros_ind_vec)))
            
            # test the sampled samples
            if (!(all(Xi[index_Xi_ones] == 1)) | !(all(Xi[index_Xi_zeros] == 0))) 
                message("the hierachical sampling does not work")
            
            index_Xi_test <- c(index_Xi_ones, index_Xi_zeros)
        } else {
            non_NaN_mat <- 1 - is.na(Xi)
            non_NaN_ind_vec <- full_ind_vec[non_NaN_mat > 0]
            index_Xi_test <- sample(non_NaN_ind_vec, round(ratio_mis * length(non_NaN_ind_vec)))
        }
        
        # generate the train set
        Xi_train <- Xi
        Xi_train[index_Xi_test] <- NA
        trainSets[[i]] <- Xi_train
        
        # generate the test set
        Xi_test <- Xi[index_Xi_test]
        testSets[[i]] <- Xi_test
        indexSets[[i]] <- index_Xi_test
    }
    
    # return
    result <- list()
    result$trainSets <- trainSets
    result$testSets <- testSets
    result$indexSets <- indexSets
    return(result)
}

#' Compute CV errors
#'
#' This function will compute CV errors for a specific model
#'
#' @param splitedData output of function \code{dataSplit}
#' @param dataTypes the data types for each data set
#' @param alphas dispersion parameters for each data set
#' @param ThetaHat estimated Theta
#' @param d a numeric vector contains the number of variables of data sets
#'
#' @return This function returns a vector contains CV errors
#'
#' @examples
#' \dontrun{cvError_comput(splitedData,dataTypes,alphas,ThetaHat,d)}
cvError_comput <- function(splitedData, dataTypes, alphas, ThetaHat, d) {
    
    nDataSets <- length(d)
    testSets <- splitedData$testSets
    indexSets <- splitedData$indexSets
    
    testError_vec <- rep(0, nDataSets)
    
    for (i in 1:nDataSets) {
        # index out ThetaHat_Xi
        columns_Xi <- index_Xi(i, d)
        ThetaHat_Xi <- ThetaHat[, columns_Xi]
        
        # compute the CV error
        index_Xi_test <- indexSets[[i]]
        Xi_test <- testSets[[i]]
        dataType_Xi <- dataTypes[i]
        if (dataType_Xi == "G") {
            testError_Xi <- (1/alphas[i]) * 0.5 * norm(Xi_test - ThetaHat_Xi[index_Xi_test], "2")^2
        } else if (dataType_Xi == "B") {
            testError_Xi <- (1/alphas[i]) * obj_logistic(Xi_test, ThetaHat_Xi[index_Xi_test])
        }
        testError_vec[i] <- testError_Xi
    }
    
    # return
    return(testError_vec)
}

#' A function to compute the trace of two matrices product
#'
#' This function will compute the trace of two matrices
#'
#' @param X a numerical matrix
#' @param Y a numerical matrix with the same size as \code{X}
#'
#' @return This function returns a scalar contains the trace
#'
#' @examples
#' \dontrun{trace_fast(X, Y)}
trace_fast <- function(X, Y) {
    # fast trace function if n>p, trace(X,Y) = trace(X'Y); if n<p, trace(X,Y) = trace(YX');
    
    stopifnot(all(dim(X) == dim(Y)))
    
    result <- fast_traceC(X, Y)
    return(result)
}
