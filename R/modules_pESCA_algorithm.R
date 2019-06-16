#' Updating loading matrix B with conave L2 norm penalty
#'
#' This is an intermediate step of the algorithm for fitting pESCA model. The
#' details of this function can be found in ref thesis.
#'
#' @param JHk An output of the majorizaiotn step
#' @param A The score matrix A during k-th iteration
#' @param B0 The loading matrix B during the previous iteration
#' @param Sigmas0 The group length during the previous iteration
#' @param d A numeric vector contains the numbers of variables in different data sets
#' @param fun_concave A string indicates the used concave function
#' @param alphas The dispersion parameters of exponential dispersion families
#' @param rhos An output of the majorizaiotn step
#' @param lambdas A numeric vector indicates the values of tuning parameters for
#' each data set.
#' @param gamma The hyper-parameter of the concave penalty
#'
#' @return This function returns the updated loading matrix B.
#'
#' @examples
#' \dontrun{
#' B <- update_B_L2(JHk,A,B0,Sigmas0,d,
#'                    fun_concave,alphas,rhos,lambdas,gamma)
#' }
update_B_L2 <- function(JHk, A, B0, Sigmas0, d, fun_concave, alphas, rhos, lambdas, gamma) {
    sumd <- sum(d)
    nDataSets <- length(d)
    n <- dim(A)[1]
    R <- dim(A)[2]
    hfun_sg <- get(paste0(fun_concave, "_sg"))  # super gradient of penalty function
    
    B <- matrix(NA, sumd, R)
    for (i in 1:nDataSets) {
        columns_Xi <- index_Xi(i, d)
        JHk_i <- JHk[, columns_Xi]
        JHkitA <- crossprod(JHk_i, A)
        alpha_i <- alphas[i]
        rho_i <- rhos[i]
        weight_i <- sqrt(d[i])  # weight for L2 norm
        lambda_i <- lambdas[i] * weight_i * alpha_i/rho_i
        
        for (r in 1:R) {
            # form weights of the penalty according to previous sigma0_ir
            sigma0_ir <- Sigmas0[i, r]
            omega_ir <- hfun_sg(sigma0_ir, gamma = gamma, lambda = 1)  # weights
            
            # proximal operator of L2 norm
            JHkitA_r <- JHkitA[, r]
            lambda_ir <- lambda_i * omega_ir
            JHkitA_r_norm <- norm(JHkitA_r, "2")
            
            B_ir <- max(0, 1 - (lambda_ir/JHkitA_r_norm)) * JHkitA_r
            B[columns_Xi, r] <- B_ir
        }
    }
    
    return(B)
}

#' Updating loading matrix B with conave L1 norm penalty
#'
#' This is an intermediate step of the algorithm for fitting pESCA model. The
#' details of this function can be found in ref thesis.
#'
#' @inheritParams update_B_L2
#'
#' @return This function returns the updated loading matrix B.
#'
#' @examples
#' \dontrun{
#' B <- update_B_L1(JHk,A,B0,Sigmas0,d,
#'                    fun_concave,alphas,rhos,lambdas,gamma)
#' }
update_B_L1 <- function(JHk, A, B0, Sigmas0, d, fun_concave, alphas, rhos, lambdas, gamma) {
    sumd <- sum(d)
    nDataSets <- length(d)
    n <- dim(A)[1]
    R <- dim(A)[2]
    hfun_sg <- get(paste0(fun_concave, "_sg"))  # super gradient of the concave function
    
    B <- matrix(NA, sumd, R)
    for (i in 1:nDataSets) {
        columns_Xi <- index_Xi(i, d)
        JHk_i <- JHk[, columns_Xi]
        JHkitA <- crossprod(JHk_i, A)
        alpha_i <- alphas[i]
        rho_i <- rhos[i]
        weight_i <- d[i]  # weight for L1 norm
        lambda_i <- lambdas[i] * weight_i * alpha_i/rho_i
        
        for (r in 1:R) {
            # form weights of the penalty according to previous sigma0_ir
            sigma0_ir <- Sigmas0[i, r]
            omega_ir <- hfun_sg(sigma0_ir, gamma = gamma, lambda = 1)  # weights
            
            # proximal operator of L1 norm
            JHkitA_r <- JHkitA[, r]
            lambda_ir <- lambda_i * omega_ir
            JHkitA_r <- JHkitA[, r]
            
            B_ir_pre <- abs(JHkitA_r) - lambda_ir
            B_ir_pre[B_ir_pre < 0] <- 0
            B_ir <- sign(JHkitA_r) * B_ir_pre
            B[columns_Xi, r] <- B_ir
        }
    }
    
    return(B)
}

#' Updating loading matrix B with the composite concave penalty
#'
#' This is an intermediate step of the algorithm for fitting pESCA model. The
#' details of this function can be found in ref thesis.
#'
#' @inheritParams update_B_L2
#'
#' @return This function returns the updated loading matrix B.
#'
#' @examples
#' \dontrun{
#' B <- update_B_composite(JHk,A,B0,Sigmas0,d,
#'                    fun_concave,alphas,rhos,lambdas,gamma)
#' }
update_B_composite <- function(JHk, A, B0, Sigmas0, d, fun_concave, alphas, rhos, lambdas, gamma) {
    sumd <- sum(d)
    nDataSets <- length(d)
    n <- dim(A)[1]
    R <- dim(A)[2]
    hfun_sg <- get(paste0(fun_concave, "_sg"))  # super gradient of the concave function
    
    B <- matrix(NA, sumd, R)
    for (i in 1:nDataSets) {
        columns_Xi <- index_Xi(i, d)
        JHk_i <- JHk[, columns_Xi]
        JHkitA <- crossprod(JHk_i, A)
        alpha_i <- alphas[i]
        rho_i <- rhos[i]
        weight_i <- d[i]  # weight for composite concave function
        lambda_i <- lambdas[i] * weight_i * alpha_i/rho_i
        
        for (r in 1:R) {
            # form weights of the penalty according to previous sigma0_ir
            sigma0_ir <- Sigmas0[i, r]
            sigma0_ir_vec <- abs(B0[columns_Xi, r])
            
            omega_ir_outer <- hfun_sg(sigma0_ir, gamma = gamma, lambda = 1)
            omega_ir_inner <- hfun_sg(sigma0_ir_vec, gamma = gamma, lambda = 1)
            omega_ir_vec <- omega_ir_outer * omega_ir_inner
            
            # proximal operator of L1 norm
            JHkitA_r <- JHkitA[, r]
            lambda_ir_vec <- lambda_i * omega_ir_vec
            JHkitA_r <- JHkitA[, r]
            
            B_ir_pre <- abs(JHkitA_r) - lambda_ir_vec
            B_ir_pre[B_ir_pre < 0] <- 0
            B_ir <- sign(JHkitA_r) * B_ir_pre
            B[columns_Xi, r] <- B_ir
        }
    }
    
    return(B)
}

#' Updating loading matrix B with the element-wise concave penalty
#'
#' This is an intermediate step of the algorithm for fitting pESCA model. The
#' details of this function can be found in ref thesis.
#'
#' @inheritParams update_B_L2
#'
#' @return This function returns the updated loading matrix B.
#'
#' @examples
#' \dontrun{
#' B <- update_B_element(JHk,A,B0,Sigmas0,d,
#'                    fun_concave,alphas,rhos,lambdas,gamma)
#' }
update_B_element <- function(JHk, A, B0, Sigmas0, d, fun_concave, alphas, rhos, lambdas, gamma) {
    sumd <- sum(d)
    nDataSets <- length(d)
    n <- dim(A)[1]
    R <- dim(A)[2]
    hfun_sg <- get(paste0(fun_concave, "_sg"))  # super gradient of the concave function
    
    B <- matrix(NA, sumd, R)
    for (i in 1:nDataSets) {
        columns_Xi <- index_Xi(i, d)
        JHk_i <- JHk[, columns_Xi]
        JHkitA <- crossprod(JHk_i, A)
        alpha_i <- alphas[i]
        rho_i <- rhos[i]
        weight_i <- 1  # weight element-wise penalty
        lambda_i <- lambdas[i] * weight_i * alpha_i/rho_i
        
        for (r in 1:R) {
            # form weights of the penalty according to previous sigma0_ir
            sigma0_ir_vec <- abs(B0[columns_Xi, r])
            omega_ir_vec <- hfun_sg(sigma0_ir_vec, gamma = gamma, lambda = 1)  # weights
            
            # proximal operator of L1 norm
            JHkitA_r <- JHkitA[, r]
            lambda_ir_vec <- lambda_i * omega_ir_vec
            
            B_ir_pre <- abs(JHkitA_r) - lambda_ir_vec
            B_ir_pre[B_ir_pre < 0] <- 0
            B_ir <- sign(JHkitA_r) * B_ir_pre
            B[columns_Xi, r] <- B_ir
        }
    }
    
    return(B)
}

#' Group-wise conave L2 norm penalty
#'
#' This is an intermediate step of the algorithm for fitting pESCA model. The
#' details of this function can be found in ref thesis.
#'
#' @param B_i The loading matrix for the ith data set
#' @param fun_concave A string indicates the used concave function
#' @param gamma The hyper-parameter of the concave penalty
#' @param R The number of PCs
#'
#' @return This function returns the value of the
#' group-wise conave L2 norm penalty for the pESCA model
#'
#' @examples
#' \dontrun{
#' penalty_concave_L2(B_i, fun_concave, gamma, R)
#' }
penalty_concave_L2 <- function(B_i, fun_concave, gamma, R) {
    weight_i <- sqrt(dim(B_i)[1])  # weight when L2 norm is used
    hfun <- get(fun_concave)  # name of concave function
    
    out <- 0
    sigmas <- matrix(data = 0, 1, R)
    for (r in 1:R) {
        sigma_ir <- norm(B_i[, r], "2")  # sigma_{lr} = ||b_{lr}||_2
        sigmas[1, r] <- sigma_ir
        
        out <- out + hfun(sigma_ir, gamma = gamma, lambda = 1)
    }
    out <- weight_i * out
    
    result <- list()
    result$sigmas <- sigmas
    result$out <- out
    
    return(result)
}

#' Group-wise conave L1 norm penalty
#'
#' This is an intermediate step of the algorithm for fitting pESCA model. The
#' details of this function can be found in ref thesis.
#'
#' @inheritParams penalty_concave_L2
#'
#' @return This function returns the value of the
#' group-wise conave L1 norm penalty for the pESCA model
#'
#' @examples
#' \dontrun{
#' penalty_concave_L1(B_i, fun_concave, gamma, R)
#' }
penalty_concave_L1 <- function(B_i, fun_concave, gamma, R) {
    weight_i <- dim(B_i)[1]  # weight when L1 norm is used
    hfun <- get(fun_concave)  # concave penalty function
    
    out <- 0
    sigmas <- matrix(data = 0, 1, R)
    for (r in 1:R) {
        sigma_ir <- sum(abs(B_i[, r]))  # sigma_{lr} = ||b_{lr}||_1
        sigmas[1, r] <- sigma_ir
        
        out <- out + hfun(sigma_ir, gamma = gamma, lambda = 1)
    }
    out <- weight_i * out
    
    result <- list()
    result$sigmas <- sigmas
    result$out <- out
    
    return(result)
}

#' Composition of group-wise and element-wise conave penalty
#'
#' This is an intermediate step of the algorithm for fitting pESCA model. The
#' details of this function can be found in ref thesis.
#'
#' @inheritParams penalty_concave_L2
#'
#' @return This function returns the value of the composition
#' of group-wise and element-wise conave penalty for the pESCA model
#'
#' @examples
#' \dontrun{
#' penalty_concave_composite(B_i, fun_concave, gamma, R)
#' }
penalty_concave_composite <- function(B_i, fun_concave, gamma, R) {
    weight_i <- dim(B_i)[1]  # weight when L1 norm is used
    hfun <- get(fun_concave)  # concave penalty function
    
    out <- 0
    sigmas <- matrix(data = 0, 1, R)
    for (r in 1:R) {
        sigma_ir_vec <- hfun(abs(B_i[, r]), gamma = gamma, lambda = 1)  # inner penalty
        sigma_ir <- sum(sigma_ir_vec)
        sigmas[1, r] <- sigma_ir
        
        out <- out + hfun(sigma_ir, gamma = gamma, lambda = 1)  # outer penalty
    }
    out <- weight_i * out
    
    result <- list()
    result$sigmas <- sigmas
    result$out <- out
    
    return(result)
}

#' Element-wise conave penalty
#'
#' This is an intermediate step of the algorithm for fitting pESCA model. The
#' details of this function can be found in ref thesis.
#'
#' @inheritParams penalty_concave_L2
#'
#' @return This function returns the value of the
#' element-wise conave penalty for the pESCA model
#'
#' @examples
#' \dontrun{
#' penalty_concave_element(B_i, fun_concave, gamma, R)
#' }
penalty_concave_element <- function(B_i, fun_concave, gamma, R) {
    weight_i <- 1  # weight when element-wise penalty is used
    hfun <- get(fun_concave)  # concave penalty function
    
    out <- 0
    sigmas <- matrix(data = 0, 1, R)
    for (r in 1:R) {
        sigma_ir_vec <- abs(B_i[, r])
        sigma_ir <- sum(sigma_ir_vec)  # sigma_{lr} = ||b_{lr}||_1
        sigmas[1, r] <- sigma_ir
        
        out <- out + sum(hfun(sigma_ir_vec, gamma = gamma, lambda = 1))
    }
    out <- weight_i * out
    
    result <- list()
    result$sigmas <- sigmas
    result$out <- out
    
    return(result)
}
