#' Simulate multiple data sets with element-wise sparse patterns
#'
#' This function can generate multiple data sets of mixed data types.
#' Details of the simultion process can be found in the paper
#' \url{https://arxiv.org/abs/1902.06241}. In this function, only the
#' element-wise sparsity is included into the loading matrix.
#'
#' @inheritParams dataSimu_group_sparse
#' @param R the number of simulated PCs
#' @param SNRs the SNRs for the simulation of the multiple data sets
#'
#' @return Refer the return of thefunction \code{dataSimu_group_sparse}.
#'
#' @importFrom stats rbeta
#' @importFrom stats rlogis
#' @importFrom stats rnorm
#' @importFrom stats runif
#'
#' @examples
#' \dontrun{
#' dataSimulation <- dataSimu_element_sparse(n,ds,R,
#'                                    dataTypes='GGG',
#'                                    noises=rep(1,3),
#'                                    margProb=0.1,sparse_ratio=0.5,
#'                                    SNRs=rep(1,3))
#'}
#'
#' @export
dataSimu_element_sparse <- function(n, ds, R = 3, 
                                    dataTypes = "GGG", 
                                    noises = rep(1, 3), 
                                    margProb = 0.1, 
                                    sparse_ratio = 0.5, 
                                    SNRs = rep(1, 3)) {
    # change dataTypes from a string to a vector of string
    if (length(dataTypes) == 1) {
        dataTypes <- unlist(strsplit(dataTypes, split = ""))
    }
    if (length(dataTypes) != length(ds)) {
        stop("The number of specified dataTypes are not equal to the number of data sets")
    }
    
    # parameters used in the simulation
    sumd <- sum(ds)  # sum of variables
    
    # simulate the offset term
    mu_simu <- matrix(data = 0, nrow = 1, ncol = sumd)
    for (i in 1:3) {
        # generate index for the ith data set
        columns_Xi <- index_Xi(i, ds)
        
        # simulate the offset for the ith data set
        if (dataTypes[i] == "G") {
            # from normal distribution
            mu_simu[1, columns_Xi] <- stats::rnorm(ds[i])
        } else if (dataTypes[i] == "B") {
            # logit transform of beta distribution
            beta_a <- round(n * margProb) + 1
            beta_b <- n - round(n * margProb) + 1
            mu_simu[1, columns_Xi] <- logit(stats::rbeta(ds[i], beta_a, beta_b))
        }
    }
    
    # the simulation of U_simu, D_pre, V_pre
    U_pre <- matrix(stats::rnorm(n * R), nrow = n, ncol = R)
    U_simu <- svd(scale(U_pre, center = TRUE, scale = FALSE), R)$u
    D_pre <- abs(1 + 0.5 * stats::rnorm(R))
    V_pre <- matrix(stats::rnorm(sumd * R), nrow = sumd, ncol = R)
    V_pre <- svd(V_pre, R)$u
    
    # modify V_pre to have element-wise sparsity
    V_simu <- V_pre
    for (i in 1:3) {
        # generate index for the ith data set
        columns_Xi <- index_Xi(i, ds)
        
        # induce element-wise saprsity on the rth column of the ith data set
        for (r in 1:R) {
            V_ir <- V_simu[columns_Xi, r]
            
            # element-wise sparsity: randomly set sparse_ratio% elements to be 0
            index_tmp <- sample(1:ds[i], round(sparse_ratio * ds[i]))
            V_ir[index_tmp] <- 0
            V_simu[columns_Xi, r] <- V_ir
        }
    }
    
    # simulate noise terms
    E_simu <- matrix(data = 0, nrow = n, ncol = sumd)
    for (i in 1:3) {
        columns_Xi <- index_Xi(i, ds)
        
        if (dataTypes[i] == "G") {
            # from normal distribution
            E_simu[1:n, columns_Xi] <- matrix(noises[i] * stats::rnorm(n * ds[i]), n, ds[i])
        } else if (dataTypes[i] == "B") {
            # from logistic distribution
            E_simu[1:n, columns_Xi] <- matrix(noises[i] * stats::rlogis(n * ds[i]), n, ds[i])
        }
    }
    
    # Modify the D_pre as C.*D to satisfy the pre-defined SNR
    C_factors <- rep(0, R)
    Theta_pre <- U_simu %*% diag(D_pre) %*% t(V_simu)
    
    for (i in 1:3) {
        columns_Xi <- index_Xi(i, ds)
        Theta_pre_i <- Theta_pre[, columns_Xi]
        E_i <- E_simu[, columns_Xi]
        
        # compute the proper scale for SNR
        C_factors[i] <- sqrt(SNRs[i]) * norm(E_i, "F")/norm(Theta_pre_i, "F")
    }
    
    # simulate D_simu
    D_simu <- C_factors * D_pre
    
    # simulate Theta_simu
    mu_simu <- t(mu_simu)  # tansform mu_simu to column vector form
    Theta_simu <- ones(n) %*% t(mu_simu) + U_simu %*% (diag(D_simu)) %*% t(V_simu)
    
    # simulate quantitative or binary data sets
    X <- list()
    
    for (i in 1:3) {
        columns_Xi <- index_Xi(i, ds)
        
        if (dataTypes[i] == "G") {
            tmp_Xi <- Theta_simu[, columns_Xi] + E_simu[, columns_Xi]
        } else if (dataTypes[i] == "B") {
            X_tmp <- sign(inver_logit(Theta_simu[, columns_Xi]) - matrix(stats::runif((n * ds[i])), n, 
                ds[i]))
            tmp_Xi <- 0.5 * (X_tmp + 1)
        }
        
        X[[i]] <- tmp_Xi
    }
    
    # latent data sets
    Xstar <- Theta_simu + E_simu
    
    # variance explained for each data set
    B_simu <- V_simu %*% diag(D_simu)
    varExpTotals_simu <- matrix(data = 0, nrow = 1, ncol = 4)
    varExpPCs_simu <- matrix(data = 0, nrow = 4, ncol = R)
    
    for (i in 1:3) {
        columns_Xi <- index_Xi(i, ds)
        Xi <- Xstar[, columns_Xi]
        mu_Xi <- mu_simu[columns_Xi, 1]
        B_Xi <- B_simu[columns_Xi, ]
        W_Xi <- matrix(data = 1, nrow = n, ncol = ds[i])
        
        varExp_tmp <- varExp_Gaussian(Xi, as.vector(mu_Xi), U_simu, B_Xi, W_Xi)
        varExpTotals_simu[i] <- varExp_tmp$varExp_total
        varExpPCs_simu[i, ] <- varExp_tmp$varExp_PCs
    }
    
    # variance explained ratio combine all the data sets noise levels are assumed to be equal here
    W = matrix(data = 1, nrow = n, ncol = sumd)
    varExp_tmp <- varExp_Gaussian(Xstar, as.vector(mu_simu), U_simu, B_simu, W)
    varExpTotals_simu[4] <- varExp_tmp$varExp_total
    varExpPCs_simu[4, ] <- varExp_tmp$varExp_PCs
    names(varExpTotals_simu) <- paste0("X_", c(1:3, "full"))
    rownames(varExpPCs_simu) <- paste0("X_", c(1:3, "full"))
    colnames(varExpPCs_simu) <- paste0("PC", 1:R)
    
    # save the results
    out <- list()
    
    # save X, Theta, E
    out$X <- X
    out$Theta_simu <- Theta_simu
    out$E_simu <- E_simu
    out$dataTypes <- dataTypes
    out$SNRs <- SNRs
    
    # save mu, U, D, V
    out$mu_simu <- mu_simu
    out$U_simu <- U_simu
    out$D_simu <- D_simu
    out$V_simu <- V_simu
    
    # save true variation explained
    out$varExpTotals_simu <- varExpTotals_simu
    out$varExpPCs_simu <- varExpPCs_simu
    
    return(out)
}

#' Test data simulation only with element-wise sparse pattern
#'
#' This function will test if the simulated data sets have the proper
#' properties, such as correct data types, orghogonal score matrix,
#' score matrix has 0 column means, element-wise sparsity, with correct SNRs.
#'
#' @inheritParams test_dataSimu_group
#' @param R number of simulated PCs
#'
#' @return This function returns a list contains all the test results
#'
#' @examples
#' \dontrun{
#' test_results_GGG <- test_dataSimu_element(n=n, ds=ds,R=3
#'                                 dataTypes = 'GGG')
#' test_results_BBB <- test_dataSimu_element(n=n, ds=ds,R=3
#'                                 dataTypes = 'BBB')
#' }
test_dataSimu_element <- function(n, ds, R, dataTypes) {
    
    # data simulation
    simulatedData <- dataSimu_element_sparse(n, ds, R, dataTypes)
    
    # test simulated data types are correct
    dataTypes_object <- simulatedData$dataTypes
    dataTypes_test <- rep("G", length(ds))
    
    for (i in 1:length(ds)) {
        Xi <- simulatedData$X[[i]]
        unique_values <- unique(as.vector(Xi))
        if (all(unique_values %in% c(1, 0))) {
            dataTypes_test[i] <- "B"
        }
    }
    
    test_dataTypes <- all(dataTypes_test == dataTypes_object)
    
    # test the orghogality of simulated score matrix
    sumR <- dim(simulatedData$U_simu)[2]
    UtU <- t(simulatedData$U_simu) %*% simulatedData$U_simu
    test_ortho <- all.equal(diag(UtU), rep(1, sumR)) & all.equal(UtU - diag(diag(UtU)), matrix(data = 0, 
        sumR, sumR))
    
    # test if the column mean of the score matrix is 0
    test_means <- all.equal(colMeans(simulatedData$U_simu), rep(0, sumR))
    
    # test if desiered sparsity are included
    S_test <- colMeans((simulatedData$V_simu) == 0)
    test_element_sparse <- all.equal(S_test, rep(0.5, 3))
    
    result <- list()
    result$test_dataTypes <- test_dataTypes
    result$test_ortho <- test_ortho
    result$test_means <- test_means
    result$test_element_sparse <- test_element_sparse
    
    return(result)
}

