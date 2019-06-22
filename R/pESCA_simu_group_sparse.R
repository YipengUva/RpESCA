#' Simulate multiple data sets with group and element-wise sparse patterns
#'
#' This function can generate multiple data sets of mixed data types.
#' Details of the simultion process can be found in the paper
#' \url{https://arxiv.org/abs/1902.06241}. In the current version, both
#' group and element-wise sparsity are included into the loading matrix.
#'
#' @param n the number of objects
#' @param ds a vector for the number of variables in each data set
#' @param dataTypes a string indicates the data type of each data set,
#' possible options include 'G': Gaussian, 'B': Bernoulli.
#' @param noises noise levels of simulated data sets
#' @param margProb desired marginal probability for binary data simulation,
#' used to simulate imbalanced binary data.
#' @param sparse_ratio controls the sparse level of element-wise sparsity
#' @param SNRgc SNR for global common structure
#' @param SNRlc SNRs of the local common structures
#' @param SNRd SNRs of the distinct structures
#'
#' @return A list contains the simulated data sets and the simulated parameters. \itemize{
#' \item X: a list contains the simulated multiple data sets;
#' \item Theta_simu: simulated natural parameter matrix Theta;
#' \item mu_simu: simulated offset term mu;
#' \item U_simu: simulated U;
#' \item D_simu: simulated D;
#' \item V_simu: simulated V;
#' \item E_simu: simulated noise term E;
#' \item dataTypes: vector form simulated data types;
#' \item S_simu: desired group sparse pattern;
#' \item SNRs: used SNRs in simulating common (global, local) and distinct structures;
#' \item varExpTotals_simu: variation explained ratios for each data set computed
#' using the simulated parameters;
#' \item varExpPCs_simu: variation explained raitos for each PC for each data set
#' computed using the simulated parameters;
#' }
#'
#' @importFrom stats rbeta
#' @importFrom stats rlogis
#' @importFrom stats rnorm
#' @importFrom stats runif
#'
#' @examples
#' \dontrun{
#' dataSimulation <- dataSimu_group_sparse(n=200,
#'                                        ds=c(400,200,100),
#'                                        dataTypes='GGB')
#' }
#'
#' @export
dataSimu_group_sparse <- function(n, ds, dataTypes = "GGG", noises = rep(1, 3), margProb = 0.1, sparse_ratio = 0, 
    SNRgc = 1, SNRlc = rep(1, 3), SNRd = rep(1, 3)) {
    
    # parameters used in the simulation
    if (length(dataTypes) == 1) {
        dataTypes <- unlist(strsplit(dataTypes, split = ""))
    }
    if (length(dataTypes) != length(ds)) {
        stop("The number of specified dataTypes are not equal to the number of data sets")
    }
    Rgc <- 3  # number of PCs for the global common structure
    Rlc <- rep(3, 3)  # number of PCs for the local common structures
    Rd <- rep(3, 3)  # number of PCs for the distint structures
    sumR <- Rgc + sum(Rlc) + sum(Rd)  # sum of PCs
    sumd <- sum(ds)  # sum of variables
    
    # noise levels used in the data simulation process
    SNRs <- c(SNRgc, SNRlc, SNRd)
    
    # simulate the offset term
    mu_simu <- matrix(data = 0, nrow = sumd, ncol = 1)
    for (i in 1:3) {
        # generate index for the ith data set
        columns_Xi <- index_Xi(i, ds)
        
        # simulate the offset for the ith data set
        if (dataTypes[i] == "G") {
            # from normal distribution
            mu_simu[columns_Xi, 1] <- stats::rnorm(ds[i])
        } else if (dataTypes[i] == "B") {
            # logit transform of beta distribution
            beta_a <- round(n * margProb) + 1
            beta_b <- n - round(n * margProb) + 1
            mu_simu[columns_Xi, 1] <- logit(stats::rbeta(ds[i], beta_a, beta_b))
        }
    }
    
    # the simulation of U_simu, D_pre, V_pre
    U_pre <- matrix(stats::rnorm(n * sumR), nrow = n, ncol = sumR)
    U_simu <- svd(scale(U_pre, center = TRUE, scale = FALSE), sumR)$u
    D_pre <- abs(1 + 0.5 * stats::rnorm(sumR))
    V_pre <- matrix(stats::rnorm(sumd * sumR), nrow = sumd, ncol = sumR)
    V_pre <- svd(V_pre, sumR)$u
    
    # modify V_pre to have the predifed structure the target group sparsity
    S_simu <- cbind(matrix(rep(c(1, 1, 1), Rgc), 3, Rgc), matrix(rep(c(1, 1, 0), Rlc[1]), 3, Rlc[1]), 
        matrix(rep(c(1, 0, 1), Rlc[2]), 3, Rlc[2]), matrix(rep(c(0, 1, 1), Rlc[3]), 3, Rlc[3]), matrix(rep(c(1, 
            0, 0), Rd[1]), 3, Rd[1]), matrix(rep(c(0, 1, 0), Rd[2]), 3, Rd[2]), matrix(rep(c(0, 0, 1), 
            Rd[3]), 3, Rd[3]))
    
    # modify V_pre to have the target group sparsity
    V_simu <- matrix(data = 0, nrow = sumd, ncol = sumR)
    for (i in 1:3) {
        # generate index for the ith data set
        columns_Xi <- index_Xi(i, ds)
        
        # induce group saprsity on the rth column of the ith data set
        for (r in 1:sumR) {
            V_simu[columns_Xi, r] <- S_simu[i, r] * V_pre[columns_Xi, r]
        }
    }
    
    # impose element-wise sparsity
    for (i in 1:3) {
        # generate index for the ith data set
        columns_Xi <- index_Xi(i, ds)
        
        # induce element-wise saprsity on the rth column of the ith data set
        for (r in 1:sumR) {
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
    C_factors <- rep(0, sumR)
    
    for (i in 1:length(SNRs)) {
        index_factors <- (3 * (i - 1) + 1):(3 * i)
        
        Theta_pre <- U_simu[, index_factors] %*% diag(D_pre[index_factors]) %*% t(V_simu[, index_factors])
        index_variables <- colSums(abs(Theta_pre)) > 0
        Theta_pre <- Theta_pre[, index_variables]
        E_pre <- E_simu[, index_variables]
        
        # compute the proper scale for SNR
        C_factors[index_factors] <- sqrt(SNRs[i]) * norm(E_pre, "F")/norm(Theta_pre, "F")
    }
    
    # simulate D_simu
    D_simu <- C_factors * D_pre
    
    # simulate Theta_simu
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
    varExpPCs_simu <- matrix(data = 0, nrow = 4, ncol = sumR)
    
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
    colnames(varExpPCs_simu) <- paste0("PC", 1:sumR)
    
    # save the results
    out <- list()
    
    # save X, Theta, E
    out$X <- X
    out$Theta_simu <- Theta_simu
    out$E_simu <- E_simu
    out$dataTypes <- dataTypes
    out$S_simu <- S_simu
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

#' Test data simulation with group sparse pattern
#'
#' This function will test if the simulated data sets have the proper
#' properties, such as correct data types, orghogonal score matrix,
#' score matrix has 0 column means, group sparsity, with correct SNRs.
#'
#' @param n the number of objects
#' @param ds a vector for the number of variables in each data set
#' @param dataTypes a string indicates the data type of each data set,
#' possible options include 'G': Gaussian, 'B': Bernoulli.
#'
#' @return This function returns a list contains all the test results
#'
#' @examples
#' \dontrun{
#' test_results_GGG <- test_dataSimu_group(n=n, ds=ds,
#'                                 dataTypes = 'GGG')
#' test_results_BBB <- test_dataSimu_group(n=n, ds=ds,
#'                                 dataTypes = 'BBB')
#' }
test_dataSimu_group <- function(n, ds, dataTypes) {
    
    # data simulation
    simulatedData <- dataSimu_group_sparse(n, ds, dataTypes)
    
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
    S_object <- simulatedData$S_simu
    S_test <- matrix(data = 0, 3, sumR)
    S_test[simulatedData$varExpPCs_simu[1:3, ] > 0] <- 1
    test_group_sparse <- all(S_test == S_object)
    
    # test if correct SNR is defined
    SNRs_object <- simulatedData$SNRs
    SNRs_test <- rep(0, length(SNRs_object))
    
    for (i in 1:length(SNRs_object)) {
        index_factors <- (3 * (i - 1) + 1):(3 * i)
        
        Theta_factors <- simulatedData$U_simu[, index_factors] %*% diag(simulatedData$D_simu[index_factors]) %*% 
            t(simulatedData$V_simu[, index_factors])
        
        index_variables <- colSums(abs(Theta_factors)) > 0
        Theta_factors <- Theta_factors[, index_variables]
        E_factors <- simulatedData$E_simu[, index_variables]
        
        # compute the proper scale for SNR
        SNRs_test[i] <- norm(Theta_factors, "F")/norm(E_factors, "F")
    }
    
    test_SNRs <- all.equal(SNRs_test, SNRs_object)
    
    result <- list()
    result$test_dataTypes <- test_dataTypes
    result$test_ortho <- test_ortho
    result$test_means <- test_means
    result$test_group_sparse <- test_group_sparse
    result$test_SNRs <- test_SNRs
    
    return(result)
}







