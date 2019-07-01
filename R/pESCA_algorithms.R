#' Penalized exponential family simultaneous component analysis (pESCA) model
#'
#' This is the main function for construncting a pESCA model on multiple data
#' sets. The potential different data types in these data sets are tackled by
#' the assumption of exponential family distribution. Gaussian for quantitative
#' data, Bernoulli for binary data and Poisson for count data. Although the option
#' for count data using Poisson distribution is included in the algorithm, we recommend
#' to do variance stabilizing transformation on the count data, such as Next-Gen
#' sequencing data, and then use the transformed data as quantitative data sets. The
#' details of the developed algorithm can be found in \url{https://arxiv.org/abs/1902.06241}.
#'
#' @param dataSets a list contains multiple matrices with same number of rows.
#' Each matrix (\code{samples * variables}) indicates a data set.
#' @param dataTypes a string indicates the data types of the multiple data sets.
#' @param lambdas a numeric vector indicates the values of tuning parameters for
#' each data set.
#' @param penalty The name of the penalty you want to used. \itemize{
#' \item "L2": group-wise concave L2 norm penalty;
#' \item "L1": group-wise concave L1 norm penalty;
#' \item "element": element-wise concave penalty;
#' \item "composite": the composition of group- and element-wise penalty.
#' }
#' @param fun_concave a string indicates the used concave function. Three options
#' are included in the algorithm. \itemize{
#' \item "gdp": GDP penalty;
#' \item "lq": Lq penalty;
#' \item "scad": SCAD penalty.
#' }
#' @param opts a list contains the options of the algorithms. \itemize{
#' \item tol_obj: tolerance for relative change of joint loss function, default:1E-6;
#' \item maxit: max number of iterations, default: 1000;
#' \item gamma: hyper-parameter of the concave penalty, default: 1;
#' \item R: the initial number of PCs, default: 0.5 \code{0.5*min(I,J)};
#' \item rand_start: initilization method, random (1), SCA(0), default: 0;
#' \item alphas: dispersion parameters of exponential dispersion families, default: 1.
#' \item thr_path: the option to generate thresholding path, default: 0;
#' \item quiet: quiet==1, not show the progress when running the algorithm, default: 0.
#' }
#'
#' @return This function returns a list contains the results of a pESCA mdoel. \itemize{
#' \item mu: offset term;
#' \item A: score matrix;
#' \item B: loading matrix;
#' \item S: group sparse pattern of B;
#' \item varExpTotals: total variation explained of each data set and the full data set;
#' \item varExpPCs: variation explained of each data set and each component;
#' \item Sigmas: the group length (the definition depends on the used input type). Only meaningful
#' for group-wise sparisty;
#' \item iter: number of iterations;
#' \item diagnose$hist_objs: the value of loss function pluse penalty at each iteration;
#' \item diagnose$f_objs: the value of loss function at each iteration;
#' \item diagnose$g_objs: the value of penalty function at each iteration;
#' \item diagnose$rel_objs: relative change of diagnose$hist_obj at each iteration;
#' \item diagnose$rel_Thetas: relative change of Theta at each iteration.
#' }
#'
#' @importFrom RSpectra svds
#'
#' @examples
#' \dontrun{
#' # Suppose we have three data sets X1, X2, X3
#' # They are quantitative, quantitative and binary matrices
#' pESCA(dataSets = list(X1, X2, X3),
#'               dataTypes = 'GGB',
#'               lambdas = c(20, 20, 10),
#'               penalty = 'L2',
#'               fun_concave = 'gdp',
#'               opts = list())
#' }
#'
#' @export
pESCA <- function(dataSets, dataTypes,
                          lambdas, penalty='L2', fun_concave='gdp', opts=list()){
  # check if data sets are in a list
  stopifnot(class(dataSets) == "list") 
  
  # check if the used penalties are included in the algorithm
  possible_penalties <- c("L2", "L1", "composite", "element")
  if ( !(penalty %in% possible_penalties) )
    stop("The used penalty is not included in the algorithm")
  
  # check if the used concave function are included in the algorithm
  possible_funs <- c("gdp", "lq", "scad")
  if ( !(fun_concave %in% possible_funs) )
    stop("The used concave function is not included in the algorithm")
  
  # check if the numer of data sets are equal to the number of lambdas
  if ( !(length(dataSets) == length(lambdas) ) )
    stop("the number of data sets are not equal to the number of lambdas")
  
  # check if the used data Types are included in the algorithm
  if(length(dataTypes) == 1) dataTypes <- unlist(strsplit(dataTypes, split=""))
  if(length(dataTypes) != length(dataSets) )
    stop("The number of specified dataTypes are not equal to the number of data sets")

  # default parameters
  if(exists('tol_obj', where = opts)){tol_obj <- opts$tol_obj} else{tol_obj <- 1E-6};
  if(exists('maxit', where = opts)){maxit <- opts$maxit} else{maxit <- 1000};
  if(exists('gamma', where = opts)){gamma <- opts$gamma} else{gamma <- 1};
  if(exists('rand_start', where = opts)){rand_start <- opts$rand_start
  }else{rand_start <- 0};
  if(exists('thr_path', where = opts)){thr_path <- opts$thr_path} else{thr_path <- 0};
  if(exists('quiet', where = opts)){quiet <- opts$quiet} else{quiet <- 0};

  # number of data sets, size of each data set
  nDataSets <- length(dataSets) # number of data sets
  n <- rep(0, nDataSets)  # number of samples
  d <- rep(0, nDataSets)  # numbers of variables in different data sets
  for(i in 1:nDataSets){
    n[i] <- dim(dataSets[[i]])[1]
    d[i] <- dim(dataSets[[i]])[2]
  }
  if(length(unique(n))!=1)
    stop("multiple data sets have unequal sample size")
  n <- n[1]
  sumd <- sum(d) # total number of variables

  # form full data set X, X = [X1,...Xl,...XL]
  # form full weighting matrix W, W = [W1,...Wl,...WL]
  # test dataTypes
  X <- matrix(data=NA, nrow=n, ncol=sumd)
  W <- matrix(data=NA, nrow=n, ncol=sumd)
  test_dataTypes <- rep('G', nDataSets)

  for(i in 1:nDataSets){
    columns_Xi <- index_Xi(i, d)
    X_i <- as.matrix(dataSets[[i]])
    W_i <- matrix(data=1, nrow=n, ncol=d[i])
    W_i[is.na(X_i)] <- 0
    X_i[is.na(X_i)] <- 0

    X[,columns_Xi] <- X_i
    W[,columns_Xi] <- W_i

    unique_values <- unique(as.vector(X_i))
    if((length(unique_values) == 2) & all(unique_values %in% c(1,0))){
      test_dataTypes[i] <- 'B'}
  }
  if(!all(dataTypes == test_dataTypes))
    warning("The specified data types may not correct")

  # default dispersion parameters alphas and number of PCs
  if(exists('alphas', where=opts)){alphas<-opts$alphas} else{alphas<-rep(1,nDataSets)};
  if(exists('R', where=opts)){R <- opts$R} else {R <- round(0.5*min(c(n,d)))};

  # initialization
  if(exists('A0', where=opts)){ # use imputted initialization
    mu0 <- opts$mu0; mu0 <- t(mu0);
    A0 <- opts$A0
    B0 <- opts$B0
    R <- dim(B0)[2]
  } else if(rand_start == 1){ # use random initialization
    mu0 <- matrix(data=0, 1, sumd)
    A0 <- matrix(rnorm(n*R), n, R)
    B0 <- matrix(rnorm(sumd*R), sumd, R)
  } else if(rand_start == 0){ # use SCA model as initialization
    if(!('P' %in% dataTypes)){ # Possition distribution is not used
      mu0 <- matrix(colMeans(X), 1, sumd)
      X_svd <- RSpectra::svds(scale(X,center=TRUE,scale=FALSE),R,nu=R,nv=R)
      A0 <- X_svd$u
      B0 <- X_svd$v %*% diag(X_svd$d[1:R])
    } else{
      X_tmp <- X
      for(i in 1:nDataSets){ #log transformation applied to Possion data
        dataType <- dataTypes[i]
        if(dataType == 'P'){
          columns_Xi <- index_Xi(i,d)
          X_tmp[,columns_Xi] <- log(X[,columns_Xi] + 1)
        }}
      mu0 <- matrix(colMeans(X_tmp), 1, sumd)
      X_svd <- RSpectra::svds(scale(X_tmp,center=TRUE,scale=FALSE),R,nu=R,nv=R)
      A0 <- X_svd$u
      B0 <- X_svd$v %*% diag(X_svd$d[1:R])
    }
  }
  Theta0 = ones(n) %*% mu0 + tcrossprod(A0, B0)

  # initial value of loss function
  # specify penalty function
  penalty_fun <- get(paste0("penalty_concave_",penalty))
  update_B_fun <- get(paste0("update_B_",penalty))

  f_obj0 <- 0
  g_obj0 <- 0
  Sigmas0 <- matrix(data=0, nDataSets, R)
  for(i in 1:nDataSets){
    dataType <- dataTypes[i]
    columns_Xi <- index_Xi(i,d)
    X_i <- X[, columns_Xi]
    W_i <- W[, columns_Xi]
    Theta0_i <- Theta0[, columns_Xi]
    alpha_i <-  alphas[i]
    lambda_i <- lambdas[i]

    # loss function for fitting ith data set
    # specify log-partiton function for ith data set
    log_partition <- get(paste0("log_part_", dataType))
    f_obj0 <- {f_obj0 + (1/alpha_i)*
        (trace_fast(W_i, log_partition(Theta0_i)) - trace_fast(Theta0_i, X_i))}

    # penalty for the ith loading matrix B_l
    B0_i <- B0[columns_Xi, ]
    g_penalty <- penalty_fun(B0_i, fun_concave, gamma, R)
    g_obj0 <- g_obj0 + lambda_i*g_penalty$out
    Sigmas0[i,] <- g_penalty$sigmas
  }

  # record loss and penalty for diagnose purpose
  diagnose <- list()
  diagnose$hist_objs <- vector()
  diagnose$f_objs <- vector()
  diagnose$g_objs <- vector()
  diagnose$rel_objs <- vector()
  diagnose$rel_Thetas <- vector()

  #inital values for loss function and penalty
  obj0 <- f_obj0 + g_obj0   # objective + penalty
  diagnose$f_objs[1] <- f_obj0 # objective
  diagnose$g_objs[1] <- g_obj0 # penalty
  diagnose$hist_objs[1] <- obj0

  # iterations
  for (k in 1:maxit){
    if(quiet == 0) print(paste(k, 'th iteration'))

    # majorizaiton step for p_ESCA model
    #--- form Hk ---
    #--- update mu ---
    #--- form JHk ---
    JHk <- matrix(data=NA, n, sumd)
    mu <- matrix(data=NA, 1, sumd)
    rhos <- rep(NA, nDataSets) # form rhos, the Lipshize constant for each data types
    cs <- rep(NA, sumd) # scaling factors

    for(i in 1:nDataSets){
      dataType   <- dataTypes[i]
      columns_Xi <- index_Xi(i,d)
      X_i <- X[, columns_Xi]
      W_i <- W[, columns_Xi]
      Theta0_i <- Theta0[, columns_Xi]

      # specify the gradient of the log-partiton function
      log_partition_g <- get(paste0("log_part_",dataType,"_g"))

      # form rhos, the Lipshize constant for each data types
      if(dataType == 'G'){
        rho_i <- 1
      }else if(dataType == 'B'){
        rho_i <- 0.25
      }else if(dataType == 'P'){
        rho_i <- max(exp(Theta0_i))
      }
	     rhos[i] <- rho_i

      # form Hk_i
      Hk_i <- Theta0_i - (1/rho_i) * (W_i * (log_partition_g(Theta0_i) - X_i))

      # update mu_i
      mu_i <- colMeans(Hk_i)
      mu[1, columns_Xi] <- mu_i

      # form JHk_i
      JHk[, columns_Xi] <- scale(Hk_i, center=TRUE, scale=FALSE)

      # form scaling factors for scaled_JHk_i, scaled_Bk_i
      alpha_i <- alphas[i]
      cs[columns_Xi]<- rep(sqrt(rho_i/alpha_i), d[i])
    }
 
    # update A
    #A_tmp <- JHk %*% (diag(cs^2) %*% B0)
    A_tmp <- mat_vec_matC(JHk, cs^2, B0)
    A_svd <- svd(A_tmp, nu=R, nv=R)
    A <- tcrossprod(A_svd$u, A_svd$v)

    # update B
    B <- update_B_fun(JHk,A,B0,Sigmas0,d,
                      fun_concave,alphas,rhos,lambdas,gamma)

    # diagnostics
    Theta = ones(n)%*%mu + tcrossprod(A,B)

    f_obj <- 0
    g_obj <- 0
    Sigmas <- matrix(data=0, nDataSets, R)
    for(i in 1:nDataSets){
      dataType <- dataTypes[i]
      columns_Xi <- index_Xi(i,d)
      X_i <- X[, columns_Xi]
      W_i <- W[, columns_Xi]
      Theta_i <- Theta[, columns_Xi]
      alpha_i <-  alphas[i]
      weight_i  <- sqrt(d[i]) # weight when L2 norm is used
      lambda_i  <- lambdas[i]

      # loss function for fitting ith data set
      # specify log-partiton function for ith data set
      log_partition <- get(paste0("log_part_",dataType))
      f_obj <- {f_obj + (1/alpha_i)*
          (trace_fast(W_i,  log_partition(Theta_i))-trace_fast(Theta_i,X_i))}

      # penalty for the ith loading matrix B_l
      B_i <- B[columns_Xi,]
      g_penalty <- penalty_fun(B_i, fun_concave, gamma, R)
      g_obj <- g_obj + lambda_i*g_penalty$out
      Sigmas[i,] <- g_penalty$sigmas
      }
    obj <- f_obj + g_obj   # objective + penalty

    # reporting
    diagnose$hist_objs[k+1] <- obj   # objective + penalty
    diagnose$f_objs[k+1]    <- f_obj; # objective
    diagnose$g_objs[k+1]   <- g_obj; # penalty

    diagnose$rel_objs[k] <- (obj0-obj)/abs(obj0) # relative change of loss function
    # relative change of parameters
    diagnose$rel_Thetas[k]  <- norm(Theta0-Theta,'F')^2/norm(Theta0,'F')^2

    # remove the all zeros columns to simplify the computation and save memory
    if(thr_path == 0){
      nonZeros_index <- (colMeans(Sigmas) > 0)
      if(sum(nonZeros_index) > 3){
        A <- A[,nonZeros_index]
        B <- B[,nonZeros_index]
        Sigmas <- Sigmas[,nonZeros_index]
        R <- dim(B)[2]
      }
    }

    # stopping checks
    if((k>1) & ((diagnose$rel_objs[k] < tol_obj))) break

    # save previous results
    B0 <- B
    Theta0  <- Theta
    Sigmas0 <- Sigmas
    obj0    <- obj
  }

  # output
  result <- list()
  result$mu <- t(mu)
  result$A <- A
  result$B <- B
  result$Sigmas <- Sigmas
  result$iter <- k
  result$diagnose <- diagnose
  
  if (!is.null(R)){
    # variation explained ratios
    varExpTotals <- rep(NA, nDataSets + 1)  # +1 is for the full data set
    varExpPCs    <- matrix(NA, nDataSets + 1, R) # +1 is for the full data set
    X_full <- matrix(NA, n, sumd) # combine the quantitative data and the pseudo data of nonGaussian data
    for(i in 1:nDataSets){
      columns_Xi <- index_Xi(i,d)
      Xi    <- X[,columns_Xi]
      mu_Xi <- mu[1,columns_Xi]
      if(!(dataTypes[i]=='G')){
        Xi <- JHk[,columns_Xi] + ones(n)%*%mu_Xi
      }
      X_full[,columns_Xi] <- Xi
      
      B_Xi  <- B[columns_Xi,]
      W_Xi  <- W[,columns_Xi]
      
      varExp_tmp <- varExp_Gaussian(Xi,mu_Xi,A,B_Xi,W_Xi)
      varExpTotals[i] <- varExp_tmp$varExp_total
      varExpPCs[i,] <- varExp_tmp$varExp_PCs
    }
    
    # variance explained ratio of the full data set
    if(length(unique(alphas))==1){
      varExp_tmp <- varExp_Gaussian(X_full,as.vector(mu),A,B,W)
    } else{
      weighted_X  <- matrix(NA, n, sumd)
      weighted_mu <- matrix(NA, 1, sumd)
      weighted_B  <- matrix(NA, sumd, R)
      for (i in 1:nDataSets){
        columns_Xi <- index_Xi(i,d)
        Xi <- X_full[,columns_Xi]
        mu_Xi <- mu[1,columns_Xi]
        B_Xi  <- B[columns_Xi,]
        alpha_i <- alphas[i]
        
        weighted_X[,columns_Xi] <- (1/sqrt(alpha_i))*Xi
        weighted_mu[1,columns_Xi] <- (1/sqrt(alpha_i))*mu_Xi
        weighted_B[columns_Xi,]  <- (1/sqrt(alpha_i))*B_Xi
      }
      varExp_tmp <- varExp_Gaussian(weighted_X,as.vector(weighted_mu),A,weighted_B,W)
    }
    varExpTotals[nDataSets+1] <- varExp_tmp$varExp_total
    varExpPCs[nDataSets+1,] <- varExp_tmp$varExp_PCs
    dataSets_names <- paste0(rep("X_"), c(as.character(1:nDataSets), "full"))
    PCs_names <- paste0("PC", c(as.character(1:R)))
    names(varExpTotals) <- dataSets_names
    rownames(varExpPCs) <- dataSets_names
    colnames(varExpPCs) <- PCs_names
    
    # extract the strcuture index
    S <- matrix(data=0, nDataSets, R)
    S[varExpPCs[1:nDataSets,] > 0] <- 1
    
    # save the results
    result$S <- S
    result$varExpTotals <- varExpTotals
    result$varExpPCs <- varExpPCs
  }
  # output
  return(result)
}
