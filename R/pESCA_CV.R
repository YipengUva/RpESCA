#' pESCA model selection based on cross validation error
#'
#' This function implements a missing value based CV model selection
#' approach for the pESCA model on mutliple data sets with the same data type.
#' The details can be found in  \url{https://arxiv.org/abs/1902.06241}.
#'
#' @inheritParams pESCA
#' @param lambdas_CV a vector cotains a sequence of values of lambda
#'
#' @return This function returns a list contains the results of a pESCA mdoel. \itemize{
#' \item cvErrors_mat: a matrix contains the CV errors for the full data set and each
#' single data set;
#' \item inits: a list contains the initilizations of all the constructed models;
#' \item outs: a list contains the outputs of all the constructed models;
#' }
#'
#' @examples
#' \dontrun{
#' result_CV <- pESCA_CV(dataSets, dataTypes,
#'                             lambdas_CV, penalty='L2', fun_concave='gdp', opts=opts)
#' }
#'
#' @export
pESCA_CV <-function(dataSets, dataTypes,
                            lambdas_CV=NULL, penalty='L2', fun_concave='gdp', opts=list()){
  # check if the inputs statisfy the requirements
  stopifnot(class(dataSets) == "list")
  stopifnot(class(penalty) == "character")
  stopifnot(class(fun_concave) == "character")
  if(length(dataTypes)==1){dataTypes <- unlist(strsplit(dataTypes, split=""))}
  if(exists('quiet', where=opts)){quiet <- opts$quiet} else{quiet<-0};
  if(exists('thr_path', where=opts)){thr_path <- opts$thr_path} else{thr_path<-0};

  # number of data sets, size of each data set
  nTries <- length(lambdas_CV)
  nDataSets <- length(dataSets) # number of data sets
  n <- rep(0,nDataSets)  # number of samples
  d <- rep(0,nDataSets)  # numbers of variables in different data sets
  for(i in 1:nDataSets){
    n[i] <- dim(dataSets[[i]])[1]
    d[i] <- dim(dataSets[[i]])[2]
  }
  if(length(unique(as.factor(n)))!=1)
    stop("multiple data sets have unequal sample size")
  n <- n[1]
  sumd <- sum(d) # total number of variables

  # default dispersion parameters alphas
  if(exists('alphas', where=opts)){alphas<-opts$alphas} else{alphas<-rep(1,nDataSets)};

  # create zero matrix to hold results
  cvErrors_mat <- matrix(data=0, # +1 is used for the sum of all the Xi
                         nrow=nTries, ncol=nDataSets+1)

  # model selection process
  opts_inner <- opts

  # split data sets into training set and test set
  splitedData <- dataSplit(dataSets=dataSets,
                           dataTypes=dataTypes,
                           ratio_mis=0.1)
  trainSets <- splitedData$trainSets

  # save the parameters during the model selection
  inits <- as.list(1:nTries)
  if(thr_path == 1){outs <- as.list(1:nTries)}

  # model selection
  for(j in 1:nTries){
    lambda <- lambdas_CV[j]

    # using the training set to construct a ESCA model
    lambdas <- lambda*rep(1,nDataSets)

    trainModel <- pESCA(dataSets = trainSets,
                              dataTypes = dataTypes,
                              lambdas = lambdas,
                              penalty = penalty,
                              fun_concave= fun_concave,
                              opts=opts_inner)
    if((trainModel$iter <= 2) & (quiet==0)){
      print("less than 3 iteration is used.")
    }

    # warm start
    mu <- trainModel$mu; A <- trainModel$A; B <- trainModel$B
    opts_inner$mu0 <- mu
    opts_inner$A0 <- A
    opts_inner$B0 <- B

    inits[[j]] <- opts_inner
    if(thr_path == 1){outs[[j]] <- trainModel$Sigmas}

    # compute the test error
    ThetaHat <- ones(n) %*% t(mu) + A %*% t(B)
    testError_vec <- cvError_comput(splitedData,dataTypes,alphas,ThetaHat,d)

    cvErrors_tmp <- c(sum(testError_vec), testError_vec)
    cvErrors_mat[j,] <- cvErrors_tmp
  }

  colnames(cvErrors_mat) <- paste0(rep("X_"), c("full", as.character(1:nDataSets)))

  result_CV <- list()
  result_CV$cvErrors_mat <- cvErrors_mat
  result_CV$inits <- inits
  if(thr_path == 1){result_CV$outs <- outs}

  return(result_CV)
}

