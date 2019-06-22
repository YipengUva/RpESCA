## test the CV procedure for pESCA model with two tuning parameters
set.seed(123)

# simulate proper data sets
n <- 100
ds <- c(200, 100, 50)
simulatedData <- dataSimu_group_sparse(n = n, ds = ds, 
                                       dataTypes = "GBB", 
                                       noises = rep(1, 3), 
                                       margProb = 0.1, 
                                       sparse_ratio = 0, 
                                       SNRgc = 1, 
                                       SNRlc = rep(1, 3), 
                                       SNRd = rep(1, 3))

# parameters of a pESCA with concave L2norm penalty model
alphas <- rep(1, 3)
gamma <- 1
opts <- list()
opts$gamma <- gamma  # hyper-parameter for the used penalty
opts$rand_start <- 0
opts$tol_obj <- 1e-06  # stopping criteria
opts$maxit <- 500
opts$alphas <- alphas
opts$R <- 30  # components used
opts$thr_path <- 0  # generaint thresholding path or not
opts$quiet <- 1

# specify data sets, data types and used concave functions
dataSets <- simulatedData$X
dataTypes <- "GBB"
fun_concave <- "gdp"

# model selection pESCA conave L2
nTries <- 15
lambdas_g <- log10_seq(from = 1, to = 500, length.out = nTries)
lambdas_b <- log10_seq(from = 5, to = 100, length.out = nTries)
penalty = "L2"

result_CV <- pESCA_CV_twoSteps(dataSets, 
                               dataTypes, 
                               lambdas_g, 
                               lambdas_b, 
                               penalty = penalty, 
                               fun_concave = "gdp", 
                               opts = opts)

# fit the final model
lambdas_opt <- result_CV$lambdas_opt
opts_opt <- result_CV$opts_opt
opts_opt$tol_obj <- 1e-06  # using high precision model

pESCA_L2 <- pESCA(dataSets = dataSets, 
                  dataTypes = dataTypes, 
                  lambdas = lambdas_opt, 
                  penalty = penalty, 
                  fun_concave = fun_concave, 
                  opts = opts_opt)

# test if CV procedure and the final model has been succefully constructed
expect_true(class(pESCA_L2) == "list")
