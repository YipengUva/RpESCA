## test the CV procedure for pESCA model when simulated parameters are avaliable
set.seed(123)

# simulate proper data sets
n <- 100
ds <- c(200, 100, 50)
simulatedData <- dataSimu_group_sparse(n = n, ds = ds, 
                                       dataTypes = "GGG", 
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
dataTypes <- "GGG"
fun_concave <- "gdp"


# model selection pESCA conave L2
nTries <- 15
lambdas_CV <- log10_seq(from = 1, to = 500, length.out = nTries)
penalty = "L2"

result_CV_full <- pESCA_CV_fullInfo(simulatedData, 
                                    lambdas_CV = lambdas_CV, 
									penalty = penalty, 
									fun_concave = fun_concave, 
									opts = opts)

# select the model with minimum CV error
index_min_cv <- which.min(result_CV_full$cvErrors_mat[, 1])

# evaluate the performance of the selected model during model selection
expect_true(all(result_CV_full$RVs_mat[index_min_cv, ] > 0.9))
expect_true(all(result_CV_full$RMSEs_mat[index_min_cv, ] < 0.1))








