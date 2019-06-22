## test the CV procedure for pESCA model
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

result_CV <- pESCA_CV(dataSets, dataTypes, 
                      lambdas_CV, 
                      penalty = penalty, 
                      fun_concave = fun_concave, 
                      opts = opts)

# select the model with minimum CV error
index_min_cv <- which.min(result_CV$cvErrors_mat[, 1])

# fit the final model
lambdas_opt <- rep(lambdas_CV[index_min_cv], length(dataSets))
opts_opt <- result_CV$inits[[index_min_cv]]
opts_opt$tol_obj <- 1e-08  # using high precision model

pESCA_L2 <- pESCA(dataSets = dataSets,
                  dataTypes = dataTypes, 
                  lambdas = lambdas_opt, 
                  penalty = penalty, 
                  fun_concave = fun_concave, 
                  opts = opts_opt)

# evaluate the final model
mu <- pESCA_L2$mu
A <- pESCA_L2$A
B <- pESCA_L2$B
S <- pESCA_L2$S
pESCA_L2_eval <- eval_metrics_simu_group(mu, A, B, S, ds, simulatedData)

expect_true(all(pESCA_L2_eval$RVs_structs > 0.95))
expect_true(all(pESCA_L2_eval$RMSEs_params < 0.05))

