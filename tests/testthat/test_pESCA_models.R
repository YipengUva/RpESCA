## test pESCA models
set.seed(123)

# simulate proper data sets
n <- 200
ds <- c(400, 200, 100)
simulatedData <- dataSimu_group_sparse(n = n, ds = ds, dataTypes = "GGG", noises = rep(1, 3), margProb = 0.1, 
    sparse_ratio = 0, SNRgc = 1, SNRlc = rep(1, 3), SNRd = rep(1, 3))

# parameters of a pESCA with concave L2norm penalty model
alphas <- rep(1, 3)
gamma <- 1
opts <- list()
opts$gamma <- gamma  # hyper-parameter for the used penalty
opts$rand_start <- 0
opts$tol_obj <- 1e-06  # stopping criteria
opts$maxit <- 3
opts$alphas <- alphas
opts$R <- 30  # components used
opts$thr_path <- 0  # generaint thresholding path or not
opts$quiet <- 1

# set up the pESCA model
dataSets <- simulatedData$X
dataTypes <- "GGG"  # or   rep('G',3)
lambdas <- rep(30, 3)
fun_concave <- "gdp"

## pESCA with concave L2 norm penalty
pESCA_L2 <- pESCA(dataSets = dataSets, 
                  dataTypes = dataTypes, 
                  lambdas = lambdas, 
                  penalty = "L2", 
                  opts = opts)
# test if pESCA_L2 model has been succefully constructed
expect_true(class(pESCA_L2) == "list")


## pESCA with concave L1 norm penalty
lambdas <- rep(10, 3)
pESCA_L1 <- pESCA(dataSets = dataSets, 
                  dataTypes = dataTypes, 
                  lambdas = lambdas,
                  penalty = "L1", 
                  fun_concave = fun_concave, 
                  opts = opts)
# test if pESCA_L1 model has been succefully constructed
expect_true(class(pESCA_L1) == "list")


## pESCA with composite concave penalty
lambdas <- rep(10, 3)
penalty <- "composite"

pESCA_composite <- pESCA(dataSets = dataSets, 
                         dataTypes = dataTypes, 
                         lambdas = lambdas, 
                         penalty = penalty, 
                         fun_concave = fun_concave, 
                         opts = opts)
# test if pESCA_composite model has been succefully constructed
expect_true(class(pESCA_composite) == "list")

## pESCA with concave element-wise penalty
lambdas <- rep(5, 3)
penalty <- "element"
opts$R <- 10

pESCA_element <- pESCA(dataSets = dataSets, 
                       dataTypes = dataTypes, 
                       lambdas = lambdas, 
                       penalty = penalty, 
                       fun_concave = fun_concave, 
                       opts = opts)
# test if pESCA_element model has been succefully constructed
expect_true(class(pESCA_element) == "list")

