## test logit function
prob_mat <- matrix(rbeta(3 * 4, 1, 1), nrow = 3, ncol = 4)
expect_true(all(log(prob_mat/(1 - prob_mat)) == logit(prob_mat)))

## test inverse logit function
real_mat <- matrix(rnorm(3 * 4), nrow = 3, ncol = 4)
expect_true(all((1/(1 + exp(-real_mat))) == inver_logit(real_mat)))

## test index genrating function
expect_true(all(1:400 == index_Xi(1, c(400, 200, 100))) & 
               all(401:600 == index_Xi(2, c(400, 200, 100))) & 
               all(601:700 == index_Xi(3, c(400, 200, 100))))

## test varExp_Gaussian function fix a data set into the package and use it to test the function

## test log10_seq function
from <- 1
to <- 500
length.out <- 30
expect_true(all(10^(seq(log10(from), log10(to), length.out = length.out)) == log10_seq(from, to, length.out)))

## test obj_logistic function Theta <- matrix(rnorm(3*4),3,4) X <- matrix(data=0,3,4) X[Theta>0] <- 1
## obj_logistic(X,Theta)




