## test penalty functions we mainly test the propertiy of a conave penalty function the function is
## monotonically non-decreasing the super-gradient function is monotonically non-increasing
x <- 0:100
lambda <- 20

# GDP and its super-gradient
penalty_name <- "gdp"
penalty_fun <- get(penalty_name)
penalty_fun_sg <- get(paste0(penalty_name, "_sg"))

x_gdp <- penalty_fun(x = x, gamma = 1, lambda = lambda)
x_gdp_sg <- penalty_fun_sg(x = x, gamma = 1, lambda = lambda)

test_eixst_gdp <- all(c(class(penalty_fun), class(penalty_fun_sg)) == rep("function", 2))
test_increa_gdp <- all(x_gdp == sort(x_gdp, decreasing = FALSE))
test_decrea_gdp <- all(x_gdp_sg == sort(x_gdp_sg, decreasing = TRUE))

expect_true(test_eixst_gdp & test_increa_gdp & test_decrea_gdp)

# Lq and its super-gradient
penalty_name <- "lq"
penalty_fun <- get(penalty_name)
penalty_fun_sg <- get(paste0(penalty_name, "_sg"))

x_lq <- penalty_fun(x = x, gamma = 1, lambda = lambda)
x_lq_sg <- penalty_fun_sg(x = x, gamma = 1, lambda = lambda)

test_eixst_lq <- all(c(class(penalty_fun), class(penalty_fun_sg)) == rep("function", 2))
test_increa_lq <- all(x_lq == sort(x_lq, decreasing = FALSE))
test_decrea_lq <- all(x_lq_sg == sort(x_lq_sg, decreasing = TRUE))

expect_true(test_eixst_lq & test_increa_lq & test_decrea_lq)

# SCAD and its super-gradient
penalty_name <- "scad"
penalty_fun <- get(penalty_name)
penalty_fun_sg <- get(paste0(penalty_name, "_sg"))

x_scad <- penalty_fun(x = x, gamma = 1, lambda = lambda)
x_scad_sg <- penalty_fun_sg(x = x, gamma = 1, lambda = lambda)

test_eixst_scad <- all(c(class(penalty_fun), class(penalty_fun_sg)) == rep("function", 2))
test_increa_scad <- all(x_scad == sort(x_scad, decreasing = FALSE))
test_decrea_scad <- all(x_scad_sg == sort(x_scad_sg, decreasing = TRUE))

expect_true(test_eixst_scad & test_increa_scad & test_decrea_scad)


