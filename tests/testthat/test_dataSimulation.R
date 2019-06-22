## parameter for testing data simulation with group-wise sparsity
n <- 200
ds <- c(400, 200, 100)

# test GGG simulation
test_results_GGG <- test_dataSimu_group(n = n, ds = ds, dataTypes = "GGG")
expect_true(all(unlist(test_results_GGG)))

# test GGB simulation
test_results_GGB <- test_dataSimu_group(n = n, ds = ds, dataTypes = "GGB")
expect_true(all(unlist(test_results_GGB)))

# test GBB simulation
test_results_GBB <- test_dataSimu_group(n = n, ds = ds, dataTypes = "GBB")
expect_true(all(unlist(test_results_GBB)))

# test BBB simulation
test_results_BBB <- test_dataSimu_group(n = n, ds = ds, dataTypes = "BBB")
expect_true(all(unlist(test_results_BBB)))

## testing data simulation only with element-wise sparsity test GGG simulation
test_results_GGG <- test_dataSimu_element(n = n, ds = ds, R = 3, dataTypes = "GGG")
expect_true(all(unlist(test_results_GGG)))

# test GGB simulation
test_results_GGB <- test_dataSimu_element(n = n, ds = ds, R = 3, dataTypes = "GGB")
expect_true(all(unlist(test_results_GGB)))

# test GBB simulation
test_results_GBB <- test_dataSimu_element(n = n, ds = ds, R = 3, dataTypes = "GBB")
expect_true(all(unlist(test_results_GBB)))

# test BBB simulation
test_results_BBB <- test_dataSimu_element(n = n, ds = ds, R = 3, dataTypes = "BBB")
expect_true(all(unlist(test_results_BBB)))








