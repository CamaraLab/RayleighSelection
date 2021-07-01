#test rayleigh_selection
load("../testdata/mnist_complex.RData")

test_that("rayleigh_selection works with no optimizations",{
  gy <- nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))
  result <- rayleigh_selection(gy,c(0,1,1,0,0,0,0,0,0,1), num_perms = 1000, one_forms = TRUE)

  expect_equal(result[1,"R0"], 0.9584665, tolerance = 1e-5)
  expect_equal(result[1,"R1"], 0.5536398, tolerance = 1e-5)

  expect_equal(result[1,"p0"], 0.42, tolerance = 1e-1)
  expect_equal(result[1,"p1"], 0.148, tolerance = 1e-1)

  funcs <- matrix(rep(c(0,1,1,0,0,0,0,0,0,1), 2), nrow = 2, byrow = T)
  result_mc <- rayleigh_selection(gy, funcs, num_perms = 1000, one_forms = TRUE, num_cores = 2)
  expect_equal(as.numeric(result_mc[2,]), as.numeric(result[1,]),  tolerance = 2e-1)
})

test_that("rayleigh_selection works with permutation optimization",{
  out <- rayleigh_selection(gg, mnist[201:203,], num_perms = 400, seed = 10,
                            one_forms = TRUE, optimize.p = "perm", min_perms = 50)
  expect_identical(ncol(out), 8L)
  out$n0.conv <- NULL
  out$n1.conv <- NULL
  rand.cols <- c("p0", "q0", "p1", "q1")
  expect_true(all(out[, names(out) %in% rand.cols] <= 1  &  out[, names(out) %in% rand.cols] >=0))
  expect_snapshot_value(out[, !(names(out) %in% rand.cols)], style = "serialize")
})

test_that("rayleigh_selection works with GPD optimization",{
  out <- rayleigh_selection(gg, mnist[201:203,], num_perms = 400, seed = 10,
                            one_forms = TRUE, optimize.p = "gpd", min_perms = 50)
  expect_identical(ncol(out), 8L)
  out$n0.conv <- NULL
  out$n1.conv <- NULL
  rand.cols <- c("p0", "q0", "p1", "q1")
  expect_true(all(out[, names(out) %in% rand.cols] <= 1  &  out[, names(out) %in% rand.cols] >=0))
  expect_snapshot_value(out[, !(names(out) %in% rand.cols)], style = "serialize")
})

test_that("rayleigh_selection works with covariates",{
  gy <- nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))
  result <- rayleigh_selection(gy, c(0,1,1,0,0,0,0,0,0,1), num_perms = 1000, seed = 10,
                            covariates = c(1,1,0,0,0,0,0,1,1,1), one_forms = TRUE)

  expect_equal(result[1,"R0"], 0.9584665, tolerance = 1e-5)
  expect_equal(result[1,"R1"], 0.5536398, tolerance = 1e-5)

  expect_equal(result[1,"p0"], 0.383, tolerance = 1e-1)
  expect_equal(result[1,"p1"], 0.153, tolerance = 1e-1)
})

test_that("rayleigh_selection works with an ensamble",{
  gy <- nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))
  result <- rayleigh_selection(list(gy,gy), c(0,1,1,0,0,0,0,0,0,1), num_perms = 1000,
                               seed = 10, one_forms = TRUE)

  expect_equal(result$R0[1,1], 0.9584665, tolerance = 1e-5)
  expect_equal(result$R0[1,2], 0.9584665, tolerance = 1e-5)

  expect_equal(result$R1[1,1], 0.5536398, tolerance = 1e-5)
  expect_equal(result$R1[1,2], 0.5536398, tolerance = 1e-5)

  expect_equal(result$p0, 0.42, tolerance = 1e-1)
  expect_equal(result$p1, 0.148, tolerance = 1e-1)

})
