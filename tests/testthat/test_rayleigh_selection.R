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

test_that("rayleigh_selection works with covariates",{
  # only testing for 0-forms because 1-forms generate an expected warning
  gy <- nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))
  result <- rayleigh_selection(gy, c(0,1,1,0,0,0,0,0,0,1), num_perms = 1000, seed = 10,
                            covariates = c(1,1,0,0,0,0,0,1,1,1), one_forms = FALSE)

  expect_equal(result[1,"R0"], 0.9584665, tolerance = 1e-5)
  expect_equal(result[1,"p0"], 0.383, tolerance = 1e-1)
})

test_that("rayleigh_selection works with an ensamble",{
  gy <- nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))
  result <- rayleigh_selection(list(gy,gy), c(0,1,1,0,0,0,0,0,0,1), num_perms = 1000,
                               seed = 10, one_forms = TRUE)

  expect_equal(result$R0.1[1], 0.9584665, tolerance = 1e-5)
  expect_equal(result$R0.2[1], 0.9584665, tolerance = 1e-5)

  expect_equal(result$R1.1[1], 0.5536398, tolerance = 1e-5)
  expect_equal(result$R1.2[1], 0.5536398, tolerance = 1e-5)

  expect_equal(result$combined.p0, 0.42, tolerance = 1e-1)
  expect_equal(result$combined.p1, 0.148, tolerance = 1e-1)

})

test_that("rayleigh_selection works with an ensamble, covarites and gpd optimization",{
  out <- rayleigh_selection(list(gg, gg), mnist[201:203,], num_perms = 1000,
                            covariates = mnist[200,], seed = 10, num_cores = 4,
                            one_forms = F, optimize.p = "gpd", min_perms = 100)
  out$p0.1 <- NULL
  out$p0.2 <- NULL
  out$q0 <- NULL
  out$n0.conv <- NULL
  out$combined.p0 <- NULL

  expect_snapshot_value(round(out, 6), style = "json2")
})
