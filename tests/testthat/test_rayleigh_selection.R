#test rayleigh_selection

test_that("rayleigh_selection works with no optimizations",{
  gy <- nerve_complex(list(c(1,2,6), c(2,3,4), c(4,5,6))) #this simplicial complex is a hollow triangle
  result <- rayleigh_selection(gy,t(as.matrix(c(3, 0 , 3, 0, 0, 3))), num_perms = 1000, one_forms = TRUE)

  #these can be verified by hand
  expect_equal(result[1,"R0"], 1.5, tolerance = 1e-9)
  expect_equal(result[1,"R1"], 1/6, tolerance = 1e-9)
  #checking p-values
  expect_equal(result[1,"p0"], 0.9, tolerance = 1e-1)
  expect_equal(result[1,"p1"], 0.3, tolerance = 1e-1)

  funcs <- rbind(t(as.matrix(c(3, 0 , 3, 0, 0, 3))), t(as.matrix(c(3, 0 , 3, 0, 0, 3))) )
  result_mc <- rayleigh_selection(gy, funcs, num_perms = 1000, one_forms = TRUE, num_cores = 2)
  expect_equal(as.numeric(result_mc[2,]), as.numeric(result[1,]),  tolerance = 2e-1)
})

test_that("rayleigh_selection works with permutation optimization",{
  gy <- nerve_complex(list(c(1,2,6), c(2,3,4), c(4,5,6)))
  result_perm <- rayleigh_selection(gy,t(as.matrix(c(3, 0 , 3, 0, 0, 3))), num_perms = 1000, one_forms = TRUE, optimize.p = "perm")

  expect_equal(result_perm[1,"R0"], 1.5, tolerance = 1e-9)
  expect_equal(result_perm[1,"R1"], 1/6, tolerance = 1e-9)
  expect_equal(result_perm[1,"p0"], 0.9, tolerance = 1e-1)
  expect_equal(result_perm[1,"p1"], 0.3, tolerance = 1e-1)

  funcs <- rbind(t(as.matrix(c(3, 0 , 3, 0, 0, 3))), t(as.matrix(c(3, 0 , 3, 0, 0, 3))) )
  result_perm_mc <- rayleigh_selection(gy, funcs, num_perms = 1000, one_forms = TRUE, num_cores = 2, optimize.p = "perm")
  expect_equal(as.numeric(result_perm_mc[2,]), as.numeric(result_perm[1,]), tolerance = 2e-1)
})

test_that("rayleigh_selection works with GPD optimization",{
  gy <- nerve_complex(list(c(1,2,6), c(2,3,4), c(4,5,6)))
  result_gpd <- rayleigh_selection(gy,t(as.matrix(c(3, 0 , 3, 0, 0, 3))), num_perms = 1000, one_forms = TRUE, optimize.p = "gpd")

  expect_equal(result_gpd[1,"R0"], 1.5, tolerance = 1e-9)
  expect_equal(result_gpd[1,"R1"], 1/6, tolerance = 1e-9)
  expect_equal(result_gpd[1,"p0"], 0.9, tolerance = 1e-1)
  expect_equal(result_gpd[1,"p1"], 0.3, tolerance = 1e-1)


  funcs <- rbind(t(as.matrix(c(3, 0 , 3, 0, 0, 3))), t(as.matrix(c(3, 0 , 3, 0, 0, 3))) )
  result_gpd_mc <- rayleigh_selection(gy, funcs, num_perms = 1000, one_forms = TRUE, num_cores = 2, optimize.p = "gpd")
  expect_equal(as.numeric(result_gpd_mc[2,]), as.numeric(result_gpd[1,]), tolerance = 2e-1)
})
