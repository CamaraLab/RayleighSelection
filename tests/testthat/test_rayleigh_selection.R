#test rayleigh_selection

test_that("R0 and R1 are correct",{
  gy <- nerve_complex(list(c(1,2,6), c(2,3,4), c(4,5,6))) #this simplicial complex is a hollow triangle

  result <- rayleigh_selection(gy,t(as.matrix(c(3, 0 , 3, 0, 0, 3))), num_perms = 100)

  #these can be verified by hand
  expect_equal(result[1,"R0"], 1.5, tolerance = 1e-9)
  expect_equal(result[1,"R1"], 1/6, tolerance = 1e-9)
})
