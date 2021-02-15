#test rayleigh_selection

test_that("R0 and R1 are correct",{
  gy <- nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))

  result <- rayleigh_selection(gy,t(as.data.frame(c(0,1,1,0,0,0,0,0,0,1))), num_perms = 100)

  expect_true(is.data.frame(result))
  expect_equal(result[1,"R0"], 0.958466, tolerance = 10**(-6))
  expect_equal(result[1,"R1"], 0.553640, tolerance = 10**(-6))

})
