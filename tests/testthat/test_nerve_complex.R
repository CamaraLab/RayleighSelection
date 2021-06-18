
#Tests nerve complex


test_that("Correct complex in built", {
  gy <- nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))

  expect_equal(gy$points_in_vertex, list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))

  adj <- c(0,0,0,1,0,
           1,0,1,0,0,
           0,0,0,1,0,
           0,0,0,0,0,
           1,0,0,1,0)
  expect_equal(as.matrix(gy$adjacency), matrix(adj, nrow = 5, ncol = 5, byrow = TRUE), ignore_attr = TRUE)

  expect_true(all(as.logical(lapply(gy$two_simplices, function(x) all(x>= 0) ) )))
  expect_equal(lapply(lapply(gy$two_simplices, sum), sum), list(1, 0,0,0,0))
  expect_equal(gy$two_simplices[[1]][4,5], 1)

  expect_equal(gy$order, c(2,3,5,1,4))
})
