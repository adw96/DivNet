test_that("conversions work as expected", {
  
  # to_composition
  yy <- c(2, 1, -1)
  
  for (i in 1:4) {
    # check inverse
    tcm <- to_composition_matrix(yy, base = i)
    expect_equal(sum(tcm), 1)
    expect_equal(yy, to_log_ratios(tcm, base = i))
    
    # check ranks
    expected_ranks <- R.utils::insert(c(4, 3, 1), ats = i, values = 2)
    expect_equal(rank(tcm), expected_ranks)
    
  }
  
  # sanity check
  expect_equal(to_composition_matrix(c(0,0,0), base = i), 
               matrix(rep(0.25, 4), nrow = 1))
  
})
# test error messages 
test_that("errors work as expected", {
  expect_error(to_composition_matrix(Y = list(1)), 
               "to_composition_matrix was passed a list and doesn't know what to do")
  expect_error(to_composition_matrix(Y = array(0, c(2, 2, 2))),
               "to_composition_matrix generated an error")
  expect_error(to_log_ratios(W = NULL, base = NULL),
               "to_log_ratios needs which taxon is base taxon")
})
# test warning message 
test_that("warning works as expected", {
  expect_warning(to_composition_matrix(Y = matrix(0, 5, 5), 
                                       base = NULL),
                 "base is NULL, setting to last taxon")
})
# test log ratio conversion from vector to matrix 
test_that("to_log_ratio converts from vector to matrix", {
  vec <- rexp(10)
  matrix_res <- to_log_ratios(t(as.matrix(vec)), base = 1)
  expect_equal(to_log_ratios(vec, base = 1), matrix_res)
})
