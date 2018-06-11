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
