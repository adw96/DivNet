test_that("make_w returns error", {
  mu <- matrix(0, nrow = 4, ncol = 4)
  Sigma <- matrix(0, nrow = 3, ncol = 3)
  expect_error(make_w(mu, Sigma, mm = 4),
               "mu and Sigma are not the right dimensions")
})

test_that("make_w mm fixed", {
  mu <- matrix(rnorm(16), nrow = 4, ncol = 4)
  Sigma <- matrix(1, nrow = 4, ncol = 4)
  set.seed(1)
  res_vec <- make_w(mu, Sigma, rep(4, nrow(mu)), base = 1)
  set.seed(1)
  expect_equal(make_w(mu, Sigma, mm = 4, base = 1), res_vec)
})
