#' YtoW
#'
#' This function transforms logratio matrix Y to counts W
#'
#' @param Y logratio matrix
#' @param M vector of observed counts by row
#' @param base base value used to calculate logratios
#'
#' @export
YtoW <- function(Y, M, base) {
    N <- nrow(Y)
    Q <- ncol(Y) + 1
    exp_Y <- exp(Y)
    sum_exp_Y <- apply(exp_Y, 1, sum)
    X <- W <- matrix(0, N, Q)
    for (i in 1:N) {
        X[i, -base] <- exp_Y[i, ]/(sum_exp_Y[i] + 1)
        X[i, base] <- 1/(sum_exp_Y[i] + 1)
        # W[i, ] <- M[i] * X[i, ]
        W[i, ] <- rmultinom(n = 1, size = M[i], prob = X[i, ])
    }
    return(W)
}
