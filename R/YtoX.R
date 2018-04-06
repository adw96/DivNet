#' YtoX
#'
#' This function transforms logratio matrix Y to composition X
#'
#' @param Y logratio matrix
#' @param base base value used to calculate logratios
#'
#' @export
YtoX <- function(Y, base) {
    N <- nrow(Y)
    Q <- ncol(Y) + 1
    exp_Y <- exp(Y)
    sum_exp_Y <- apply(exp_Y, 1, sum)
    X <- matrix(0, N, Q)
    for (i in 1:N) {
        X[i, -base] <- exp_Y[i, ]/(sum_exp_Y[i] + 1)
        X[i, base] <- 1/(sum_exp_Y[i] + 1)
    }
    return(X)
}
