#include <RcppEigen.h>
#include <Eigen/Core>

using namespace Eigen;

double
Eq5pt1(const Eigen::VectorXd& Wi,
       const Eigen::VectorXd& Yi,
       const Eigen::VectorXd& Yi_star)
{

  // Eq5pt1 <- sum(Wi) * (log(sum(exp(Yi)) + 1) - log(sum(exp(Yi.star)) + 1))
  return Wi.sum() * (log(Yi.array().exp().sum() + 1) -
                     log(Yi_star.array().exp().sum() + 1));
}

double
Eq5pt2(const Eigen::VectorXd& Wi_no_base,
       const Eigen::VectorXd& Yi,
       const Eigen::VectorXd& Yi_star)
{

  // Eq5pt2 <- sum(Wi[-base] * (Yi.star - Yi))
  return (Wi_no_base * (Yi_star - Yi)).sum();
}


double
Eq5pt3(const Eigen::VectorXd& Yi_star,
       const Eigen::VectorXd& eYi,
       const Eigen::MatrixXd& sigma_inverse)
{
  VectorXd vec = Yi_star - eYi;

  // Eq5pt3 <- -0.5 * crossprod((Yi.star - eYi), sigInv) %*% (Yi.star - eYi)
  return -0.5 * ((vec.transpose() * sigma_inverse) * vec)(0,0); // first element!
}

double
Eq5pt4(const Eigen::VectorXd& Yi,
       const Eigen::VectorXd& eYi,
       const Eigen::MatrixXd& sigma_inverse)
{
  VectorXd vec = Yi - eYi;

  // Eq5pt3 <- -0.5 * crossprod((Yi.star - eYi), sigInv) %*% (Yi.star - eYi)
  return -0.5 * ((vec.transpose() * sigma_inverse) * vec)(0,0); // first element!
}

double
mcrow_full_ratio(const Eigen::VectorXd& Wi,
                 const Eigen::VectorXd& Wi_no_base,
                 const Eigen::VectorXd& Yi,
                 const Eigen::VectorXd& Yi_star,
                 const Eigen::VectorXd& eYi,
                 const Eigen::MatrixXd& sigma_inverse)
{
  return Eq5pt1(Wi, Yi, Yi_star) +
    Eq5pt2(Wi_no_base, Yi, Yi_star) +
    Eq5pt3(Yi_star, eYi, sigma_inverse) +
    Eq5pt4(Yi, eYi, sigma_inverse);
}

// acceptance <- min(1, exp(fullRat))
double
mcrow_acceptance(const double full_ratio)
{
  // TODO exp can fail!
  double full_ratio_exp = exp(full_ratio);

  if (full_ratio_exp < 1) {
    return full_ratio_exp;
  } else {
    return 1;
  }
}

Eigen::VectorXd
mcrow_make_yi_star(const double stepsize,
                   const Eigen::VectorXd& Yi)
{
  VectorXd Yi_star(Yi.size());

  // This class calls GetRNGState in it's constructor and PutRNGState
  // in its destructure, which update .Random.seed properly.  Since we
  // don't export this function, we must do this manually.
  Rcpp::RNGScope rcpp_rngScope_gen;

  for (int i = 0; i < Yi_star.size(); ++i) {
    Yi_star(i) = Yi(i) + R::rnorm(0, stepsize);
  }

  return Yi_star;
}

Eigen::MatrixXd
mcrow_mc_iteration(const int num_iters,
                   const double stepsize,
                   const Eigen::VectorXd& Wi,
                   const Eigen::VectorXd& Wi_no_base,
                   const Eigen::VectorXd& Yi,
                   const Eigen::VectorXd& eYi,
                   const Eigen::MatrixXd& sigma_inverse)
{

  // Note that Wi length is the number of OTUs.
  const int notus = Wi.size();

  // # extra column for acceptance indicator
  // Yi.MH <- matrix(0, MCiter, Q)
  MatrixXd Yi_MH = MatrixXd(num_iters, notus).setZero();

  // Set up the random number gen
  Rcpp::RNGScope rcpp_rngScope_gen;

  for (int i = 0; i < num_iters; ++i) {
    VectorXd Yi_star = mcrow_make_yi_star(stepsize, Yi);

    double full_ratio = mcrow_full_ratio(Wi,
                                         Wi_no_base,
                                         Yi,
                                         Yi_star,
                                         eYi,
                                         sigma_inverse);

    double acceptance = mcrow_acceptance(full_ratio);

    // Note: The original code checked for is.nan for acceptance.....
    if (R::runif(0, 1) < acceptance) {
      Yi_MH(i, 0) = 1; // accepted!

      for (int j = 1; j < notus; ++j) {
        Yi_MH(i, j) = Yi_star(j - 1); // Yi_star has one fewer item!
      }
    } else {
      Yi_MH(i, 0) = 0; // not accepted

      // If we're on the first iteration, use the original Yi vals,
      // else use the last iterations.
      if (i == 0) {
        for (int j = 1; j < notus; ++j) {
          Yi_MH(i, j) = Yi(j - 1); // Yi has one fewer item!
        }
      } else {
        for (int j = 1; j < notus; ++j) {
          Yi_MH(i, j) = Yi_MH(i - 1, j);
        }
      }
    }
  }

  return Yi_MH;
}

// This is a zero based index!!!!!!!!!
Eigen::VectorXd
remove_element(const int idx_to_remove,
               const Eigen::VectorXd v)
{
  VectorXd new_v(v.size() - 1);

  int new_idx = 0;
  for (int i = 0; i < v.size(); ++i) {
    if (i != idx_to_remove) {
      new_v(new_idx++) = v(i);
    }
  }

  return new_v;
}

//' MCrow
//'
//' This function simulates MC step for a single row.
//'
//' @author Ryan Moore
//'
//' @param Yi row of logratio matrix
//' @param Wi corresponding row of count matrix
//' @param eYi current expected value of logratio matrix
//' @param base OTU index used for base (1-based index)
//' @param sigInv current estimate of sigma inverse
//' @param MCiter number of MC samples to generate
//' @param stepsize variance used for MH samples. Tweak to adjust acceptance ratio
// [[Rcpp::export]]
Eigen::MatrixXd
mcrow_MCrow(const Eigen::VectorXd& Yi,
            const Eigen::VectorXd& Wi,
            const Eigen::VectorXd& eYi,
            const int base,
            const Eigen::MatrixXd& sigInv,
            const int MCiter,
            const double stepsize)
{
  // Remove base from Wi
  VectorXd Wi_no_base = remove_element(base - 1, Wi);

  MatrixXd Yi_MH = mcrow_mc_iteration(MCiter,
                                      stepsize,
                                      Wi,
                                      Wi_no_base,
                                      Yi,
                                      eYi,
                                      sigInv);

  return Yi_MH;
}
