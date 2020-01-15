#include <RcppEigen.h>
#include <Eigen/Core>

using namespace Eigen;

// Convenience type for mapped MatrixXd
typedef Eigen::Map<Eigen::MatrixXd> MappedMatrixXd;
typedef Eigen::Map<Eigen::VectorXd> MappedVectorXd;

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

// [[Rcpp::export]]
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

// TODO exp can fail!
// acceptance <- min(1, exp(fullRat))
double
mcrow_acceptance(const double full_ratio)
{
  double full_ratio_exp = exp(full_ratio);

  if (full_ratio_exp < 1) {
    return full_ratio_exp;
  } else {
    return 1;
  }
}

// [[Rcpp::export]]
Eigen::MatrixXd
mcrow_mc_iteration(const int num_iters,
                   const Eigen::VectorXd& unif_rand_vals,
                   const Eigen::VectorXd& Wi,
                   const Eigen::VectorXd& Wi_no_base,
                   const Eigen::VectorXd& Yi,
                   const Eigen::MatrixXd& Yi_star_all,
                   const Eigen::VectorXd& eYi,
                   const Eigen::MatrixXd& sigma_inverse)
{

  // Note that Wi length is the number of OTUs.
  const int notus = Wi.size();

  // # extra column for acceptance indicator
  // Yi.MH <- matrix(0, MCiter, Q)
  MatrixXd Yi_MH = MatrixXd(num_iters, notus).setZero();

  for (int i = 0; i < num_iters; ++i) {
    VectorXd Yi_star = Yi_star_all.col(i);

    double full_ratio = mcrow_full_ratio(Wi,
                                         Wi_no_base,
                                         Yi,
                                         Yi_star,
                                         eYi,
                                         sigma_inverse);

    double acceptance = mcrow_acceptance(full_ratio);

    // Note: The original code checked for is.nan for acceptance.....
    if (unif_rand_vals(i) < acceptance) {
      Yi_MH(i, 0) = 1; // accepted!

      // TODO switch to Yi_MH.ncols
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
