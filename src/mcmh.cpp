#include <RcppEigen.h>
#include <Eigen/Core>

using namespace Eigen;

// Like the R function crossprod(m).
MatrixXd
self_crossprod(const Eigen::MatrixXd& m)
{
  const int nrow = m.rows();
  const int ncol = m.cols();

  MatrixXd res = MatrixXd(ncol, ncol).
    setZero().
    selfadjointView<Lower>().
    rankUpdate(m.adjoint());

  return res;
}

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
       const Eigen::MatrixXd& sigInv)
{
  VectorXd vec = Yi_star - eYi;

  // Eq5pt3 <- -0.5 * crossprod((Yi.star - eYi), sigInv) %*% (Yi.star - eYi)
  return -0.5 * ((vec.transpose() * sigInv) * vec)(0,0); // first element!
}

double
Eq5pt4(const Eigen::VectorXd& Yi,
       const Eigen::VectorXd& eYi,
       const Eigen::MatrixXd& sigInv)
{
  VectorXd vec = Yi - eYi;

  // Eq5pt3 <- -0.5 * crossprod((Yi.star - eYi), sigInv) %*% (Yi.star - eYi)
  return -0.5 * ((vec.transpose() * sigInv) * vec)(0,0); // first element!
}

double
get_full_ratio(const Eigen::VectorXd& Wi,
               const Eigen::VectorXd& Wi_no_base,
               const Eigen::VectorXd& Yi,
               const Eigen::VectorXd& Yi_star,
               const Eigen::VectorXd& eYi,
               const Eigen::MatrixXd& sigInv)
{
  return Eq5pt1(Wi, Yi, Yi_star) +
    Eq5pt2(Wi_no_base, Yi, Yi_star) +
    Eq5pt3(Yi_star, eYi, sigInv) -
    Eq5pt4(Yi, eYi, sigInv);
}

// acceptance <- min(1, exp(fullRat))
double
get_acceptance(const double full_ratio)
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
make_yi_star(const double stepsize,
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


//' get_Y_new_and_sigSum
//'
//' Calculate Y_new and sigSum for the current EM iteration.
//'
//' This function avoids having to keep a 3D array in memory and do
//' matrix math along weird strides by updating sigSum after each
//' non-burnt MC iteration, and by keeping a running total of estimated
//' Y values after each non-burnt MC iteration for each sample
//' (converted to a mean at the end of MC iterations).  This also gives
//' the benefit of saving a lot of memory vs. the straight R versions
//' which needed to apply and reduce across a couple different
//' dimensions.  In theory, this should also allow more optimal storage
//' for faster calculations and less memory churn, with the downside
//' being a more complicated function.
//'
//' @author Ryan Moore
//'
//' @param num_samples Number of samples (nrow(Y))
//' @param Y Logratio matrix
//' @param W Corresponding count matrix
//' @param eY Current expected value of logratio matrix
//' @param base OTU index used for base (1-based index)
//' @param sigInv Output of *_network(sigma)
//' @param mciters Number of MC iterations
//' @param iters_to_burn How many MC iters to burn?
//' @param stepsize Variance used for MH samples.
//'   Tweak to adjust acceptance ratio.
//'
//' @return Returns a List containing Y_new and sigSum.
// [[Rcpp::export]]
Rcpp::List
get_Y_new_and_sigSum(const int num_samples,
                     const Eigen::MatrixXd& Y,
                     const Eigen::MatrixXd& W,
                     const Eigen::MatrixXd& eY,
                     const int base,
                     const Eigen::MatrixXd& sigInv,
                     const int mciters,
                     const int iters_to_burn,
                     const double stepsize)
{

  const int notus = W.cols();
  MatrixXd sigSum = MatrixXd(notus - 1, notus - 1).setZero();
  MatrixXd Y_MH_sums = MatrixXd(notus - 1, num_samples).setZero();

  // Track previous lr estimates for all samples.
  MatrixXd previous = MatrixXd(notus - 1, num_samples).setZero();

  for (int mciter = 0; mciter < mciters; ++mciter) {
    // log ratio estimates for this sample.  Normally will have more
    // OTUs than samples.  So save it that way to have longer bits of
    // data in memory.
    MatrixXd sample_lr_estimates = MatrixXd(notus - 1, num_samples);

    for (int smpl = 0; smpl < num_samples; ++smpl) {
      VectorXd Yi = Y.row(smpl);
      VectorXd Wi = W.row(smpl);
      VectorXd eYi = eY.row(smpl);

      if (mciter == 0) {
        previous.col(smpl) = Yi;
      }

      // Remove base from Wi
      VectorXd Wi_no_base = remove_element(base - 1, Wi);

      // Calculate estimated Yi values (Yi_star)
      VectorXd Yi_star = make_yi_star(stepsize, Yi);

      double full_ratio = get_full_ratio(Wi,
                                         Wi_no_base,
                                         Yi,
                                         Yi_star,
                                         eYi,
                                         sigInv);

      double acceptance = get_acceptance(full_ratio);

      // Add estimated Yi vals to the matrix.
      if (R::runif(0, 1) < acceptance) {
        // accepted!
        if (mciter >= iters_to_burn) {
          for (int i = 0; i < notus - 1; ++i) {
            // TODO is it off by 1?
            Y_MH_sums(i, smpl) += Yi_star(i);
            sample_lr_estimates(i, smpl) = Yi_star(i);
          }
        }

        previous.col(smpl) = Yi_star;
      } else { // not accepted!
        if (mciter >= iters_to_burn) {
          for (int i = 0; i < notus - 1; ++i) {
            Y_MH_sums(i, smpl) += previous(i, smpl);
            sample_lr_estimates(i, smpl) = previous(i, smpl);
          }
        }
      }
    }

    if (mciter >= iters_to_burn) {
      // Now we have to update sigSum.  Sadly this will create tmp
      // matrix....
      sigSum += self_crossprod(sample_lr_estimates.transpose() - eY);
    }
  }

  // Finally, need to take mean of estimated lrs across mciters.
  for (int j = 0; j < Y_MH_sums.cols(); ++j) {
    for (int i = 0; i < Y_MH_sums.rows(); ++i) {
      Y_MH_sums(i, j) /= (mciters - iters_to_burn);
    }
  }

  MatrixXd Y_new = Y_MH_sums.transpose();

  return Rcpp::List::create(Rcpp::Named("Y_new") = Rcpp::wrap(Y_new),
                            Rcpp::Named("sigSum") = Rcpp::wrap(sigSum));
}
