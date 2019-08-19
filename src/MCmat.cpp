#include <RcppEigen.h>
#include <random>


using namespace std;
using namespace Rcpp;
using namespace Eigen;

typedef Eigen::Map<Eigen::MatrixXd> MappedMatrixXd;

// nrow logratios = no. samples
// ncol logratios = no. otus - 1

// nrow counts = no. samples
// ncol couts = no. otus

// Return value is a list of mc_iters x notus matrices

// [[Rcpp::export]]
Rcpp::List
eigen_mc_array(const Rcpp::NumericMatrix r_logratios,
               const Rcpp::NumericMatrix r_counts,
               const Rcpp::NumericMatrix r_expected_logratios,
               const int base_otu,
               const Rcpp::NumericMatrix r_sigma_inverse,
               const int mc_iters,
               const int iters_to_burn,
               const double stepsize)
{

  MappedMatrixXd logratios(as<MappedMatrixXd>(r_logratios));
  MappedMatrixXd counts(as<MappedMatrixXd>(r_counts));
  MappedMatrixXd expected_logratios(as<MappedMatrixXd>(r_expected_logratios));
  MappedMatrixXd sigma_inverse(as<MappedMatrixXd>(r_sigma_inverse));

  // Original function, N is num samples aka counts.rows()
  // Original function, Q is num OTUs aka counts.cols()

  const int nsamples = counts.rows();
  const int notus = counts.cols();

  int mci = 0;
  std::vector<NumericMatrix> mc;
  vector<MatrixXd> per_iter_lrs;
  for (int i = 0; i < mc_iters; ++i) {
    per_iter_lrs.push_back(MatrixXd(notus - 1, nsamples).setZero());
  }


  for (int sidx = 0; sidx < nsamples; ++sidx) {
    const Eigen::VectorXd Yi = logratios.row(sidx);
    const Eigen::VectorXd Wi = counts.row(sidx);
    const Eigen::VectorXd eYi = expected_logratios.row(sidx);

    // Wi minus the base OTU
    Eigen::VectorXd Wi_no_base(Wi.size() - 1);
    int cur_idx = 0;
    for (int i = 0; i < Wi.size(); ++i) {
      if (i != (base_otu - 1)) {
        Wi_no_base(cur_idx++) = Wi(i);
      }
    }

    NumericMatrix Yi_MH(mc_iters, notus);

    Eigen::VectorXd Yi_star(Yi.size());

    double Eq5pt1 = 0.0,
      Eq5pt2 = 0.0,
      Eq5pt3 = 0.0,
      Eq5pt4 = 0.0,
      fullRat = 0.0,
      fullRat_exp = 0.0,
      acceptance = 0.0;

    // The main MCiter loop
    for (int iter = 0; iter < mc_iters; ++iter) {
      if (iter == 0) {
        // The first iteration, take the original Yi values.
        for (int i = 0; i < Yi.size(); ++i) {
          Yi_star(i) = Yi(i) + R::rnorm(0, stepsize);
        }
      } else {
        // For subsequent iterations, pull Yi_MH values from the last
        // MC iteration.
        for (int i = 1; i < notus; ++i) {
          Yi_star(i - 1) = Yi_MH(iter - 1, i) + R::rnorm(0, stepsize); // first column is the acceptance
        }
      }

      // Denominator
      Eq5pt1 = Wi.sum() * (log(Yi     .array().exp().sum() + 1) -
                           log(Yi_star.array().exp().sum() + 1));

      // Numerator
      Eq5pt2 = (Wi_no_base * (Yi_star - Yi)).sum();

      Eq5pt3 = -0.5 * ((Yi_star - eYi).transpose() * sigma_inverse) * (Yi_star - eYi);

      Eq5pt4 = -0.5 * ((Yi      - eYi).transpose() * sigma_inverse) * (Yi      - eYi);

      fullRat = Eq5pt1 + Eq5pt2 + Eq5pt3 - Eq5pt4;

      fullRat_exp = exp(fullRat);

      if (fullRat_exp < 1) {
        acceptance = fullRat_exp;
      } else {
        acceptance = 1;
      }

      if (isnan(acceptance) || R::runif(0, 1) < acceptance) {
        Yi_MH(iter, 0) = 1;
        for (int j = 1; j < notus; ++j) {
          Yi_MH(iter, j) = Yi_star(j - 1);  // Yi_star has otus-1 items
        }
      } else {
        Yi_MH(iter, 0) = 0;
        for (int j = 1; j < notus; ++j) {
          Yi_MH(iter, j) = Yi(j - 1); // Yi_star has otus-1 items
        }
      }

      // Add the logratios to the other data frame.
      MatrixXd & these_lrs = per_iter_lrs[iter];
      for (int j = 1; j < notus - 1; ++j) {
        these_lrs(j - 1, sidx) = Yi_MH(iter, j);
      }
    }

    mc.push_back(Yi_MH);
  }

  MatrixXd sigma = MatrixXd(notus - 1, notus - 1).setZero();
  for (int iter = iters_to_burn; iter < mc_iters; ++iter) { // per_iter_lrs
    MatrixXd & these_lrs = per_iter_lrs[iter];

    MatrixXd tmp = these_lrs.transpose() - expected_logratios;
    MatrixXd cprod = tmp.transpose() * tmp;

    sigma += cprod;
  }

  sigma /= (nsamples * (mc_iters - iters_to_burn));

  return List::create(Named("mc") = wrap(mc),
                      Named("sigma") = wrap(sigma));
}

// [[Rcpp::export]]
std::vector<Eigen::MatrixXd>
thing()
{
  int mc_iters = 10;
  int notus = 5;
  int nsamples = 2;
  int n = 0;

  vector<MatrixXd> per_iter_lrs;
  for (int i = 0; i < mc_iters; ++i) {
    per_iter_lrs.push_back(MatrixXd(nsamples, notus - 1).setZero());
  }

  for (int i = 0; i < mc_iters; ++i) {
    MatrixXd & tmp = per_iter_lrs[i];

    for (int sidx = 0; sidx < nsamples; ++sidx) {
      for (int oidx = 0; oidx < notus - 1; ++oidx) {
        tmp(sidx, oidx) = n++;
      }
    }
  }

  return per_iter_lrs;
}
