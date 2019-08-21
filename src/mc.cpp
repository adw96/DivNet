#include <RcppEigen.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

using namespace std;
using namespace Rcpp;
using namespace Eigen;
using namespace RcppParallel;

// Convenience type for mapped MatrixXd
typedef Eigen::Map<Eigen::MatrixXd> MappedMatrixXd;

//' remove_base_otu
//'
//' Returns a new VectorXd with the `base_otu` removed.
//'
//' @author Ryan Moore
//'
//' @param base_otu The index of the otu to remove from the input vector
//' @param Wi A single row/sample of the count matrix
//'
//' @return A count vector minus the `base_otu`
VectorXd
remove_base_otu(const int base_otu,
                const Eigen::Ref<Eigen::VectorXd> Wi)
{
  VectorXd Wi_no_base(Wi.size() - 1);

  int current_idx = 0;
  for (int i = 0; i < Wi.size(); ++i) {
    if (i != (base_otu - 1)) {
      Wi_no_base(current_idx++) = Wi(i);
    }
  }

  return Wi_no_base;
}

//' mc_loop
//'
//' Runs the MC step for a single sample.
//'
//' @author Ryan Moore
//'
//' @param sidx Index of the current sample
//' @param mc_iters Total number of mc iters to run
//' @param stepsize Variance used for MH samples
//' @param Yi A single row/sample of the logratio matrix
//' @param Wi A single row/sample of the count matrix (same row/sample as for Yi)
//' @param Wi_no_base The same row/sample of the count matrix as Wi, but with the base OTU removed
//' @param eYi A single row/sample of the expected value of the logratio matrix
//' @param sigma_inverse Current estimate of sigma inverse
//' @param Yi_MH An mc_iters x Wi.size() (no. otus) matrix to hold the MC iteration results.  This param will be modified in place.
//' @param per_iter_lrs A vector that stores an (no. otus - 1 ie Yi.size()) x no. samples matrix for each MC iteration.  This param will be modified in place.
void
mc_loop(const int sidx,
        const int mc_iters,
        const int stepsize,
        const Eigen::Ref<Eigen::VectorXd> Yi,
        const Eigen::Ref<Eigen::VectorXd> Wi,
        const Eigen::Ref<Eigen::VectorXd> Wi_no_base,
        const Eigen::Ref<Eigen::VectorXd> eYi,
        const Eigen::Ref<Eigen::MatrixXd> sigma_inverse,
        Eigen::Ref<Eigen::MatrixXd> Yi_MH,
        vector<MatrixXd> per_iter_lrs)
{
  const int notus = Wi.size();

  VectorXd Yi_star(Yi.size());

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
    MatrixXd & lrs = per_iter_lrs[iter];
    for (int j = 1; j < notus - 1; ++j) {
      lrs(j - 1, sidx) = Yi_MH(iter, j);
    }
  }
}


//' run_mc
//'
//' Simulate MC step.
//'
//' @author Ryan Moore
//'
//' @param logratios The logratio matrix.  Dim: nsamples x (notus - 1).
//' @param counts The corresponding count matrix.  Dim: nsamples x notus.
//' @param expected_logratios The current expected value of the logratio matrix.  Dim: nsamples x (notus - 1).
//' @param sigma_inverse The current estimate of sigma inverse.
//' @param base_otu The index of the otu to remove from the input vector
//' @param mc_iters Total number of mc iters to run
//' @param iters_to_burn The number of MC iters to burn
//' @param stepsize Variance used for MH samples
// [[Rcpp::export]]
Rcpp::List
run_mc(const Rcpp::NumericMatrix logratios,
       const Rcpp::NumericMatrix counts,
       const Rcpp::NumericMatrix expected_logratios,
       const Rcpp::NumericMatrix sigma_inverse,
       const int base_otu,
       const int mc_iters,
       const int iters_to_burn,
       const double stepsize)
{

  MappedMatrixXd mlogratios(as<MappedMatrixXd>(logratios));
  MappedMatrixXd mcounts(as<MappedMatrixXd>(counts));
  MappedMatrixXd mexpected_logratios(as<MappedMatrixXd>(expected_logratios));
  MappedMatrixXd msigma_inverse(as<MappedMatrixXd>(sigma_inverse));

  // Original function, N is num samples aka counts.rows()
  // Original function, Q is num OTUs aka counts.cols()

  const int nsamples = mcounts.rows();
  const int notus = mcounts.cols();

  int mci = 0;
  vector<MatrixXd> mc;
  vector<MatrixXd> per_iter_lrs;
  for (int i = 0; i < mc_iters; ++i) {
    per_iter_lrs.push_back(MatrixXd(notus - 1, nsamples).setZero());
  }


  for (int sidx = 0; sidx < nsamples; ++sidx) {
    VectorXd Yi = mlogratios.row(sidx);
    VectorXd Wi = mcounts.row(sidx);
    VectorXd eYi = mexpected_logratios.row(sidx);

    // Wi minus the base OTU
    VectorXd Wi_no_base = remove_base_otu(base_otu, Wi);

    MatrixXd Yi_MH = MatrixXd(mc_iters, notus);

    // Modifies Yi_MH and per_iter_lrs.
    mc_loop(sidx,
            mc_iters,
            stepsize,
            Yi,
            Wi,
            Wi_no_base,
            eYi,
            msigma_inverse,
            Yi_MH,
            per_iter_lrs);

    mc.push_back(Yi_MH);
  }

  MatrixXd sigma = MatrixXd(notus - 1, notus - 1).setZero();
  for (int iter = iters_to_burn; iter < mc_iters; ++iter) { // per_iter_lrs
    MatrixXd & lrs = per_iter_lrs[iter];

    MatrixXd tmp = lrs.transpose() - mexpected_logratios;
    MatrixXd cprod = tmp.transpose() * tmp;

    sigma += cprod;
  }

  sigma /= (nsamples * (mc_iters - iters_to_burn));

  return List::create(Named("mc") = wrap(mc),
                      Named("sigma") = wrap(sigma));
}
