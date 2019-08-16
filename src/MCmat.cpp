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

// [[Rcpp::export]]
std::vector<Rcpp::NumericMatrix>
eigen_mc_array(const Rcpp::NumericMatrix r_logratios,
                const Rcpp::NumericMatrix r_counts,
                const Rcpp::NumericMatrix r_expected_logratios,
                const int base_otu,
                const Rcpp::NumericMatrix r_sigma_inverse,
                const int mc_iters,
                const double stepsize)
{

  MappedMatrixXd logratios(as<MappedMatrixXd>(r_logratios));
  MappedMatrixXd counts(as<MappedMatrixXd>(r_counts));
  MappedMatrixXd expected_logratios(as<MappedMatrixXd>(r_expected_logratios));
  MappedMatrixXd sigma_inverse(as<MappedMatrixXd>(r_sigma_inverse));

  default_random_engine generator;

  normal_distribution<double> norm_dist(0.0, stepsize);
  // This is the interval [0.0, 1.0)
  uniform_real_distribution<double> unif_dist(0.0, 1.0);

  // Original function, N is num samples aka counts.rows()
  // Original function, Q is num OTUs aka counts.cols()

  const int nsamples = counts.rows();
  const int notus = counts.cols();

  int mci = 0;
  std::vector<NumericMatrix> mc;

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
      for (int i = 0; i < Yi.size(); ++i) {
        Yi_star(i) = Yi(i) + norm_dist(generator);
      }

      // Denominator
      Eq5pt1 = Wi.sum() * (log(Yi     .array().exp().sum() + 1) -
                           log(Yi_star.array().exp().sum() + 1));

      // Numerator
      Eq5pt2 = (Wi_no_base * (Yi_star - Yi)).sum();

      // If anything any sigma_inverse is NaN, then Eq5pt3 and 4 will
      // be NaN.  Then fullRat will be too.  If this happens a lot we
      // might save time by checking for NaN before doing the actual
      // matrix multiplication.

      Eq5pt3 = -0.5 * ((Yi_star - eYi).transpose() * sigma_inverse) * (Yi_star - eYi);

      Eq5pt4 = -0.5 * ((Yi      - eYi).transpose() * sigma_inverse) * (Yi      - eYi);

      fullRat = Eq5pt1 + Eq5pt2 + Eq5pt3 - Eq5pt4;

      fullRat_exp = exp(fullRat);

      if (fullRat_exp < 1) {
        acceptance = fullRat_exp;
      } else {
        acceptance = 1;
      }

      if (isnan(acceptance) || unif_dist(generator) < acceptance) {
        Yi_MH(iter, 0) = 1;
        for (int j = 1; j < Yi.size(); ++j) {
          Yi_MH(iter, j) = Yi_star(j);
        }
      } else {
        Yi_MH(iter, 0) = 0;
        for (int j = 1; j < Yi.size(); ++j) {
          Yi_MH(iter, j) = Yi(j);
        }
      }
    }

    mc.push_back(Yi_MH);
  }

  return mc;
}
