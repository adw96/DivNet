#include <RcppEigen.h>
#include <random>


using namespace std;
using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::export]]
Eigen::MatrixXd eigen_MCrow(const Eigen::VectorXd Yi, // length is Wi - 1
                            const Eigen::VectorXd Wi, // length is Yi + 1
                            const Eigen::VectorXd eYi,
                            const int Q, // no. OTUs aka length of Wi
                            const int base, // passed in as 1-based
                            const Eigen::MatrixXd sigInv,
                            const int MCiter,
                            const double stepsize)
{

  // Manually set the seed.
  default_random_engine generator;

  normal_distribution<double> norm_dist(0.0, stepsize);
  // This is the interval [0.0, 1.0)
  uniform_real_distribution<double> unif_dist(0.0, 1.0);

  // Wi minus the base OTU
  Eigen::VectorXd Wi_no_base(Wi.size() - 1);
  int cur_idx = 0;
  for (int i = 0; i < Wi.size(); ++i) {
    if (i != (base - 1)) {
      Wi_no_base(cur_idx++) = Wi(i);
    }
  }

  // Extra column for acceptance indicator.
  // MCiter rows and Q columns.
  Eigen::MatrixXd Yi_MH = Eigen::MatrixXd(MCiter, Q).setZero();

  Eigen::VectorXd Yi_star(Yi.size());
  double Eq5pt1 = 0.0,
    Eq5pt2 = 0.0,
    Eq5pt3 = 0.0,
    Eq5pt4 = 0.0,
    fullRat = 0.0,
    fullRat_exp = 0.0,
    acceptance = 0.0;

  // The main MCiter loop
  for (int iter = 0; iter < MCiter; ++iter) {
    for (int i = 0; i < Yi.size(); ++i) {
      Yi_star(i) = Yi(i) + norm_dist(generator);
    }

    // Denominator
    Eq5pt1 = Wi.sum() * (log(Yi     .array().exp().sum() + 1) -
                         log(Yi_star.array().exp().sum() + 1));

    // Numerator
    Eq5pt2 = (Wi_no_base * (Yi_star - Yi)).sum();

    Eq5pt3 = -0.5 * ((Yi_star - eYi).transpose() * sigInv) * (Yi_star - eYi);

    Eq5pt4 = -0.5 * ((Yi      - eYi).transpose() * sigInv) * (Yi      - eYi);

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

  return Yi_MH;
}
