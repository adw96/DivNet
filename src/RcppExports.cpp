// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// mcrow_full_ratio
double mcrow_full_ratio(const Eigen::VectorXd& Wi, const Eigen::VectorXd& Wi_no_base, const Eigen::VectorXd& Yi, const Eigen::VectorXd& Yi_star, const Eigen::VectorXd& eYi, const Eigen::MatrixXd& sigma_inverse);
RcppExport SEXP _DivNet_mcrow_full_ratio(SEXP WiSEXP, SEXP Wi_no_baseSEXP, SEXP YiSEXP, SEXP Yi_starSEXP, SEXP eYiSEXP, SEXP sigma_inverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Wi(WiSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Wi_no_base(Wi_no_baseSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Yi(YiSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Yi_star(Yi_starSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type eYi(eYiSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type sigma_inverse(sigma_inverseSEXP);
    rcpp_result_gen = Rcpp::wrap(mcrow_full_ratio(Wi, Wi_no_base, Yi, Yi_star, eYi, sigma_inverse));
    return rcpp_result_gen;
END_RCPP
}
// mcrow_MCrow
Eigen::MatrixXd mcrow_MCrow(const Eigen::VectorXd& Yi, const Eigen::VectorXd& Wi, const Eigen::VectorXd& eYi, const int base, const Eigen::MatrixXd& sigInv, const int MCiter, const double stepsize);
RcppExport SEXP _DivNet_mcrow_MCrow(SEXP YiSEXP, SEXP WiSEXP, SEXP eYiSEXP, SEXP baseSEXP, SEXP sigInvSEXP, SEXP MCiterSEXP, SEXP stepsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Yi(YiSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Wi(WiSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type eYi(eYiSEXP);
    Rcpp::traits::input_parameter< const int >::type base(baseSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type sigInv(sigInvSEXP);
    Rcpp::traits::input_parameter< const int >::type MCiter(MCiterSEXP);
    Rcpp::traits::input_parameter< const double >::type stepsize(stepsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(mcrow_MCrow(Yi, Wi, eYi, base, sigInv, MCiter, stepsize));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DivNet_mcrow_full_ratio", (DL_FUNC) &_DivNet_mcrow_full_ratio, 6},
    {"_DivNet_mcrow_MCrow", (DL_FUNC) &_DivNet_mcrow_MCrow, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_DivNet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
