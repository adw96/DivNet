// See http://minimallysufficient.github.io/r/programming/c++/2018/02/16/profiling-rcpp-packages.html

#include <Rcpp.h>
#include <gperftools/profiler.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP start_profiler(SEXP str) {
  ProfilerStart(as<const char*>(str));
  return R_NilValue;
}

// [[Rcpp::export]]
SEXP stop_profiler() {
  ProfilerStop();
  return R_NilValue;
}
