#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix poisson_deviance_by_folder(const NumericMatrix& log_hazard,
                                         const NumericVector& tij,
                                         const IntegerVector& delta,
                                         const IntegerVector& folder, // 1..F
                                         const int F,
                                         const double eps = 1e-15) {
  const int n = log_hazard.nrow();
  const int L = log_hazard.ncol();

  if (tij.size() != n) stop("tij length must equal nrow(log_hazard).");
  if (delta.size() != n) stop("delta length must equal nrow(log_hazard).");
  if (folder.size() != n) stop("folder length must equal nrow(log_hazard).");
  if (F <= 0) stop("F must be positive.");

  NumericMatrix dev(F, L);
  std::fill(dev.begin(), dev.end(), 0.0);

  for (int i = 0; i < n; ++i) {
    const int f = folder[i];
    if (f < 1 || f > F) stop("folder must be in 1..F.");
    const int y = delta[i];
    if (y != 0 && y != 1) stop("delta must be 0/1.");
    const double ti = tij[i];
    if (!R_finite(ti) || ti < 0.0) stop("tij must be finite and nonnegative.");

    const int fi = f - 1; // 0-based row
    for (int l = 0; l < L; ++l) {
      const double eta = log_hazard(i, l);
      if (!R_finite(eta)) continue; // or stop(...)
      double mu = std::exp(eta) * ti;
      if (!R_finite(mu) || mu < eps) mu = eps;

      // contribution to deviance/2 for y in {0,1}
      double c = (y == 0) ? mu : (-std::log(mu) - 1.0 + mu);
      dev(fi, l) += 2.0 * c;
    }
  }

  return dev; // F x L (already deviance, not dev/2)
}
