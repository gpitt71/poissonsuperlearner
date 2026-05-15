#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix poisson_deviance_by_folder_hazard_cols(
    const DataFrame& data,
    const CharacterVector& hazard_cols,
    const NumericVector& tij,
    const IntegerVector& delta,
    const IntegerVector& fold,   // 1..nfold
    const int nfold,
    const double eps = 1e-15
) {
  const int n = tij.size();
  const int L = hazard_cols.size();

  if (delta.size() != n) stop("delta length must equal tij length.");
  if (fold.size() != n) stop("fold length must equal tij length.");
  if (nfold <= 0) stop("nfold must be positive.");

  NumericMatrix dev(nfold, L);
  std::fill(dev.begin(), dev.end(), 0.0);

  for (int i = 0; i < n; ++i) {
    if (fold[i] < 1 || fold[i] > nfold) {
      stop("fold must be in 1..nfold.");
    }

    if (delta[i] != 0 && delta[i] != 1) {
      stop("delta must be 0/1.");
    }

    if (!R_finite(tij[i]) || tij[i] < 0.0) {
      stop("tij must be finite and nonnegative.");
    }
  }

  for (int l = 0; l < L; ++l) {
    const std::string colname = as<std::string>(hazard_cols[l]);

    if (!data.containsElementNamed(colname.c_str())) {
      stop("hazard column not found: " + colname);
    }

    NumericVector h_col = data[colname];

    if (h_col.size() != n) {
      stop("hazard column has wrong length: " + colname);
    }

    for (int i = 0; i < n; ++i) {
      const int fi = fold[i] - 1;

      double current = dev(fi, l);

      // If this learner-fold already failed, keep it as NA.
      if (NumericVector::is_na(current)) {
        continue;
      }

      const double h = h_col[i];

      // Invalid hazard prediction poisons this learner-fold.
      if (!R_finite(h) || h < 0.0) {
        dev(fi, l) = NA_REAL;
        continue;
      }

      double mu = h * tij[i];

      if (!R_finite(mu) || mu < eps) {
        mu = eps;
      }

      const int y = delta[i];

      const double c = (y == 0)
        ? mu
      : (-std::log(mu) - 1.0 + mu);

      dev(fi, l) = current + 2.0 * c;
    }
  }

  colnames(dev) = hazard_cols;

  return dev;
}
