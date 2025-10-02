#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pch_survival(IntegerVector id,
                           NumericVector dt,
                           NumericMatrix haz,
                           const bool na_is_zero = false) {
  const R_xlen_t n = id.size();
  if (dt.size() != n || haz.nrow() != n)
    stop("Length mismatch: id, dt, and haz rows must match.");
  if (haz.ncol() < 1)
    stop("haz must have at least one column (one cause).");

  NumericVector S(n);
  double cumhaz = 0.0;
  int prev_id = n > 0 ? id[0] : 0;

  for (R_xlen_t i = 0; i < n; ++i) {
    if (id[i] != prev_id) {
      cumhaz = 0.0;
      prev_id = id[i];
    }
    const double dti = dt[i];
    if (!R_finite(dti) || dti < 0.0)
      stop("dt must be finite and nonnegative at row %lld.",
           static_cast<long long>(i + 1));

    double hsum = 0.0;
    bool bad = false;
    for (int j = 0; j < haz.ncol(); ++j) {
      const double hij = haz(i, j);
      if (R_finite(hij)) {
        if (hij < 0.0)
          stop("hazards must be >= 0 at row %lld, col %d.",
               static_cast<long long>(i + 1), j + 1);
        hsum += hij;
      } else {
        if (!na_is_zero) { bad = true; }
        // else treat NA/Inf as 0
      }
    }
    if (bad) { S[i] = NA_REAL; cumhaz = NA_REAL; continue; }
    if (!R_finite(cumhaz)) { S[i] = NA_REAL; continue; }

    cumhaz += dti * hsum;
    S[i] = std::exp(-cumhaz);
  }
  return S;
}
