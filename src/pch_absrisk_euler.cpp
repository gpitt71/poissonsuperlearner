#include <Rcpp.h>
#include <limits>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pch_absolute_risk_euler(IntegerVector id,
                                      NumericVector dt,
                                      NumericMatrix haz,
                                      int cause_idx,
                                      const bool one_based = true,
                                      const bool na_is_zero = false) {
  const R_xlen_t n = id.size();
  if (dt.size() != n || haz.nrow() != n) stop("Length mismatch.");
  if (haz.ncol() < 1) stop("haz must have >= 1 column.");
  if (one_based) {
    if (cause_idx < 1 || cause_idx > haz.ncol()) stop("cause_idx out of range.");
    --cause_idx;
  } else {
    if (cause_idx < 0 || cause_idx >= haz.ncol()) stop("cause_idx out of range.");
  }

  NumericVector F(n);     // cumulative incidence per row
  double logS = 0.0;      // log survival at start of interval (log(1)=0)
  double Fcum = 0.0;
  int prev_id = n ? id[0] : 0;

  const double LOG_DBL_MIN = std::log(std::numeric_limits<double>::min()); // ~ -708

  for (R_xlen_t i = 0; i < n; ++i) {
    // restart for a new individual
    if (i == 0 || id[i] != prev_id) {
      logS  = 0.0;
      Fcum  = 0.0;
      prev_id = id[i];
    }

    const double dti = dt[i];
    if (!R_finite(dti) || dti < 0.0)
      stop("dt must be finite and >= 0 at row %lld.", static_cast<long long>(i + 1));

    // sum of hazards in the interval and target-cause hazard
    double hsum = 0.0, hc = 0.0;
    bool bad = false;
    for (int j = 0; j < haz.ncol(); ++j) {
      const double hij = haz(i, j);
      if (R_finite(hij)) {
        if (hij < 0.0)
          stop("hazards must be >= 0 at row %lld, col %d.",
               static_cast<long long>(i + 1), j + 1);
        hsum += hij;
        if (j == cause_idx) hc = hij;
      } else {
        if (!na_is_zero) bad = true; // NA/Inf not allowed unless treated as 0
      }
    }
    if (bad) { F[i] = NA_REAL; logS = NA_REAL; Fcum = NA_REAL; continue; }
    if (!R_finite(logS)) { F[i] = NA_REAL; continue; }

    // Euler increment: S(t_{k-1}) * hc * dt
    double incr = 0.0;
    if (dti > 0.0 && hc > 0.0) {
      const double log_incr = logS + std::log(hc) + std::log(dti);
      incr = (log_incr > LOG_DBL_MIN) ? std::exp(log_incr) : 0.0;
    }

    // update survival using piece-wise constant total hazard over the interval
    if (hsum > 0.0 && dti > 0.0) logS -= hsum * dti;

    Fcum += incr;
    F[i] = Fcum;
  }

  return F;
}
