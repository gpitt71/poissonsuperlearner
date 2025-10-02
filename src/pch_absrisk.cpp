#include <Rcpp.h>
#include <limits>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pch_absolute_risk(IntegerVector id,
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

  NumericVector F(n);                 // cumulative incidence per row
  double logS = 0.0;                  // log survival at start of interval; log(1)=0
  double Fcum = 0.0;
  int prev_id = n ? id[0] : 0;

  const double LOG_DBL_MIN = std::log(std::numeric_limits<double>::min()); // ~ -708

  for (R_xlen_t i = 0; i < n; ++i) {
    if (id[i] != prev_id) {
      logS = 0.0;
      Fcum = 0.0;
      prev_id = id[i];
    }

    const double dti = dt[i];
    if (!R_finite(dti) || dti < 0.0)
      stop("dt must be finite and >= 0 at row %lld.", static_cast<long long>(i + 1));

    // Sum hazards and extract target-cause hazard
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
        if (na_is_zero) {
          // treat NA/Inf as 0
        } else {
          bad = true;
        }
      }
    }
    if (bad) { F[i] = NA_REAL; logS = NA_REAL; Fcum = NA_REAL; continue; }
    if (!R_finite(logS)) { F[i] = NA_REAL; continue; }

    // Compute increment in a stable way:
    double incr = 0.0;
    if (dti > 0.0 && hc > 0.0 && hsum > 0.0) {
      const double hsum_dt = hsum * dti;

      // stable 1 - exp(-x) as -expm1(-x)
      const double one_minus_e = -std::expm1(-hsum_dt);

      // log of the multiplicative step (excluding S): log(hc/hsum) + log(1 - e^{-hsum*dt})
      double log_step;
      if (hsum > 0.0) {
        // log(hc) - log(hsum) + log(one_minus_e)
        log_step = std::log(hc) - std::log(hsum) + std::log(one_minus_e);
      } else {
        log_step = -INFINITY; // unreachable since hsum>0 in this branch
      }

      const double log_incr = logS + log_step;

      // avoid underflow when exponentiating tiny values
      incr = (log_incr > LOG_DBL_MIN) ? std::exp(log_incr) : 0.0;

      // update survival in log-space
      logS -= hsum_dt;  // S <- S * exp(-hsum*dt)
    } else {
      // zero time or zero hazard or zero hc => no increment, survival unchanged if hsum==0
      if (hsum > 0.0 && dti > 0.0) logS -= hsum * dti;
    }

    Fcum += incr;
    F[i] = Fcum;
  }

  return F;
}
