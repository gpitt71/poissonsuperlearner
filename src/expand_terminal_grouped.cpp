#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

static inline bool is_bad_double(double x) {
  return Rcpp::NumericVector::is_na(x) || !std::isfinite(x);
}

static inline bool is_kept_cell(double risk, double events, double nterm) {
  return (risk != 0.0) || (events != 0.0) || (nterm != 0.0);
}

static void validate_widths(const NumericVector& widths) {
  const int K = widths.size();
  if (K < 1) stop("'widths' must contain at least one interval width.");

  for (int k = 0; k < K; ++k) {
    const double w = widths[k];
    if (is_bad_double(w) || w <= 0.0) {
      stop("All 'widths' must be finite and strictly positive.");
    }
  }
}

static void validate_inputs(const IntegerVector& gid,
                            const IntegerVector& node_ix,
                            const NumericVector& tij,
                            const NumericVector& deltaij,
                            const NumericVector& N,
                            const NumericVector& widths) {
  const int n = gid.size();
  const int K = widths.size();

  if (node_ix.size() != n || tij.size() != n ||
      deltaij.size() != n || N.size() != n) {
    stop("All input vectors must have the same length.");
  }

  validate_widths(widths);

  int last_gid = 0;

  for (int i = 0; i < n; ++i) {
    const int g = gid[i];
    const int j = node_ix[i];
    const double y = tij[i];
    const double d = deltaij[i];
    const double nn = N[i];

    if (g == NA_INTEGER || g < 1) {
      stop("'gid' must be positive and non-missing.");
    }

    if (g < last_gid) {
      stop("Input must be sorted by 'gid'. In R, run: data.table::setorder(x, gid, node_ix).");
    }

    last_gid = g;

    if (j == NA_INTEGER || j < 1 || j > K) {
      stop("'node_ix' must be an interval index in 1:length(widths).");
    }

    if (is_bad_double(y) || y < 0.0) {
      stop("'tij' must be finite and non-negative.");
    }

    if (is_bad_double(d) || d < 0.0) {
      stop("'deltaij' must be finite and non-negative.");
    }

    if (is_bad_double(nn) || nn < 0.0) {
      stop("'N' must be finite and non-negative.");
    }
  }
}

static std::size_t count_output_cells(const IntegerVector& gid,
                                      const IntegerVector& node_ix,
                                      const NumericVector& tij,
                                      const NumericVector& deltaij,
                                      const NumericVector& N,
                                      const NumericVector& widths) {
  const int n = gid.size();
  const int K = widths.size();

  std::vector<double> risk(K, 0.0);
  std::vector<double> events(K, 0.0);
  std::vector<double> nterm(K, 0.0);
  std::vector<double> diffN(K + 1, 0.0);

  std::size_t m = 0;
  int i = 0;

  while (i < n) {
    const int g = gid[i];

    std::fill(risk.begin(), risk.end(), 0.0);
    std::fill(events.begin(), events.end(), 0.0);
    std::fill(nterm.begin(), nterm.end(), 0.0);
    std::fill(diffN.begin(), diffN.end(), 0.0);

    while (i < n && gid[i] == g) {
      const int jj = node_ix[i] - 1;
      const double nn = N[i];

      risk[jj]   += tij[i];
      events[jj] += deltaij[i];
      nterm[jj]  += nn;

      if (jj > 0 && nn != 0.0) {
        diffN[0]  += nn;
        diffN[jj] -= nn;
      }

      ++i;
    }

    double runningN = 0.0;

    for (int jj = 0; jj < K; ++jj) {
      runningN += diffN[jj];

      if (runningN != 0.0) {
        risk[jj] += runningN * widths[jj];
      }

      if (is_kept_cell(risk[jj], events[jj], nterm[jj])) {
        ++m;
      }
    }
  }

  return m;
}

// [[Rcpp::export]]
List expand_terminal_grouped_cpp(IntegerVector gid,
                                 IntegerVector node_ix,
                                 NumericVector tij,
                                 NumericVector deltaij,
                                 NumericVector N,
                                 NumericVector widths) {
  validate_inputs(gid, node_ix, tij, deltaij, N, widths);

  const int n = gid.size();
  const int K = widths.size();

  const std::size_t m = count_output_cells(
    gid, node_ix, tij, deltaij, N, widths
  );

  IntegerVector out_gid(m);
  IntegerVector out_node_ix(m);
  NumericVector out_deltaij(m);
  NumericVector out_tij(m);
  NumericVector out_N_terminal(m);

  std::vector<double> risk(K, 0.0);
  std::vector<double> events(K, 0.0);
  std::vector<double> nterm(K, 0.0);
  std::vector<double> diffN(K + 1, 0.0);

  std::size_t pos = 0;
  int i = 0;

  while (i < n) {
    const int g = gid[i];

    std::fill(risk.begin(), risk.end(), 0.0);
    std::fill(events.begin(), events.end(), 0.0);
    std::fill(nterm.begin(), nterm.end(), 0.0);
    std::fill(diffN.begin(), diffN.end(), 0.0);

    while (i < n && gid[i] == g) {
      const int jj = node_ix[i] - 1;
      const double nn = N[i];

      risk[jj]   += tij[i];
      events[jj] += deltaij[i];
      nterm[jj]  += nn;

      if (jj > 0 && nn != 0.0) {
        diffN[0]  += nn;
        diffN[jj] -= nn;
      }

      ++i;
    }

    double runningN = 0.0;

    for (int jj = 0; jj < K; ++jj) {
      runningN += diffN[jj];

      if (runningN != 0.0) {
        risk[jj] += runningN * widths[jj];
      }

      if (is_kept_cell(risk[jj], events[jj], nterm[jj])) {
        out_gid[pos]        = g;
        out_node_ix[pos]    = jj + 1;
        out_deltaij[pos]    = events[jj];
        out_tij[pos]        = risk[jj];
        out_N_terminal[pos] = nterm[jj];
        ++pos;
      }
    }
  }

  return List::create(
    _["gid"]        = out_gid,
    _["node_ix"]    = out_node_ix,
    _["deltaij"]    = out_deltaij,
    _["tij"]        = out_tij,
    _["N_terminal"] = out_N_terminal
  );
}
