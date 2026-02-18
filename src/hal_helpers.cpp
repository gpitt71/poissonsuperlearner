#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector inter2_cpp(IntegerVector a, IntegerVector b) {
  const int la = a.size(), lb = b.size();
  if (la == 0 || lb == 0) return IntegerVector(0);

  IntegerVector out(std::min(la, lb));
  int i = 0, j = 0, k = 0;

  while (i < la && j < lb) {
    int ai = a[i], bj = b[j];
    if (ai < bj) {
      ++i;
    } else if (ai > bj) {
      ++j;
    } else {
      out[k++] = ai;
      ++i; ++j;
    }
  }

  out.erase(k, out.size());
  return out;
}

// [[Rcpp::export]]
IntegerVector interN_cpp(List lst) {
  const int m = lst.size();
  if (m == 0) return IntegerVector(0);
  if (m == 1) return as<IntegerVector>(lst[0]);

  IntegerVector acc = inter2_cpp(as<IntegerVector>(lst[0]),
                                 as<IntegerVector>(lst[1]));
  if (acc.size() == 0 || m == 2) return acc;

  for (int i = 2; i < m; ++i) {
    acc = inter2_cpp(acc, as<IntegerVector>(lst[i]));
    if (acc.size() == 0) break;
  }
  return acc;
}

// Helper: unique on sorted vector (keeps order)
static inline void unique_sorted_in_place(std::vector<double>& v) {
  v.erase(std::unique(v.begin(), v.end()), v.end());
}

// [[Rcpp::export]]
List mk_main_numeric_cpp(NumericVector x, int K) {
  const int n = x.size();

  // Collect finite, non-NA values for quantiles
  std::vector<double> vals;
  vals.reserve(n);
  for (int i = 0; i < n; ++i) {
    double xi = x[i];
    if (!NumericVector::is_na(xi) && std::isfinite(xi)) vals.push_back(xi);
  }

  if (vals.empty()) {
    return List::create(
      _["cutpoints"] = NumericVector(0),
      _["idxs"]      = List(0)
    );
  }

  std::sort(vals.begin(), vals.end());
  const int m = (int)vals.size();
  const double min_x = vals.front();

  // Type=1 quantiles at probs seq(0,1,length=K+2)
  std::vector<double> cps;
  cps.reserve(std::max(0, K));

  if (K > 0) {
    for (int j = 0; j < K + 2; ++j) {
      double p = (double)j / (double)(K + 1);
      // type=1: index = max(1, ceil(p*m))
      int idx = (int)std::ceil(p * m);
      if (idx < 1) idx = 1;
      if (idx > m) idx = m;
      cps.push_back(vals[idx - 1]);
    }
  } else {
    // K=0 -> endpoints only, then you drop them -> no cutpoints
    cps.push_back(vals.front());
    cps.push_back(vals.back());
  }

  // Remove non-finite (shouldn't happen due to vals filter, but keep consistent)
  cps.erase(std::remove_if(cps.begin(), cps.end(),
                           [](double z){ return !std::isfinite(z); }),
                           cps.end());

  // Drop endpoints
  if ((int)cps.size() >= 2) {
    cps.erase(cps.begin());
    cps.pop_back();
  } else {
    cps.clear();
  }

  if (cps.empty()) {
    return List::create(
      _["cutpoints"] = NumericVector(0),
      _["idxs"]      = List(0)
    );
  }

  std::sort(cps.begin(), cps.end());
  unique_sorted_in_place(cps);

  // Standard HAL primitive is I(x >= c).
  // If c == min(x) then the column is all-1; drop it (and anything below, for safety under ties).
  cps.erase(std::remove_if(cps.begin(), cps.end(),
                           [&](double c){ return c <= min_x; }),
                                                    cps.end());

  if (cps.empty()) {
    return List::create(
      _["cutpoints"] = NumericVector(0),
      _["idxs"]      = List(0)
    );
  }

  // Order indices by x (ascending), excluding NA / non-finite
  std::vector<int> ord;
  ord.reserve(m);
  for (int i = 0; i < n; ++i) {
    double xi = x[i];
    if (!NumericVector::is_na(xi) && std::isfinite(xi)) ord.push_back(i + 1); // 1-based
  }

  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
    // a,b are 1-based
    return x[a - 1] < x[b - 1];
  });

  // Precompute sorted xs aligned with ord
  std::vector<double> xs;
  xs.reserve(ord.size());
  for (int idx : ord) xs.push_back(x[idx - 1]);

  // For each cutpoint c: take indices with xs >= c (I(x >= c))
  List idxs((int)cps.size());
  const int N = (int)xs.size();

  for (int j = 0; j < (int)cps.size(); ++j) {
    double c = cps[j];
    auto it = std::lower_bound(xs.begin(), xs.end(), c); // first >= c
    int pos = (int)(it - xs.begin());
    int len = N - pos;

    if (len <= 0) {
      idxs[j] = IntegerVector(0);
      continue;
    }

    // Collect row indices, then SORT ASCENDING by row index for inter2_cpp compatibility
    IntegerVector out(len);
    for (int t = 0; t < len; ++t) out[t] = ord[pos + t];

    std::sort(out.begin(), out.end());  // critical for correct intersections

    idxs[j] = out;
  }

  return List::create(
    _["cutpoints"] = wrap(cps),
    _["idxs"]      = idxs
  );
}

// [[Rcpp::export]]
List mk_main_factor_cpp(IntegerVector x_int, int n_levels) {
  const int n = x_int.size();
  std::vector< std::vector<int> > bins(n_levels);
  for (int l = 0; l < n_levels; ++l) bins[l].reserve(n / std::max(1, n_levels));

  for (int i = 0; i < n; ++i) {
    int code = x_int[i];
    if (code == NA_INTEGER) continue;
    if (code < 1 || code > n_levels) continue;
    bins[code - 1].push_back(i + 1); // 1-based row index
  }

  List idxs(n_levels);
  for (int l = 0; l < n_levels; ++l) {
    IntegerVector out((int)bins[l].size());
    for (int t = 0; t < (int)bins[l].size(); ++t) out[t] = bins[l][t];
    idxs[l] = out;
  }

  return List::create(_["idxs"] = idxs);
}


// [[Rcpp::export]]
List add_cols_cpp(List idxs_list, int p_start) {
  const int k = idxs_list.size();

  LogicalVector keep(k);
  int ncol = 0;
  R_xlen_t nnz = 0;

  // pass 1: keep + sizes
  for (int t = 0; t < k; ++t) {
    SEXP el = idxs_list[t];
    if (el == R_NilValue) { keep[t] = false; continue; }

    IntegerVector idx(el);
    if (idx.size() == 0) { keep[t] = false; continue; }

    keep[t] = true;
    ++ncol;
    nnz += idx.size();
  }

  if (ncol == 0 || nnz == 0) {
    return List::create(
      _["i"]    = IntegerVector(0),
      _["j"]    = IntegerVector(0),
      _["keep"] = keep,
      _["ncol"] = 0
    );
  }

  IntegerVector i_out(nnz);
  IntegerVector j_out(nnz);

  // pass 2: fill
  R_xlen_t pos = 0;
  int col_offset = 0;

  for (int t = 0; t < k; ++t) {
    if (!keep[t]) continue;
    IntegerVector idx = idxs_list[t];

    const int col_id = p_start + col_offset + 1; // 1-based col index
    const R_xlen_t m = idx.size();

    for (R_xlen_t r = 0; r < m; ++r) {
      i_out[pos] = idx[r];
      j_out[pos] = col_id;
      ++pos;
    }
    ++col_offset;
  }

  return List::create(
    _["i"]    = i_out,
    _["j"]    = j_out,
    _["keep"] = keep,
    _["ncol"] = ncol
  );
}
