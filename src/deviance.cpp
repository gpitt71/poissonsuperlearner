#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix poisson_deviance_by_folder_cols(List log_hazard_cols,
                                              NumericVector tij,
                                              IntegerVector delta,
                                              IntegerVector folder,
                                              int nfold,
                                              double eps = 1e-15) {

  const int n = tij.size();
  const int L = log_hazard_cols.size();

  if (delta.size() != n || folder.size() != n)
    stop("Length mismatch among tij/delta/folder.");

  if (nfold <= 0)
    stop("nfold must be positive.");

  if (L == 0)
    stop("No learner columns supplied.");

  // rows = folds, cols = learners
  NumericMatrix out(nfold, L);

  // Convert folds to 0-based indexing
  IntegerVector fold0(n);
  for (int i = 0; i < n; ++i) {
    int f = folder[i];
    if (f < 1 || f > nfold)
      stop("folder values must be in 1..nfold.");
    fold0[i] = f - 1;
  }

  // Loop learners
  for (int j = 0; j < L; ++j) {

    NumericVector lh = log_hazard_cols[j];
    if (lh.size() != n)
      stop("Learner column length mismatch.");

    for (int i = 0; i < n; ++i) {

      const double mu = std::exp(lh[i]) * tij[i];
      const int y = delta[i];
      const double mu_safe = (mu > eps) ? mu : eps;

      // y*log(y/mu) - (y - mu), with y ∈ {0,1}
      double contrib = 0.0;

      if (y == 1)
        contrib += std::log(1.0 / mu_safe);

      contrib -= ( (double)y - mu );

      out(fold0[i], j) += contrib;
    }
  }

  return out;
}
