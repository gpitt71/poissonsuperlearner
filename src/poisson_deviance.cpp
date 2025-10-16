// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <unordered_map>
#include <string>
using namespace Rcpp;

// safe deviance term for a single Poisson observation y ~ Poisson(mu)
// returns: y*log(y/mu) - (y - mu), with the convention 0*log(0)=0
inline double dev_term(int y, double mu) {
  if (y == 0) return mu;                 // 0*log(0/mu) - (0 - mu) = mu
  if (mu <= 0.0) return R_PosInf;        // positive count but zero/neg mean -> infinite deviance
  return y * std::log(static_cast<double>(y) / mu) - (y - mu);
}

// [[Rcpp::export]]
DataFrame poisson_deviance_cpp(IntegerVector folder,
                               CharacterVector learner,
                               NumericMatrix lambda,   // columns = hazards
                               NumericVector tij,
                               IntegerMatrix delta)    // columns = hazards
{
  const int n = lambda.nrow();
  const int H = lambda.ncol();

  if (delta.nrow() != n || delta.ncol() != H)
    stop("lambda and delta must have the same nrow and ncol.");
  if (tij.size() != n || folder.size() != n || learner.size() != n)
    stop("folder, learner, tij must have length equal to nrow(lambda).");

  // Aggregate sum of contributions per (learner, folder)
  // key format: <learner>\t<folder>
  std::unordered_map<std::string, double> sum_by_lrn_fold;
  sum_by_lrn_fold.reserve(n);

  for (int i = 0; i < n; ++i) {
    const double t = tij[i];
    double row_sum = 0.0;

    for (int h = 0; h < H; ++h) {
      const double mu = lambda(i, h) * t;
      const int y = delta(i, h);
      const double term = dev_term(y, mu);
      if (R_finite(term)) {
        row_sum += term;
      } else {
        row_sum = R_PosInf; // any infinite term makes the row/group infinite
        break;
      }
    }

    // build key
    std::string key = std::string(learner[i]) + '\t' + std::to_string(folder[i]);
    auto it = sum_by_lrn_fold.find(key);
    if (it == sum_by_lrn_fold.end()) {
      sum_by_lrn_fold.emplace(std::move(key), row_sum);
    } else {
      // associative: sum over rows within learner×folder
      it->second += row_sum;
    }
  }

  // Now compute deviance_v = 2 * sum(contrib) per learner×folder,
  // then deviance = mean(deviance_v) across folders per learner.
  std::unordered_map<std::string, double> sum_dev_per_learner;
  std::unordered_map<std::string, int>    n_folders_per_learner;

  for (const auto &kv : sum_by_lrn_fold) {
    // parse key -> learner and folder
    const std::string &key = kv.first;
    const double sum_contrib = kv.second;
    // extract learner up to '\t'
    const std::size_t pos = key.find('\t');
    const std::string lrn = key.substr(0, pos);

    const double dev_v = R_finite(sum_contrib) ? 2.0 * sum_contrib : R_PosInf;

    auto itS = sum_dev_per_learner.find(lrn);
    if (itS == sum_dev_per_learner.end()) {
      sum_dev_per_learner.emplace(lrn, dev_v);
      n_folders_per_learner.emplace(lrn, 1);
    } else {
      // average over folders: keep running sum and count
      itS->second += dev_v;
      n_folders_per_learner[lrn] += 1;
    }
  }

  // materialize outputs
  int L = static_cast<int>(sum_dev_per_learner.size());
  CharacterVector out_learner(L);
  NumericVector   out_deviance(L);

  int j = 0;
  for (const auto &kv : sum_dev_per_learner) {
    const std::string &lrn = kv.first;
    const double sum_dev = kv.second;
    const int cnt = n_folders_per_learner[lrn];

    out_learner[j] = lrn;
    out_deviance[j] = sum_dev / static_cast<double>(cnt); // mean over folders
    ++j;
  }

  return DataFrame::create(_["learner"] = out_learner,
                           _["deviance"] = out_deviance,
                           _["stringsAsFactors"] = false);
}
