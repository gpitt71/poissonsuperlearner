// Sparse Poisson IRLS with lasso coordinate descent for Phase 1 HAL.
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;
using Clock = std::chrono::steady_clock;

namespace {

inline double elapsed_ms(const Clock::time_point& start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

inline double soft_threshold(const double x, const double threshold) {
  if (x > threshold) return x - threshold;
  if (x < -threshold) return x + threshold;
  return 0.0;
}

struct SparseCSC {
  int n;
  int p;
  IntegerVector col_ptr;
  IntegerVector row_ind;
  NumericVector values;

  explicit SparseCSC(const S4& X) :
    col_ptr(X.slot("p")), row_ind(X.slot("i")), values(X.slot("x")) {
    IntegerVector dim = X.slot("Dim");
    if (dim.size() != 2) stop("`X` must have two dimensions.");
    n = dim[0];
    p = dim[1];
    if (col_ptr.size() != p + 1 || row_ind.size() != values.size())
      stop("Invalid dgCMatrix slots.");
  }

  void add_matvec(const std::vector<double>& beta,
                  std::vector<double>& out) const {
    for (int j = 0; j < p; ++j) {
      if (beta[j] == 0.0) continue;
      for (int k = col_ptr[j]; k < col_ptr[j + 1]; ++k)
        out[row_ind[k]] += values[k] * beta[j];
    }
  }
};

struct Timings {
  double irls = 0.0;
  double curvature = 0.0;
  double gradients = 0.0;
  double residual = 0.0;
  double checks = 0.0;
};

double objective(const SparseCSC& X,
                 const NumericVector& y,
                 const NumericVector& offset,
                 const std::vector<double>& beta,
                 const double intercept,
                 const double lambda,
                 const std::vector<double>& penalty,
                 std::vector<double>* eta_out = NULL,
                 std::vector<double>* mu_out = NULL) {
  std::vector<double> eta(X.n);
  for (int i = 0; i < X.n; ++i) eta[i] = offset[i] + intercept;
  X.add_matvec(beta, eta);

  double loss = 0.0;
  for (int i = 0; i < X.n; ++i) {
    if (!std::isfinite(eta[i]) || eta[i] > 700.0)
      return std::numeric_limits<double>::infinity();
    const double mu = std::exp(std::max(eta[i], -745.0));
    if (!std::isfinite(mu)) return std::numeric_limits<double>::infinity();
    loss += mu - y[i] * eta[i];
    if (mu_out) (*mu_out)[i] = mu;
  }
  double l1 = 0.0;
  for (int j = 0; j < X.p; ++j) l1 += penalty[j] * std::abs(beta[j]);
  if (eta_out) *eta_out = eta;
  return loss / static_cast<double>(X.n) + lambda * l1;
}

double actual_kkt(const SparseCSC& X,
                  const NumericVector& y,
                  const std::vector<double>& mu,
                  const std::vector<double>& beta,
                  const double lambda,
                  const std::vector<double>& penalty,
                  const bool intercept,
                  Timings& timing) {
  const auto start = Clock::now();
  double maximum = 0.0;
  for (int j = 0; j < X.p; ++j) {
    double gradient = 0.0;
    for (int k = X.col_ptr[j]; k < X.col_ptr[j + 1]; ++k) {
      const int i = X.row_ind[k];
      gradient += X.values[k] * (mu[i] - y[i]);
    }
    gradient /= static_cast<double>(X.n);
    const double bound = lambda * penalty[j];
    double violation;
    if (beta[j] > 1e-12) violation = std::abs(gradient + bound);
    else if (beta[j] < -1e-12) violation = std::abs(gradient - bound);
    else violation = std::max(0.0, std::abs(gradient) - bound);
    maximum = std::max(maximum, violation);
  }
  if (intercept) {
    double gradient = 0.0;
    for (int i = 0; i < X.n; ++i) gradient += mu[i] - y[i];
    maximum = std::max(maximum, std::abs(gradient / X.n));
  }
  timing.checks += elapsed_ms(start);
  return maximum;
}

} // namespace

// [[Rcpp::export]]
List hal_poisson_irls_cpp(S4 X,
                          NumericVector y,
                          NumericVector offset,
                          NumericVector lambda,
                          NumericVector penalty_factor,
                          bool intercept = true,
                          double tol = 1e-8,
                          int maxit_irls = 100,
                          int maxit_cd = 10000,
                          bool active_set = true,
                          bool rescale_penalty_factor = true,
                          int check_interrupt_every = 25) {
  if (!X.inherits("dgCMatrix")) stop("`X` must inherit from Matrix::dgCMatrix.");
  SparseCSC A(X);
  if (A.n < 1 || A.p < 1) stop("`X` must have at least one row and one column.");
  if (y.size() != A.n || offset.size() != A.n)
    stop("`y` and `offset` must have length nrow(X).");
  if (penalty_factor.size() != A.p)
    stop("`penalty_factor` must have length ncol(X).");
  if (lambda.size() < 1) stop("`lambda` must contain at least one value.");
  if (!(tol > 0.0) || maxit_irls < 1 || maxit_cd < 1)
    stop("Invalid convergence controls.");

  double penalty_sum = 0.0;
  double y_sum = 0.0;
  for (int i = 0; i < A.n; ++i) {
    if (!std::isfinite(y[i]) || y[i] < 0.0)
      stop("`y` must be finite and nonnegative.");
    if (!std::isfinite(offset[i])) stop("`offset` must be finite.");
    y_sum += y[i];
  }
  if (intercept && y_sum == 0.0)
    stop("An unpenalized intercept has no finite solution when all responses are zero.");
  for (int j = 0; j < A.p; ++j) {
    if (!std::isfinite(penalty_factor[j]) || penalty_factor[j] < 0.0)
      stop("`penalty_factor` must be finite and nonnegative.");
    penalty_sum += penalty_factor[j];
  }
  if (rescale_penalty_factor && penalty_sum <= 0.0)
    stop("At least one penalty factor must be positive when rescaling is enabled.");
  for (int k = 0; k < lambda.size(); ++k) {
    if (!std::isfinite(lambda[k]) || lambda[k] < 0.0)
      stop("`lambda` must be finite and nonnegative.");
    if (k && lambda[k] > lambda[k - 1]) stop("`lambda` must be decreasing.");
  }

  std::vector<double> pf(A.p);
  const double pf_scale = rescale_penalty_factor ? A.p / penalty_sum : 1.0;
  for (int j = 0; j < A.p; ++j) pf[j] = penalty_factor[j] * pf_scale;

  const int nlambda = lambda.size();
  NumericMatrix beta_path(A.p, nlambda);
  NumericVector intercept_path(nlambda), objective_path(nlambda), loss_path(nlambda),
    l1_path(nlambda), kkt_path(nlambda);
  IntegerVector irls_iterations(nlambda), cd_sweeps(nlambda), active_size(nlambda);
  LogicalVector converged(nlambda), cd_converged(nlambda);
  NumericMatrix timing(nlambda, 6);
  NumericMatrix detailed_ms(nlambda, 4);
  colnames(detailed_ms) = CharacterVector::create(
    "initial_full_sweep", "active_sweeps", "inactive_kkt_scans", "full_poisson_kkt");
  NumericMatrix counters(nlambda, 8);
  colnames(counters) = CharacterVector::create(
    "full_sweeps", "active_sweeps", "active_coordinate_visits",
    "full_coordinate_visits", "inactive_kkt_scans", "inactive_columns_checked",
    "coordinates_activated", "final_poisson_kkt_scans");
  colnames(timing) = CharacterVector::create(
    "irls_update", "curvature", "coordinate_gradients", "residual_updates",
    "active_and_convergence_checks", "total");

  std::vector<double> beta(A.p, 0.0);
  double b0 = 0.0;
  if (intercept) {
    double exposure_sum = 0.0;
    for (int i = 0; i < A.n; ++i) {
      if (offset[i] > 700.0) stop("Offset is too large for finite Poisson means.");
      exposure_sum += std::exp(std::max(offset[i], -745.0));
    }
    if (!(exposure_sum > 0.0) || !std::isfinite(exposure_sum))
      stop("Offsets imply a nonfinite or zero total exposure.");
    b0 = std::log(y_sum / exposure_sum);
  }

  const double n_inv = 1.0 / static_cast<double>(A.n);
  for (int path_index = 0; path_index < nlambda; ++path_index) {
    const auto total_start = Clock::now();
    Timings tm;
    double t_initial = 0, t_active = 0, t_inactive = 0, t_poisson_kkt = 0;
    double n_full = 0, n_active = 0, active_visits = 0, full_visits = 0;
    double n_inactive = 0, inactive_visits = 0, n_activated = 0, n_kkt = 0;
    const double lam = lambda[path_index];
    bool outer_ok = false;
    bool all_cd_ok = true;
    int total_sweeps = 0;
    double final_kkt = std::numeric_limits<double>::infinity();
    double previous_objective = objective(A, y, offset, beta, b0, lam, pf);

    for (int outer = 0; outer < maxit_irls; ++outer) {
      if (check_interrupt_every > 0 && outer % check_interrupt_every == 0)
        checkUserInterrupt();

      const std::vector<double> beta_old = beta;
      const double b0_old = b0;
      std::vector<double> eta(A.n), mu(A.n);
      auto step_start = Clock::now();
      const double current_objective = objective(
        A, y, offset, beta, b0, lam, pf, &eta, &mu);
      if (!std::isfinite(current_objective)) stop("Nonfinite Poisson objective.");
      std::vector<double> residual(A.n);
      for (int i = 0; i < A.n; ++i) {
        if (mu[i] == 0.0 && y[i] > 0.0)
          stop("IRLS encountered an underflowed fitted mean at a positive response.");
        residual[i] = mu[i] > 0.0 ? (y[i] - mu[i]) / mu[i] : -1.0;
        if (!std::isfinite(residual[i])) stop("Nonfinite IRLS working residual.");
      }
      tm.irls += elapsed_ms(step_start);

      step_start = Clock::now();
      std::vector<double> curvature(A.p, 0.0);
      for (int j = 0; j < A.p; ++j) {
        double value = 0.0;
        for (int k = A.col_ptr[j]; k < A.col_ptr[j + 1]; ++k) {
          const int i = A.row_ind[k];
          const double x = A.values[k];
          value += mu[i] * x * x;
        }
        curvature[j] = value * n_inv;
      }
      tm.curvature += elapsed_ms(step_start);

      std::vector<unsigned char> active(A.p, 0);
      for (int j = 0; j < A.p; ++j)
        active[j] = (beta[j] != 0.0 || pf[j] == 0.0);

      auto update_intercept = [&]() -> double {
        if (!intercept) return 0.0;
        double numerator = 0.0, denominator = 0.0;
        for (int i = 0; i < A.n; ++i) {
          numerator += mu[i] * residual[i];
          denominator += mu[i];
        }
        if (!(denominator > 0.0)) return 0.0;
        const double delta = numerator / denominator;
        if (delta != 0.0) {
          const auto start = Clock::now();
          for (int i = 0; i < A.n; ++i) residual[i] -= delta;
          tm.residual += elapsed_ms(start);
          b0 += delta;
        }
        return std::abs(delta);
      };

      auto update_coordinate = [&](const int j) -> double {
        if (!(curvature[j] > std::numeric_limits<double>::epsilon())) return 0.0;
        const auto gradient_start = Clock::now();
        double gradient = 0.0;
        for (int k = A.col_ptr[j]; k < A.col_ptr[j + 1]; ++k) {
          const int i = A.row_ind[k];
          gradient += mu[i] * A.values[k] * residual[i];
        }
        gradient *= n_inv;
        tm.gradients += elapsed_ms(gradient_start);
        const double candidate = soft_threshold(
          curvature[j] * beta[j] + gradient, lam * pf[j]) / curvature[j];
        const double delta = candidate - beta[j];
        if (delta != 0.0) {
          const auto residual_start = Clock::now();
          for (int k = A.col_ptr[j]; k < A.col_ptr[j + 1]; ++k)
            residual[A.row_ind[k]] -= A.values[k] * delta;
          tm.residual += elapsed_ms(residual_start);
          beta[j] = candidate;
          if (candidate != 0.0) active[j] = 1;
        }
        return std::abs(delta);
      };

      // One complete pass, followed by restricted active passes and inactive KKT checks.
      const auto initial_start = Clock::now();
      double max_change = update_intercept();
      for (int j = 0; j < A.p; ++j)
        max_change = std::max(max_change, update_coordinate(j));
      t_initial += elapsed_ms(initial_start);
      ++n_full;
      full_visits += A.p;
      ++total_sweeps;

      bool inner_ok = false;
      for (int sweep = 1; sweep < maxit_cd; ++sweep) {
        const auto active_start = Clock::now();
        max_change = update_intercept();
        for (int j = 0; j < A.p; ++j)
          if (!active_set || active[j]) {
            ++active_visits;
            max_change = std::max(max_change, update_coordinate(j));
          }
        t_active += elapsed_ms(active_start);
        ++n_active;
        ++total_sweeps;

        if (max_change <= tol * (1.0 + std::abs(b0))) {
          const auto check_start = Clock::now();
          bool added = false;
          if (active_set) {
            ++n_inactive;
            for (int j = 0; j < A.p; ++j) if (!active[j]) {
              ++inactive_visits;
              double gradient = 0.0;
              for (int k = A.col_ptr[j]; k < A.col_ptr[j + 1]; ++k) {
                const int i = A.row_ind[k];
                gradient += mu[i] * A.values[k] * residual[i];
              }
              gradient *= n_inv;
              if (std::abs(gradient) > lam * pf[j] + tol) {
                active[j] = 1;
                ++n_activated;
                added = true;
              }
            }
          }
          t_inactive += elapsed_ms(check_start);
          tm.checks += elapsed_ms(check_start);
          if (!added) { inner_ok = true; break; }
        }
        if (check_interrupt_every > 0 && sweep % check_interrupt_every == 0)
          checkUserInterrupt();
      }
      all_cd_ok = all_cd_ok && inner_ok;

      // Newton/IRLS can overshoot. Backtrack on the actual penalized likelihood.
      double candidate_objective = objective(A, y, offset, beta, b0, lam, pf);
      if (!std::isfinite(candidate_objective) ||
          candidate_objective > current_objective + 1e-12 * (1.0 + std::abs(current_objective))) {
        const std::vector<double> beta_candidate = beta;
        const double b0_candidate = b0;
        bool accepted = false;
        for (int ls = 1; ls <= 30; ++ls) {
          const double fraction = std::ldexp(1.0, -ls);
          for (int j = 0; j < A.p; ++j)
            beta[j] = beta_old[j] + fraction * (beta_candidate[j] - beta_old[j]);
          b0 = b0_old + fraction * (b0_candidate - b0_old);
          candidate_objective = objective(A, y, offset, beta, b0, lam, pf);
          if (std::isfinite(candidate_objective) && candidate_objective <= current_objective) {
            accepted = true;
            break;
          }
        }
        if (!accepted) { beta = beta_old; b0 = b0_old; candidate_objective = current_objective; }
      }

      std::vector<double> eta_new(A.n), mu_new(A.n);
      candidate_objective = objective(
        A, y, offset, beta, b0, lam, pf, &eta_new, &mu_new);
      const auto poisson_kkt_start = Clock::now();
      final_kkt = actual_kkt(A, y, mu_new, beta, lam, pf, intercept, tm);
      t_poisson_kkt += elapsed_ms(poisson_kkt_start);
      ++n_kkt;
      double parameter_change = std::abs(b0 - b0_old);
      for (int j = 0; j < A.p; ++j)
        parameter_change = std::max(parameter_change, std::abs(beta[j] - beta_old[j]));
      const double objective_change = std::abs(candidate_objective - current_objective) /
        (1.0 + std::abs(current_objective));
      previous_objective = candidate_objective;
      irls_iterations[path_index] = outer + 1;
      if (inner_ok && objective_change <= tol &&
          parameter_change <= std::sqrt(tol) * (1.0 + std::abs(b0)) &&
          final_kkt <= std::max(10.0 * tol, 1e-8)) {
        outer_ok = true;
        break;
      }
    }

    std::vector<double> eta(A.n), mu(A.n);
    const double final_objective = objective(A, y, offset, beta, b0, lam, pf, &eta, &mu);
    double penalty = 0.0;
    int nonzero = 0;
    for (int j = 0; j < A.p; ++j) {
      beta_path(j, path_index) = beta[j];
      penalty += pf[j] * std::abs(beta[j]);
      if (beta[j] != 0.0) ++nonzero;
    }
    double loss = 0.0;
    for (int i = 0; i < A.n; ++i) loss += mu[i] - y[i] * eta[i];
    intercept_path[path_index] = b0;
    objective_path[path_index] = final_objective;
    loss_path[path_index] = loss * n_inv;
    l1_path[path_index] = lam * penalty;
    kkt_path[path_index] = final_kkt;
    cd_sweeps[path_index] = total_sweeps;
    active_size[path_index] = nonzero;
    converged[path_index] = outer_ok;
    cd_converged[path_index] = all_cd_ok;
    timing(path_index, 0) = tm.irls;
    timing(path_index, 1) = tm.curvature;
    timing(path_index, 2) = tm.gradients;
    timing(path_index, 3) = tm.residual;
    timing(path_index, 4) = tm.checks;
    timing(path_index, 5) = elapsed_ms(total_start);
    detailed_ms(path_index, 0) = t_initial;
    detailed_ms(path_index, 1) = t_active;
    detailed_ms(path_index, 2) = t_inactive;
    detailed_ms(path_index, 3) = t_poisson_kkt;
    counters(path_index, 0) = n_full;
    counters(path_index, 1) = n_active;
    counters(path_index, 2) = active_visits;
    counters(path_index, 3) = full_visits;
    counters(path_index, 4) = n_inactive;
    counters(path_index, 5) = inactive_visits;
    counters(path_index, 6) = n_activated;
    counters(path_index, 7) = n_kkt;
  }

  NumericVector effective_pf(A.p);
  for (int j = 0; j < A.p; ++j) effective_pf[j] = pf[j];
  const double allocated_bytes =
    (9.0 * A.n + 6.0 * A.p) * sizeof(double) + A.p * sizeof(unsigned char);
  return List::create(
    _["beta"] = beta_path,
    _["a0"] = intercept_path,
    _["lambda"] = lambda,
    _["loss"] = loss_path,
    _["penalty"] = l1_path,
    _["objective"] = objective_path,
    _["max_kkt"] = kkt_path,
    _["irls_iterations"] = irls_iterations,
    _["cd_sweeps"] = cd_sweeps,
    _["active_size"] = active_size,
    _["converged"] = converged,
    _["cd_converged"] = cd_converged,
    _["timing_ms"] = timing,
    _["detailed_ms"] = detailed_ms,
    _["sweep_counters"] = counters,
    _["effective_penalty_factor"] = effective_pf,
    _["penalty_factor_scale"] = pf_scale,
    _["estimated_working_memory_bytes"] = allocated_bytes);
}
