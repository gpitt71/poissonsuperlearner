# Poisson SuperLearner

The package provides an implementation of piece-wise constant hazard models for time-to-event analysis of survival and competing risks data. The piecewise constant hazard models can be combined in an ensemble, the Poisson Superlearner, via cross-validated risk minimization for flexible hazard estimation. It enables estimation of survival functions and risk predictions.

The package provides:

* Flexible specification of base learners (e.g. penalized and unpenalized Poisson regression).
* Automatic time-splitting into piecewise-constant hazard intervals.
* Cross-validated ensemble learning.
* Absolute risk prediction at arbitrary time horizons.

---

## Installation

```r
# install.packages("devtools")
# devtools::install_github("gpitt71/poissonsuperlearner")

library(poissonsuperlearner)
```

---

## Example: Piecewise-Constant Hazard Model

Fit a single PCH model and obtain absolute risk predictions.

```r
library(poissonsuperlearner)

set.seed(42)

# Simulate synthetic survival data
d <- simulateStenoT1(
  n = 50,
  scenario = "alpha"
)

# Define an unpenalized Poisson hazard learner
l_glm <- Learner_glmnet(
  covariates = c("sex", "diabetes_duration"),
  cross_validation = FALSE,
  lambda = 0,
  intercept = TRUE,
  penalise_nodes = FALSE
)

# Fit piecewise-constant hazard model
fit_glm <- fit_learner(
  data = d,
  id = "id",
  status = "event",
  event_time = "time",
  learner = l_glm,
  number_of_nodes = 5
)

# Absolute risk prediction at time horizon t = 5
predictRisk(fit_glm, newdata = d[1, ], times = 5)
```

**What happens internally?**

* The follow-up time is split into `number_of_nodes` intervals.
* A Poisson regression is fitted to the expanded dataset.
* The estimated hazards are integrated to produce survival and absolute risk.

---

## Example: Superlearner for Piecewise Hazards

Combine multiple hazard learners into an ensemble.

```r
library(poissonsuperlearner)

set.seed(42)

d <- simulateStenoT1(
  n = 50,
  scenario = "alpha"
)

# Base learner 1: unpenalized Poisson regression
l_glm <- Learner_glmnet(
  covariates = c("sex", "diabetes_duration"),
  cross_validation = FALSE,
  lambda = 0,
  intercept = TRUE,
  penalise_nodes = FALSE
)

# Base learner 2: Lasso-penalized Poisson regression
l_lasso <- Learner_glmnet(
  covariates = c("value_Smoking", "value_LDL"),
  cross_validation = TRUE,
  alpha = 1,
  intercept = TRUE,
  penalise_nodes = FALSE
)

learners_list <- list(
  glm   = l_glm,
  lasso = l_lasso
)

# Fit the superlearner
sl_fit <- Superlearner(
  data = d,
  id = "id",
  status = "event",
  event_time = "time",
  learners = learners_list,
  number_of_nodes = 5
)

# Absolute risk prediction from the ensemble
predictRisk(sl_fit, newdata = d[1, ], times = 5)
```

**Superlearner workflow**

1. Each base learner is fitted on cross-validation folds.
2. Out-of-sample Poisson deviance is computed.
3. A meta-learner combines the base hazard predictions.
4. Final hazards are aggregated and converted to absolute risk.

---

## Main Functions

* `Learner_glmnet()` – define a Poisson hazard learner.
* `fit_learner()` – fit a single piecewise-constant hazard model.
* `Superlearner()` – fit an ensemble of hazard learners.
* `predictRisk()` – obtain absolute risk predictions.

---

## Typical Use Case

`poissonsuperlearner` is designed for:

* Survival analysis via Poisson representations.
* Flexible hazard modeling with regularization.
* Ensemble-based absolute risk prediction.
* Methodological research in time-to-event modeling.

The package focuses on modular learners, transparent cross-validation, and direct control of the piecewise hazard structure.
