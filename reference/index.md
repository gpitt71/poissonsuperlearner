# Package index

## Fit models

Define base learners and fit piecewise-constant hazard models.

- [`Learner_glmnet-class`](https://gpitt71.github.io/poissonsuperlearner/reference/Learner_glmnet-class.md)
  [`Learner_glmnet`](https://gpitt71.github.io/poissonsuperlearner/reference/Learner_glmnet-class.md)
  :

  Penalized Poisson learner via `glmnet`

- [`Learner_gam-class`](https://gpitt71.github.io/poissonsuperlearner/reference/Learner_gam-class.md)
  [`Learner_gam`](https://gpitt71.github.io/poissonsuperlearner/reference/Learner_gam-class.md)
  :

  GAM learner via [`mgcv::bam`](https://rdrr.io/pkg/mgcv/man/bam.html)

- [`Learner_hal-class`](https://gpitt71.github.io/poissonsuperlearner/reference/Learner_hal-class.md)
  [`Learner_hal`](https://gpitt71.github.io/poissonsuperlearner/reference/Learner_hal-class.md)
  : HAL learner for piecewise Poisson hazards

- [`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md)
  : Fit a single base learner

- [`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md)
  : Fit a Poisson Super Learner ensemble

## Predict and inspect

Generate predictions and inspect fitted models.

- [`predict(`*`<base_learner>`*`)`](https://gpitt71.github.io/poissonsuperlearner/reference/predict.base_learner.md)
  : Predict hazards, survival and absolute risk from a fitted base
  learner

- [`predict(`*`<poisson_superlearner>`*`)`](https://gpitt71.github.io/poissonsuperlearner/reference/predict.poisson_superlearner.md)
  : Predict hazards, survival and absolute risk from a fitted Poisson
  Super Learner

- [`predictRisk(`*`<base_learner>`*`)`](https://gpitt71.github.io/poissonsuperlearner/reference/predictRisk.base_learner.md)
  : Absolute-risk matrix predictions for a fitted base learner

- [`predictRisk(`*`<poisson_superlearner>`*`)`](https://gpitt71.github.io/poissonsuperlearner/reference/predictRisk.poisson_superlearner.md)
  : Absolute-risk matrix predictions for a fitted Poisson Super Learner

- [`coef(`*`<base_learner>`*`)`](https://gpitt71.github.io/poissonsuperlearner/reference/coef.base_learner.md)
  : Extract coefficients from a fitted base learner

- [`coef(`*`<poisson_superlearner>`*`)`](https://gpitt71.github.io/poissonsuperlearner/reference/coef.poisson_superlearner.md)
  : Extract stacking (meta-learner) coefficients from a fitted Poisson
  Super Learner

- [`print(`*`<base_learner>`*`)`](https://gpitt71.github.io/poissonsuperlearner/reference/print.base_learner.md)
  :

  Print method for `base_learner`

- [`print(`*`<poisson_superlearner>`*`)`](https://gpitt71.github.io/poissonsuperlearner/reference/print.poisson_superlearner.md)
  :

  Print method for `poisson_superlearner`

- [`summary(`*`<base_learner>`*`)`](https://gpitt71.github.io/poissonsuperlearner/reference/summary.base_learner.md)
  : Summarize a fitted base learner object

- [`summary(`*`<poisson_superlearner>`*`)`](https://gpitt71.github.io/poissonsuperlearner/reference/summary.poisson_superlearner.md)
  : Summarize a fitted Poisson Super Learner object

## Hazard and risk utilities

- [`pch_survival()`](https://gpitt71.github.io/poissonsuperlearner/reference/pch_survival.md)
  : Piecewise-constant hazards survival function
- [`pch_absolute_risk()`](https://gpitt71.github.io/poissonsuperlearner/reference/pch_absolute_risk.md)
  : Absolute risk (cumulative incidence) for a cause under
  piecewise-constant hazards
- [`pch_absolute_risk_euler()`](https://gpitt71.github.io/poissonsuperlearner/reference/pch_absolute_risk_euler.md)
  : Absolute risk (Euler approximation) for a cause under
  piecewise-constant hazards

## Example data

- [`simulateStenoT1()`](https://gpitt71.github.io/poissonsuperlearner/reference/simulateStenoT1.md)
  : Simulate time-to-event data for hypothetical type-1 diabetes
  patients

## Package

- [`poissonsuperlearner`](https://gpitt71.github.io/poissonsuperlearner/reference/poissonsuperlearner-package.md)
  [`poissonsuperlearner-package`](https://gpitt71.github.io/poissonsuperlearner/reference/poissonsuperlearner-package.md)
  : poissonsuperlearner: Poisson Super Learner
