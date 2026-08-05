# Basic Usage

## Introduction

This vignette gives fast examples of the Poisson Super Learner workflow
after the refactoring. It focuses on:

- fitting
  [`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md)
  with one learner library shared by all causes;
- fitting
  [`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md)
  with one learner library per cause;
- using [`summary()`](https://rdrr.io/r/base/summary.html) for the
  fitted super learner;
- using
  [`predictRisk()`](https://rdrr.io/pkg/riskRegression/man/predictRisk.html)
  for the same model selectors.

The examples use small simulated data, two folds, and simple `glmnet`
learners with fixed `lambda` values so the vignette remains quick to run
during package checks.

## Data

We simulate a small competing-risks data set. The observed follow-up
time is stored in `time` and the event indicator in `event`, where `0`
denotes censoring, `1` denotes cardiovascular disease, and `2` denotes
death without prior cardiovascular disease.

``` r

d <- simulateStenoT1(
  n = 45,
  scenario = "alpha",
  competing_risks = TRUE,
  seed = 1
)

d <- d[, .(
  id,
  time,
  event,
  sex,
  age,
  diabetes_duration,
  value_LDL,
  value_Smoking
)]

head(d)
#>       id      time event    sex      age diabetes_duration value_LDL
#>    <int>     <num> <int> <fctr>    <num>             <num>     <num>
#> 1:     1 0.1338473     2      1 36.03106          16.92239 3.9841195
#> 2:     2 0.1943785     1      0 42.33831          17.85834 1.5088089
#> 3:     3 0.5511647     1      0 43.04253          19.88147 0.9438390
#> 4:     4 0.6826765     2      0 45.93286          21.94420 2.9557458
#> 5:     5 1.1848861     1      1 45.95141          21.45525 4.6480322
#> 6:     6 1.2395837     1      0 39.59863          17.54727 0.8408968
#>    value_Smoking
#>           <fctr>
#> 1:             0
#> 2:             0
#> 3:             1
#> 4:             1
#> 5:             0
#> 6:             0
```

## One Learner Library For All Causes

A learner library is a list of initialized learner objects. If a single
library is supplied in a competing-risks analysis, the same library is
used for all causes.

``` r

shared_library <- list(
  simple = Learner_glmnet(
    covariates = c("sex", "diabetes_duration"),
    cross_validation = FALSE,
    lambda = 0
  ),
  shrink = Learner_glmnet(
    covariates = c("sex", "age", "value_LDL"),
    cross_validation = FALSE,
    lambda = 0.05,
    alpha = 1
  )
)

fit_shared <- Superlearner(
  data = d,
  id = "id",
  status = "event",
  event_time = "time",
  learners = shared_library,
  number_of_nodes = 3,
  nfold = 2
)
```

[`summary()`](https://rdrr.io/r/base/summary.html) gives a compact
overview of the fitted super learner, including the number of causes,
retained learner labels, cross-validated deviances, and meta-learner
coefficients when a meta-learner was fitted.

``` r

summary(fit_shared)
#> Call:
#>   Superlearner(..., metalearner = glmnet::glmnet)
#> 
#> Fitted object:
#>   Class: poisson_superlearner
#>   Number of competing risks: 2 
#>   Number of folds: 2 
#>   Maximum follow-up: 22.43464 
#>   Number of nodes: 5 
#> 
#> Retained learners by cause:
#>   cause 1: simple, shrink
#>   cause 2: simple, shrink
#> 
#> Cross-validation deviance:
#>    cause_index   cause learner_index learner  deviance
#>          <int>  <char>         <int>  <char>     <num>
#> 1:           1 cause_1             1  simple 110.34291
#> 2:           1 cause_1             2  shrink  97.91293
#> 3:           2 cause_2             1  simple 318.70522
#> 4:           2 cause_2             2  shrink  53.93232
#> 
#> Meta-learner coefficients:
#>   cause 1:
#> (Intercept)      simple      shrink 
#>   0.0000000   0.2277135   0.7138146 
#> 
#>   cause 2:
#> (Intercept)      simple      shrink 
#>   0.0000000   0.1989269   0.8821621
```

[`predictRisk()`](https://rdrr.io/pkg/riskRegression/man/predictRisk.html)
returns one row per subject and one column per requested prediction
time. The `model` argument uses the same selectors: `"sl"` for the
stacked ensemble, `"discrete_sl"` for the best cross-validated base
learner per cause, or learner labels such as `"simple"` and `"shrink"`
for models stored in the ensemble.

``` r

newdata <- d[1:2]
times <- c(1, 2)

risk_shared_sl <- predictRisk(
  fit_shared, newdata = newdata, times = times, cause = 1, model = "sl"
)

risk_shared_discrete <- predictRisk(
  fit_shared, newdata = newdata, times = times, cause = 1, model = "discrete_sl"
)

risk_shared_simple <- predictRisk(
  fit_shared, newdata = newdata, times = times, cause = 1, model = "simple"
)

risk_shared_shrink <- predictRisk(
  fit_shared, newdata = newdata, times = times, cause = 1, model = "shrink"
)

list(
  sl = risk_shared_sl,
  discrete_sl = risk_shared_discrete,
  simple = risk_shared_simple,
  shrink = risk_shared_shrink
)
#> $sl
#>            [,1]       [,2]
#> [1,] 0.03834563 0.07789423
#> [2,] 0.06632733 0.13279164
#> 
#> $discrete_sl
#>            [,1]       [,2]
#> [1,] 0.04304514 0.08340671
#> [2,] 0.05698820 0.10963701
#> 
#> $simple
#>            [,1]       [,2]
#> [1,] 0.01419019 0.02941751
#> [2,] 0.06819054 0.13891238
#> 
#> $shrink
#>            [,1]       [,2]
#> [1,] 0.04304514 0.08340671
#> [2,] 0.05698820 0.10963701
```

## One Learner Library Per Cause

For competing risks, `learners` can also be a list with one learner
library per cause. This allows different covariates, tuning parameters,
or labels for each cause.

``` r

libraries_by_cause <- list(
  cvd = list(
    cvd_simple = Learner_glmnet(
      covariates = c("sex", "diabetes_duration"),
      cross_validation = FALSE,
      lambda = 0
    ),
    cvd_shrink = Learner_glmnet(
      covariates = c("age", "value_LDL"),
      cross_validation = FALSE,
      lambda = 0.05,
      alpha = 1
    )
  ),
  death = list(
    death_simple = Learner_glmnet(
      covariates = c("sex", "age"),
      cross_validation = FALSE,
      lambda = 0
    ),
    death_shrink = Learner_glmnet(
      covariates = c("diabetes_duration", "value_Smoking"),
      cross_validation = FALSE,
      lambda = 0.05,
      alpha = 1
    )
  )
)

fit_by_cause <- Superlearner(
  data = d,
  id = "id",
  status = "event",
  event_time = "time",
  learners = libraries_by_cause,
  number_of_nodes = 3,
  nfold = 2
)
```

The fitted object can be summarized in the same way.

``` r

summary(fit_by_cause)
#> Call:
#>   Superlearner(..., metalearner = glmnet::glmnet)
#> 
#> Fitted object:
#>   Class: poisson_superlearner
#>   Number of competing risks: 2 
#>   Number of folds: 2 
#>   Maximum follow-up: 22.43464 
#>   Number of nodes: 5 
#> 
#> Retained learners by cause:
#>   cause 1: cvd_simple, cvd_shrink
#>   cause 2: death_simple, death_shrink
#> 
#> Cross-validation deviance:
#>    cause_index  cause learner_index      learner  deviance
#>          <int> <char>         <int>       <char>     <num>
#> 1:           1    cvd             1   cvd_simple 140.62531
#> 2:           1    cvd             2   cvd_shrink 100.14039
#> 3:           2  death             1 death_simple 156.00476
#> 4:           2  death             2 death_shrink  68.68571
#> 
#> Meta-learner coefficients:
#>   cause 1:
#> (Intercept)  cvd_simple  cvd_shrink 
#>   0.0000000  -0.0914743   1.0845723 
#> 
#>   cause 2:
#>  (Intercept) death_simple death_shrink 
#>    0.0000000   -0.1778133    1.3399726
```

The stacked and discrete super learner selectors still work as scalar
model selectors for prediction.

``` r

risk_by_cause_sl <- predictRisk(
  fit_by_cause, newdata = newdata, times = times, cause = 1, model = "sl"
)

risk_by_cause_discrete <- predictRisk(
  fit_by_cause, newdata = newdata, times = times, cause = 1, model = "discrete_sl"
)

list(
  sl = risk_by_cause_sl,
  discrete_sl = risk_by_cause_discrete
)
#> $sl
#>            [,1]       [,2]
#> [1,] 0.05189313 0.09760984
#> [2,] 0.05963114 0.11169165
#> 
#> $discrete_sl
#>            [,1]       [,2]
#> [1,] 0.04359312 0.08444487
#> [2,] 0.05638837 0.10851651
```

When selecting learners from cause-specific libraries, provide one
selector per cause. The first entry selects the learner for cause 1 and
the second entry selects the learner for cause 2.

``` r

cause_specific_model <- c("cvd_simple", "death_shrink")
cause_specific_model_alt <- c("cvd_shrink", "death_simple")

risk_by_cause_selected <- predictRisk(
  fit_by_cause, newdata = newdata, times = times, cause = 1,
  model = cause_specific_model
)

risk_by_cause_selected_alt <- predictRisk(
  fit_by_cause, newdata = newdata, times = times, cause = 1,
  model = cause_specific_model_alt
)

list(
  selected_learners = risk_by_cause_selected,
  selected_learners_alt = risk_by_cause_selected_alt
)
#> $selected_learners
#>            [,1]       [,2]
#> [1,] 0.01505978 0.03181764
#> [2,] 0.06979854 0.14278607
#> 
#> $selected_learners_alt
#>            [,1]     [,2]
#> [1,] 0.04175019 0.079684
#> [2,] 0.05531995 0.106062
```

Integer selectors are also supported. For cause-specific libraries, an
integer vector can choose different retained learner positions for
different causes.

``` r

risk_by_cause_indexed <- predictRisk(
  fit_by_cause, newdata = newdata, times = times, cause = 1, model = c(1, 2)
)

risk_by_cause_indexed
#>            [,1]       [,2]
#> [1,] 0.01505978 0.03181764
#> [2,] 0.06979854 0.14278607
```
