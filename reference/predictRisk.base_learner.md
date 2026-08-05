# Absolute-risk matrix predictions for a fitted base learner

Absolute-risk matrix predictions for a fitted base learner

## Usage

``` r
# S3 method for class 'base_learner'
predictRisk(object, newdata, times, cause = 1, ...)
```

## Arguments

- object:

  `base_learner`. Fitted object from
  [`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md).

- newdata:

  `data.frame`. New covariate data.

- times:

  `numeric`. Prediction times.

- cause:

  `numeric(1)`. Cause index.

- ...:

  Unused.

## Value

`numeric` matrix with `nrow(newdata)` rows and `length(times)` columns.

## Examples

``` r
d <- simulateStenoT1(30, competing_risks = TRUE)
lrn <- Learner_glmnet(
  covariates = c("sex", "value_LDL"),
  lambda = 0.01,
  cross_validation = FALSE
)
bl <- fit_learner(
  d,
  learner = lrn,
  id = "id",
  status = "status_cvd",
  event_time = "time_cvd",
  number_of_nodes = 3
)

if (requireNamespace("riskRegression", quietly = TRUE)) {
  riskRegression::predictRisk(bl, newdata = d[1:3], times = c(1, 3), cause = 1)
}
#>            [,1]      [,2]
#> [1,] 0.07174035 0.1992843
#> [2,] 0.14741794 0.3796804
#> [3,] 0.09266836 0.2521300

```
