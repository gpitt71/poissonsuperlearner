# Absolute-risk matrix predictions for a fitted Poisson Super Learner

S3 method compatible with
[`riskRegression::predictRisk`](https://rdrr.io/pkg/riskRegression/man/predictRisk.html)
returning one column per requested time.

## Usage

``` r
# S3 method for class 'poisson_superlearner'
predictRisk(object, newdata, times, cause = 1, model = "sl", ...)
```

## Arguments

- object:

  `poisson_superlearner`. Fitted object.

- newdata:

  `data.frame`. New covariate data.

- times:

  `numeric`. Prediction times.

- cause:

  `numeric(1)`. Cause index.

- model:

  Model selector. Default is `"sl"`. Allowed values are the same as in
  [`predict.poisson_superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/predict.poisson_superlearner.md),
  including `"discrete_sl"` and cause-specific vectors of base-learner
  labels or indices.

- ...:

  Unused.

## Value

`numeric` matrix with `nrow(newdata)` rows and `length(times)` columns.

## Examples

``` r
d <- simulateStenoT1(30, competing_risks = TRUE)

learners <- list(
  lasso = Learner_glmnet(
    covariates = "sex",
    alpha = 1,
    lambda = 0.01,
    cross_validation = FALSE
  ),
  ridge = Learner_glmnet(
    covariates = c("sex", "value_LDL"),
    alpha = 0,
    lambda = 0.01,
    cross_validation = FALSE
  )
)

fit <- Superlearner(
  data = d,
  id = "id",
  status = "status_cvd",
  event_time = "time_cvd",
  learners = learners,
  number_of_nodes = 3,
  nfold = 2
)

if (requireNamespace("riskRegression", quietly = TRUE)) {
  riskRegression::predictRisk(fit, newdata = d[1:3], times = c(1, 3), cause = 1)
}
#>            [,1]      [,2]
#> [1,] 0.09155551 0.2400111
#> [2,] 0.12180318 0.3137609
#> [3,] 0.11323322 0.2944297
```
