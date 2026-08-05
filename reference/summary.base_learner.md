# Summarize a fitted base learner object

Dispatches to the underlying fitted model’s
[`summary()`](https://rdrr.io/r/base/summary.html) method for the
selected cause, or returns a list of summaries for all causes.

## Usage

``` r
# S3 method for class 'base_learner'
summary(object, cause = 1, ...)
```

## Arguments

- object:

  `base_learner` returned by
  [`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md).

- cause:

  `numeric(1)` or `NULL`. Which cause to summarize. If `NULL`, returns
  one summary per cause.

- ...:

  Passed to the underlying
  [`summary()`](https://rdrr.io/r/base/summary.html) method
  (learner-dependent).

## Value

If `cause` is a single integer, returns the underlying model summary for
that cause. If `cause = NULL`, returns a list of summaries (one per
cause).

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
out <- summary(bl, cause = 1)

```
