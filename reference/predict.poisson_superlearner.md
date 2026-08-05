# Predict hazards, survival and absolute risk from a fitted Poisson Super Learner

Computes **cause-specific piecewise-constant hazards** (`pwch_k`), the
corresponding **survival function**, and **absolute risk** for a given
cause, at user-supplied prediction horizons `times`, for each row in
`newdata`.

## Usage

``` r
# S3 method for class 'poisson_superlearner'
predict(object, newdata, times, cause = 1, model = "sl", ...)
```

## Arguments

- object:

  `poisson_superlearner`. A fitted ensemble from
  [`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md).

- newdata:

  `data.frame`/`data.table`. New covariate data (one row per subject).
  If `newdata` contains the original `event_time`, `status`, or `id`
  columns used for fitting, they are ignored for prediction.

- times:

  `numeric`. Prediction horizon(s). May include `0`. Times larger than
  `object$data_info$maximum_followup` are not supported: if **all**
  requested times exceed the maximum follow-up, a warning is issued and
  `NULL` is returned; if only **some** exceed, output rows for those
  times are returned with `NA` predictions.

- cause:

  `numeric(1)`. Cause index (1, 2, ...) used for the `absolute_risk`
  calculation.

- model:

  Model selector. Default is `"sl"` for the stacked super learner.
  Allowed values are:

  `0`, `"sl"`, `"superlearner"`, or `"super_learner"`

  :   Use the stacked super learner prediction. For causes with only one
      retained learner or no fitted meta-learner, this falls back to the
      retained base learner for that cause.

  `"discrete_sl"` and aliases

  :   For each cause, use the retained base learner with the smallest
      cross-validated deviance.

  learner label

  :   Use one stored base learner by its label in
      `object$data_info$learners_labels[[k]]`.

  `"learner_j"` or character integer `"j"`

  :   Use the `j`-th stored learner.

  integer `j >= 1`

  :   Use the `j`-th stored learner.

  vector of labels or positive integer indices

  :   Use cause-specific base learners; length must equal
      `object$data_info$n_crisks`.

  Numeric positions refer to the learners actually retained for each
  cause in the fitted object.

- ...:

  Additional arguments (currently ignored).

## Value

A `data.table` with one row per `(row in newdata, time in times)` and
columns:

- (original columns):

  All columns from `newdata` (excluding ignored event columns).

- time column:

  A column with name `object$data_info$event_time` holding the requested
  horizon.

- pwch_1, pwch_2, ...:

  Predicted cause-specific piecewise hazards at the horizon.

- survival_function:

  Predicted survival probability at the horizon.

- absolute_risk:

  Predicted cumulative incidence (absolute risk) for `cause` at the
  horizon.

## Details

Internally, `newdata` is expanded to a Cartesian product with the
requested `times`, converted to long Poisson format on
`object$data_info$nodes`, and hazards are predicted either from the
stacked super learner (`model = "sl"`), the discrete super learner
(`model = "discrete_sl"`), or selected fitted base learners. Survival
and absolute risk are then computed from the predicted hazards.

**Special case `times = 0`:** when `0` is included in `times`, the
returned rows have `survival_function = 1`, `absolute_risk = 0`, and all
`pwch_k = 0` at time 0.

**Identifiers in the output:** if `newdata` contains the `id` column, it
is carried into the output. If `newdata` does not contain an id column,
an internal id is created for computation, but it is not guaranteed to
appear in the returned table unless it was present in `newdata`.

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
p <- predict(fit, newdata = d[1:3], times = c(0, 2), cause = 1)
p[, .(id, time_cvd, absolute_risk)]
#>       id time_cvd absolute_risk
#>    <int>    <num>         <num>
#> 1:     1        0    0.00000000
#> 2:     1        2    0.05184506
#> 3:     2        0    0.00000000
#> 4:     2        2    0.09291898
#> 5:     3        0    0.00000000
#> 6:     3        2    0.07529414

```
