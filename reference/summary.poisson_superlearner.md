# Summarize a fitted Poisson Super Learner object

Prints:

1.  a compact description of the fitted ensemble,

2.  cross-validated deviances for base learners (when available),

3.  cause-specific meta-learner coefficients (stacking weights).

## Usage

``` r
# S3 method for class 'poisson_superlearner'
summary(object, cause = NULL, model = "sl", ...)
```

## Arguments

- object:

  `poisson_superlearner` returned by
  [`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md).

- cause:

  `numeric(1)` or `NULL`. Which cause’s meta-learner fit to print. If
  `NULL`, prints one line per cause (classes only) instead of printing
  the full fitted objects.

- model:

  Model selector. Default is `"sl"` for the stacked super learner.
  Allowed values are:

  `0`, `"sl"`, `"superlearner"`, or `"super_learner"`

  :   Summarize the stacked meta-learner. For causes with no fitted
      meta-learner, this falls back to the retained base learner.

  `"discrete_sl"` and aliases

  :   Summarize the cause-specific base learners with the smallest
      cross-validated deviance.

  learner label

  :   Summarize one stored base learner by its label in
      `object$data_info$learners_labels[[k]]`.

  `"learner_j"` or character integer `"j"`

  :   Summarize the `j`-th stored learner.

  integer `j >= 1`

  :   Summarize the `j`-th stored learner.

  vector of labels or positive integer indices

  :   Use cause-specific base learners; length must equal
      `object$data_info$n_crisks`.

- ...:

  Passed to the underlying [`coef()`](https://rdrr.io/r/stats/coef.html)
  method for the fitted meta-learner (learner-dependent; e.g. `s` for
  `glmnet`).

## Value

Invisibly returns a `list` with elements:

- cross_validation_deviance:

  `data.table` (or `NULL`).

- meta_coefficients:

  List of length `n_crisks` with cause-specific coefficient objects (or
  `NULL`).

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

s <- summary(fit, cause = 1)
#> Call:
#>   Superlearner(..., metalearner = glmnet::glmnet)
#> 
#> Fitted object:
#>   Class: poisson_superlearner
#>   Number of competing risks: 2 
#>   Number of folds: 2 
#>   Maximum follow-up: 15.82977 
#>   Number of nodes: 5 
#> 
#> Retained learners by cause:
#>   cause 1: lasso, ridge
#>   cause 2: lasso, ridge
#> 
#> Cross-validation deviance:
#>    cause_index   cause learner_index learner  deviance
#>          <int>  <char>         <int>  <char>     <num>
#> 1:           1 cause_1             1   lasso 111.08762
#> 2:           1 cause_1             2   ridge 217.13423
#> 3:           2 cause_2             1   lasso  85.06296
#> 4:           2 cause_2             2   ridge  81.83562
#> 
#> Meta-learner coefficients:
#>   cause 1:
#> (Intercept)       lasso       ridge 
#>   0.0000000   1.7940042  -0.8157702 
#> 
names(s)
#> [1] "cross_validation_deviance" "meta_coefficients"        

```
