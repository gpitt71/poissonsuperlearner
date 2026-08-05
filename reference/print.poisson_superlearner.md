# Print method for `poisson_superlearner`

Prints a compact description of the fitted Poisson Super Learner,
including the number of base learners, the meta-learner, the time-grid
used, and competing-risk structure. Optionally prints the fitted
meta-learner for a given cause.

## Usage

``` r
# S3 method for class 'poisson_superlearner'
print(x, cause = 1, model = "sl", ...)
```

## Arguments

- x:

  `poisson_superlearner` object returned by
  [`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md).

- cause:

  `numeric(1)` or `NULL`. Which cause's meta-learner fit to print. If
  `NULL`, prints one line per cause (classes only) instead of printing
  the full fitted objects.

- model:

  Model selector. Default is `"sl"` for the stacked super learner.
  Allowed values are:

  `0`, `"sl"`, `"superlearner"`, or `"super_learner"`

  :   Print the stacked meta-learner. For causes with no fitted
      meta-learner, this falls back to the retained base learner.

  `"discrete_sl"` and aliases

  :   Print the cause-specific base learner with the smallest
      cross-validated deviance.

  learner label

  :   Print one stored base learner by its label in
      `x$data_info$learners_labels[[k]]`.

  `"learner_j"` or character integer `"j"`

  :   Print the `j`-th stored learner.

  integer `j >= 1`

  :   Print the `j`-th stored learner.

  vector of labels or positive integer indices

  :   Use cause-specific base learners; length must equal
      `x$data_info$n_crisks`.

- ...:

  Passed to the underlying fitted meta-learner
  [`print()`](https://rdrr.io/r/base/print.html) method when `cause` is
  a single integer.

## Value

Invisibly returns `x`.

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

print(fit, cause = NULL)
#> 
#> Cause 1, model = sl
#> 
#> Call:  glmnet::glmnet(x = x_meta, y = as.numeric(level_one_data[["deltaij"]][ok]),      family = "poisson", offset = log(level_one_data[["tij"]][ok]),      lambda = 0, standardize = FALSE, intercept = FALSE) 
#> 
#>   Df  %Dev Lambda
#> 1  2 82.33      0
#> 
#> Cause 2, model = sl
#> 
#> Call:  glmnet::glmnet(x = x_meta, y = as.numeric(level_one_data[["deltaij"]][ok]),      family = "poisson", offset = log(level_one_data[["tij"]][ok]),      lambda = 0, standardize = FALSE, intercept = FALSE) 
#> 
#>   Df %Dev Lambda
#> 1  2 92.2      0
```
