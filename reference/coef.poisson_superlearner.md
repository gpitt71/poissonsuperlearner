# Extract stacking (meta-learner) coefficients from a fitted Poisson Super Learner

Extracts the **meta-learner coefficients** (stacking weights) from a
fitted `poisson_superlearner` object returned by
[`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md).

## Usage

``` r
# S3 method for class 'poisson_superlearner'
coef(object, cause = NULL, model = "sl", ...)
```

## Arguments

- object:

  `poisson_superlearner`. A fitted ensemble returned by
  [`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md).

- cause:

  `numeric(1)` or `NULL`. Which cause to extract meta-learner
  coefficients for. If `NULL`, coefficients are returned for all causes.
  Causes are indexed `1, 2, ..., object$data_info$n_crisks`.

- model:

  Model selector. Default is `"sl"` for the stacked super learner.
  Allowed values are:

  `0`, `"sl"`, `"superlearner"`, or `"super_learner"`

  :   Extract coefficients from the stacked meta-learner. For causes
      with no fitted meta-learner, this falls back to the retained base
      learner.

  `"discrete_sl"` and aliases

  :   Extract coefficients from the cause-specific base learners with
      the smallest cross-validated deviance.

  learner label

  :   Extract coefficients from one stored base learner by its label in
      `object$data_info$learners_labels[[k]]`.

  `"learner_j"` or character integer `"j"`

  :   Extract coefficients from the `j`-th stored learner.

  integer `j >= 1`

  :   Extract coefficients from the `j`-th stored learner.

  vector of labels or positive integer indices

  :   Use cause-specific base learners; length must equal
      `object$data_info$n_crisks`.

- ...:

  Passed to the underlying [`coef()`](https://rdrr.io/r/stats/coef.html)
  method of the fitted meta-learner (learner-dependent; e.g., `s` for
  `glmnet`).

## Value

If `cause` is a single integer, returns the coefficient object produced
by [`coef()`](https://rdrr.io/r/stats/coef.html) for the selected
cause-specific fitted model: the meta-learner when `model = "sl"` and a
meta-learner is available, or the selected base learner when `model`
selects a base learner or no meta-learner is available.

If `cause = NULL`, returns a `list` of length
`object$data_info$n_crisks`, where element `[[k]]` contains coefficients
for the selected model for cause `k`.

If no fitted ensemble is present (`object$superlearner` is `NULL`),
signals a message and returns `invisible(object)`.

## Details

For each cause `k`, the ensemble stores a fitted meta-learner in
`object$superlearner[[k]]$meta_learner_fit`. This method dispatches to
the underlying [`coef()`](https://rdrr.io/r/stats/coef.html) method for
that fitted meta-learner.

**What coefficients represent.** These coefficients correspond to the
meta-learner regression of the outcome on the cross-validated
base-learner predictions (`Z1`, `Z2`, ...). Under the default
meta-learner, they are the stacking weights (on the scale defined by the
meta-learner).

**Learner-dependent output.** The returned coefficient object depends on
the meta-learner implementation (by default a `glmnet` fit, often
returning a sparse matrix). This method does not rename `Z*` terms or
post-process coefficients; it returns the output of
`coef(object$superlearner[[k]]$meta_learner_fit, ...)` unchanged.

**Single-learner special case.** If the ensemble was fit with only one
base learner, no meta-learner is fit and `meta_learner_fit` is `NULL`.
In that case, [`coef()`](https://rdrr.io/r/stats/coef.html) for the
`poisson_superlearner` does not have meta-learner coefficients to
return.

## Examples

``` r
d <- simulateStenoT1(50, competing_risks = TRUE)
learners <- list(
  glm = Learner_glmnet(covariates = c("age", "value_LDL"), lambda = 0, cross_validation = FALSE),
  gam = Learner_gam(covariates = c("age", "value_LDL"))
)
fit <- Superlearner(d, id="id", status="status_cvd", event_time="time_cvd",
                    learners=learners, number_of_nodes=4, nfold=2)

# meta-learner coefficients (cause 1)
coef(fit, cause = 1)
#> 3 x 1 sparse Matrix of class "dgCMatrix"
#>                     s0
#> (Intercept) .         
#> Z1          0.90743467
#> Z2          0.04019941

# meta-learner coefficients for all causes (list)
coef(fit)
#> $cause_1
#> 3 x 1 sparse Matrix of class "dgCMatrix"
#>                     s0
#> (Intercept) .         
#> Z1          0.90743467
#> Z2          0.04019941
#> 
#> $cause_2
#> 3 x 1 sparse Matrix of class "dgCMatrix"
#>                     s0
#> (Intercept)  .        
#> Z1           1.4696812
#> Z2          -0.6260814
#> 
```
