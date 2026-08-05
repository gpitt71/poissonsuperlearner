# Extract coefficients from a fitted base learner

Convenience method to extract (cause-specific) model coefficients from a
fitted `base_learner` returned by
[`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md).

## Usage

``` r
# S3 method for class 'base_learner'
coef(object, cause = NULL, ...)
```

## Arguments

- object:

  `base_learner`. A fitted object returned by
  [`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md).

- cause:

  `numeric(1)` or `NULL`. Which cause to extract coefficients for. If
  `NULL`, coefficients are returned for all causes. Causes are indexed
  `1, 2, ..., object$data_info$n_crisks` (with `0` reserved for
  censoring).

- ...:

  Passed to the underlying [`coef()`](https://rdrr.io/r/stats/coef.html)
  method of the fitted learner object (learner-dependent; e.g., `s` for
  `glmnet`).

## Value

If `cause` is a single integer, returns the coefficient object produced
by [`coef()`](https://rdrr.io/r/stats/coef.html) for that cause-specific
fitted model.

If `cause = NULL`, returns a `list` of length
`object$data_info$n_crisks`, where element `[[k]]` contains coefficients
for cause `k`.

If no fitted model is present (`object$learner_fit` is `NULL`), signals
a message and returns `invisible(object)`.

## Details

For competing risks,
[`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md)
fits one model per cause, stored in `object$learner_fit[[k]]` for
`k = 1, 2, ..., K`. This method simply dispatches to the underlying
model’s [`coef()`](https://rdrr.io/r/stats/coef.html) method for each
fitted object.

**Learner-dependent output.** The returned coefficient object depends on
the base learner used (e.g. a numeric vector, a sparse matrix, a list,
etc.). This method does not post-process or rename coefficients; it
returns the output of `coef(object$learner_fit[[k]], ...)` unchanged.

## Examples

``` r
d <- simulateStenoT1(30, competing_risks = TRUE)
lrn <- Learner_glmnet(covariates = c("age", "value_LDL"),
                      lambda = 0, cross_validation = FALSE)
bl <- fit_learner(d, learner = lrn, id = "id",
                  status = "status_cvd", event_time = "time_cvd",
                  number_of_nodes = 4)

# coefficients for cause 1
coef(bl, cause = 1)
#> 8 x 1 sparse Matrix of class "dgCMatrix"
#>                      s0
#> (Intercept) -7.27079032
#> age          0.07701094
#> value_LDL    0.20879304
#> noden2       1.04966630
#> noden3       0.98499689
#> noden4       0.60121570
#> noden5      -7.09418083
#> noden6       .         

# coefficients for all causes (list)
coef(bl)
#> $`1`
#> 8 x 1 sparse Matrix of class "dgCMatrix"
#>                      s0
#> (Intercept) -7.27079032
#> age          0.07701094
#> value_LDL    0.20879304
#> noden2       1.04966630
#> noden3       0.98499689
#> noden4       0.60121570
#> noden5      -7.09418083
#> noden6       .         
#> 
#> $`2`
#> 8 x 1 sparse Matrix of class "dgCMatrix"
#>                      s0
#> (Intercept) -296.632806
#> age            6.003270
#> value_LDL    -27.179372
#> noden2        -2.334743
#> noden3        -1.288138
#> noden4        61.936066
#> noden5        57.543521
#> noden6         .       
#> 
```
