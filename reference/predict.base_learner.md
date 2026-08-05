# Predict hazards, survival and absolute risk from a fitted base learner

Computes **cause-specific piecewise-constant hazards** (`pwch_k`), the
corresponding **survival function**, and **absolute risk** for a given
cause, at user-supplied prediction horizons `times`, using a fitted
`base_learner` object (single learner; no stacking).

## Usage

``` r
# S3 method for class 'base_learner'
predict(object, newdata, times, cause = 1, ...)
```

## Arguments

- object:

  `base_learner`. A fitted object returned by
  [`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md).
  It contains the learner specification in `object$model` and
  cause-specific fitted models in `object$learner_fit`.

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

Internally, `newdata` is expanded to a Cartesian product with `times`,
converted to long Poisson format on `object$data_info$nodes`, and the
fitted learner for each cause in `object$learner_fit` is used to predict
the cause-specific hazards. Survival and absolute risk are then computed
from the predicted hazards.

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
lrn <- Learner_glmnet(covariates = c("age", "value_LDL"), lambda = 0, cross_validation = FALSE)
bl <- fit_learner(d, learner = lrn, id="id", status="status_cvd", event_time="time_cvd",
                  number_of_nodes=8)
p <- predict(bl, newdata = d[1:5], times = c(0, 2, 5), cause = 1)
head(p)
#>       id    sex      age diabetes_duration value_SBP value_LDL value_HBA1C
#>    <int> <fctr>    <num>             <num>     <num>     <num>       <num>
#> 1:     1      0 50.75134          22.89707  137.1395  2.431660    63.32654
#> 2:     1      0 50.75134          22.89707  137.1395  2.431660    63.32654
#> 3:     1      0 50.75134          22.89707  137.1395  2.431660    63.32654
#> 4:     2      0 51.98570          24.67824  138.0197  2.349015    64.52782
#> 5:     2      0 51.98570          24.67824  138.0197  2.349015    64.52782
#> 6:     2      0 51.98570          24.67824  138.0197  2.349015    64.52782
#>    value_Smoking value_Motion value_Albuminuria      eGFR time_event_1
#>           <fctr>       <fctr>            <fctr>     <num>        <num>
#> 1:             0            1            Normal 7.3955254    0.3536507
#> 2:             0            1            Normal 7.3955254    0.3536507
#> 3:             0            1            Normal 7.3955254    0.3536507
#> 4:             0            0            Normal 0.2869973    0.8590847
#> 5:             0            0            Normal 0.2869973    0.8590847
#> 6:             0            0            Normal 0.2869973    0.8590847
#>    time_event_0 time_event_2 uncensored_time_cvd uncensored_status_cvd
#>           <num>        <num>               <num>                 <int>
#> 1:     14.58882    18.461209           0.3536507                     1
#> 2:     14.58882    18.461209           0.3536507                     1
#> 3:     14.58882    18.461209           0.3536507                     1
#> 4:     15.16226     5.201973           0.8590847                     1
#> 5:     15.16226     5.201973           0.8590847                     1
#> 6:     15.16226     5.201973           0.8590847                     1
#>         time event uncensored_time uncensored_event time_cvd    pwch_1
#>        <num> <int>           <num>            <int>    <num>     <num>
#> 1: 0.3536507     1       0.3536507                1        0 0.0000000
#> 2: 0.3536507     1       0.3536507                1        2 0.2678121
#> 3: 0.3536507     1       0.3536507                1        5 0.3629399
#> 4: 0.8590847     1       0.8590847                1        0 0.0000000
#> 5: 0.8590847     1       0.8590847                1        2 0.3699958
#> 6: 0.8590847     1       0.8590847                1        5 0.5014196
#>          pwch_2 survival_function absolute_risk
#>           <num>             <num>         <num>
#> 1: 0.000000e+00         1.0000000     0.0000000
#> 2: 1.731764e-05         0.5747452     0.4252038
#> 3: 4.498156e-01         0.1065957     0.7075264
#> 4: 0.000000e+00         1.0000000     0.0000000
#> 5: 2.523829e-05         0.4652658     0.5346641
#> 6: 6.555496e-01         0.0429922     0.7868675
```
