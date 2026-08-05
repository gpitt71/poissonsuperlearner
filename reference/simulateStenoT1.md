# Simulate time-to-event data for hypothetical type-1 diabetes patients

Simulate synthetic data inspired by the Steno Type-1 risk engine

## Usage

``` r
simulateStenoT1(
  n,
  coefficient_age = 0.05,
  coefficient_LDL = 0.1,
  value_diabetis = 0.02,
  seed = NULL,
  keep = NULL,
  scenario = c("alpha", "beta"),
  competing_risks = FALSE
)
```

## Arguments

- n:

  `numeric(1)`. Number of subjects to simulate.

- coefficient_age:

  `numeric(1)`. Log-hazard coefficient for age in the CVD model
  (`time.event.1`).

- coefficient_LDL:

  `numeric(1)`. Log-hazard coefficient for LDL in the CVD model
  (`time.event.1`).

- value_diabetis:

  `numeric(1)`. Log-hazard coefficient for diabetes duration in the CVD
  model (`time.event.1`).

- seed:

  `integer(1)` or `NULL`. Optional random seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before simulating
  the data. If `NULL`, the current RNG state is used.

- keep:

  `character` or `NULL`. Optional subset of columns to retain. If
  supplied, only those columns are returned.

- scenario:

  `character(1)`. One of `"alpha"` or `"beta"`. Scenario `"beta"`
  modifies the CVD hazard by adding nonlinear hinge-squared terms in age
  and LDL.

- competing_risks:

  `logical(1)`. If `TRUE`, simulates two event causes (CVD and death
  without CVD). Otherwise simulates CVD vs censoring.

## Value

A `data.table` with at least the following columns:

- id:

  `integer`. Subject identifier (1, ..., `n`).

- time_cvd:

  `numeric`. Observed follow-up time (minimum of event and censoring
  times; also includes competing risk time if `competing_risks = TRUE`
  in scenario `"alpha"`).

- status_cvd:

  `integer`. Observed event status: `0` = censored, `1` = CVD, and if
  `competing_risks = TRUE` in scenario `"alpha"`, `2` = death without
  prior CVD.

- time:

  `numeric`. Alias of `time_cvd` (kept for convenience).

- event:

  `integer`. Alias of `status_cvd` (kept for convenience).

- uncensored_time_cvd:

  `numeric`. Event time ignoring censoring (minimum of event causes
  only).

- uncensored_status_cvd:

  `integer`. Event cause ignoring censoring. In scenario `"alpha"` this
  is `1` (CVD) or `2` (death without CVD); in scenario `"beta"` this is
  always `1`.

- uncensored_time:

  `numeric`. Alias of `uncensored_time_cvd`.

- uncensored_event:

  `integer`. Alias of `uncensored_status_cvd`.

In addition, the returned table contains all baseline covariates listed
in **Details**. Internal latent variables used only for simulation are
removed before returning (e.g., log2 eGFR components and, in scenario
`"beta"`, the hinge-squared features).

## Details

Generates baseline covariates and event times for CVD and censoring,
with an optional competing-risks setting, for examples, benchmarks and
tests.

The simulator uses a structural equation model (via
[`lava::lvm`](https://kkholst.github.io/lava/reference/lvm.html)) to
generate realistic correlations between covariates. Event times are then
generated from cause-specific Weibull proportional hazards models, where
the linear predictor depends on the simulated covariates (and scenario).

The following baseline covariates are generated (column name, type,
interpretation):

- sex:

  `factor`. Binary sex indicator (generated Bernoulli, then stored as
  factor).

- age:

  `numeric`. Age at baseline (years).

- diabetes_duration:

  `numeric`. Duration of diabetes at baseline (years).

- value_SBP:

  `numeric`. Systolic blood pressure (SBP).

- value_LDL:

  `numeric`. LDL cholesterol.

- value_HBA1C:

  `numeric`. HbA1c.

- value_Albuminuria:

  `factor` with levels `Normal`, `Micro`, `Macro`. Albuminuria category.

- eGFR:

  `numeric`. Estimated glomerular filtration rate, constructed from
  latent age-dependent log2 eGFR components (higher values indicate
  better kidney function).

- value_Smoking:

  `factor`. Smoking indicator (generated from a logistic model, then
  stored as factor).

- value_Motion:

  `factor`. Physical activity indicator (generated from a logistic
  model, then stored as factor).

Event time variables are generated from latent Weibull PH models:
`time.event.1` (CVD), `time.event.0` (censoring), and, if
`competing_risks = TRUE`, `time.event.2` (death without prior CVD).
These latent variables are used to construct the observed outcome
variables returned by the function (see below).

## Author

Thomas A. Gerds <tag@biostat.ku.dk>

## Examples

``` r
simulateStenoT1(n = 20, scenario = "alpha", competing_risks = TRUE)
#>        id    sex      age diabetes_duration value_SBP value_LDL value_HBA1C
#>     <int> <fctr>    <num>             <num>     <num>     <num>       <num>
#>  1:     1      1 44.66874          21.43588  130.7987  3.529910    69.22495
#>  2:     2      0 33.96825          14.83521  127.6159  2.122889    63.59379
#>  3:     3      1 55.26420          26.66288  135.0435  2.848161    67.20565
#>  4:     4      1 48.06708          22.88815  131.2366  4.211202    67.09340
#>  5:     5      0 50.65504          23.46144  137.0843  2.466765    62.87616
#>  6:     6      1 42.52754          19.92477  128.7147  1.936915    63.70087
#>  7:     7      0 35.29577          17.70904  131.1816  1.884837    64.08542
#>  8:     8      1 46.21775          22.84393  129.7322  4.161552    69.79398
#>  9:     9      1 45.07253          22.92282  130.8554  3.194515    67.65601
#> 10:    10      0 35.26669          15.73616  130.0175  2.841777    64.69316
#> 11:    11      1 46.85717          22.69569  130.6185  2.761375    65.38951
#> 12:    12      1 27.75812          12.76577  119.4189  2.915309    68.87784
#> 13:    13      0 50.09513          22.53214  139.8852  4.624835    70.24427
#> 14:    14      1 45.79274          20.26410  130.4618  3.471899    69.36945
#> 15:    15      1 36.72486          18.88762  124.5884  3.304521    70.19088
#> 16:    16      0 37.62746          16.74453  130.0358  3.826638    70.34635
#> 17:    17      0 33.62173          14.72196  129.4750  1.977350    63.55149
#> 18:    18      0 45.20734          19.48598  134.4073  2.272251    63.99564
#> 19:    19      1 37.53680          19.20302  124.2261  1.741428    64.40549
#> 20:    20      1 29.46024          16.21203  120.0034  3.210640    70.50347
#>        id    sex      age diabetes_duration value_SBP value_LDL value_HBA1C
#>     <int> <fctr>    <num>             <num>     <num>     <num>       <num>
#>     value_Smoking value_Motion value_Albuminuria       eGFR time_event_1
#>            <fctr>       <fctr>            <fctr>      <num>        <num>
#>  1:             0            1             Micro 66.4025806    0.8064877
#>  2:             0            1            Normal  1.0280075    1.3440819
#>  3:             1            1             Micro  1.0013055    1.4352764
#>  4:             0            1             Micro  2.1531457    1.8217873
#>  5:             1            1            Normal  1.8323643    2.0455121
#>  6:             1            0             Macro  0.8860595    3.3951090
#>  7:             1            0            Normal 11.3830473    2.9884061
#>  8:             0            1            Normal  5.6534172    3.6098994
#>  9:             1            1            Normal  9.5790148    3.6629335
#> 10:             0            0            Normal  3.9371254   23.5155189
#> 11:             0            1            Normal 10.0279591   10.1855497
#> 12:             1            1            Normal  2.4370749   86.8539713
#> 13:             0            0            Normal  4.0919993    6.7803246
#> 14:             0            1             Micro 15.5627890    9.0773123
#> 15:             0            0            Normal  2.3105252   45.4914540
#> 16:             1            1            Normal 46.2020743   11.7161419
#> 17:             0            1            Normal  8.2270696   56.7607208
#> 18:             0            0            Normal 36.2222825   38.5529257
#> 19:             0            1            Normal  1.9420301   35.5325054
#> 20:             1            1            Normal  2.2552274   59.6181982
#>     value_Smoking value_Motion value_Albuminuria       eGFR time_event_1
#>            <fctr>       <fctr>            <fctr>      <num>        <num>
#>     time_event_0 time_event_2   time_cvd status_cvd uncensored_time_cvd
#>            <num>        <num>      <num>      <int>               <num>
#>  1:    19.999612     1.570617  0.8064877          1           0.8064877
#>  2:    16.735233    82.210661  1.3440819          1           1.3440819
#>  3:    22.435819     4.066679  1.4352764          1           1.4352764
#>  4:    14.640014    13.364725  1.8217873          1           1.8217873
#>  5:     9.912720    17.569942  2.0455121          1           2.0455121
#>  6:     9.217793     2.180126  2.1801262          2           2.1801262
#>  7:     6.446499    58.227819  2.9884061          1           2.9884061
#>  8:    10.595615    13.637685  3.6098994          1           3.6098994
#>  9:     4.650380    14.096915  3.6629335          1           3.6629335
#> 10:     4.810503   131.322662  4.8105031          0          23.5155189
#> 11:    11.109433     6.269747  6.2697466          2           6.2697466
#> 12:     6.577952    69.420007  6.5779521          0          69.4200065
#> 13:    17.956594    19.121171  6.7803246          1           6.7803246
#> 14:    14.940757    16.181654  9.0773123          1           9.0773123
#> 15:    11.262713    77.628623 11.2627134          0          45.4914540
#> 16:    11.568335    30.208000 11.5683349          0          11.7161419
#> 17:    14.733028    18.725193 14.7330279          0          18.7251927
#> 18:    16.024862    62.356943 16.0248620          0          38.5529257
#> 19:    17.814757    51.741907 17.8147570          0          35.5325054
#> 20:    17.940559    63.281092 17.9405586          0          59.6181982
#>     time_event_0 time_event_2   time_cvd status_cvd uncensored_time_cvd
#>            <num>        <num>      <num>      <int>               <num>
#>     uncensored_status_cvd       time event uncensored_time uncensored_event
#>                     <int>      <num> <int>           <num>            <int>
#>  1:                     1  0.8064877     1       0.8064877                1
#>  2:                     1  1.3440819     1       1.3440819                1
#>  3:                     1  1.4352764     1       1.4352764                1
#>  4:                     1  1.8217873     1       1.8217873                1
#>  5:                     1  2.0455121     1       2.0455121                1
#>  6:                     2  2.1801262     2       2.1801262                2
#>  7:                     1  2.9884061     1       2.9884061                1
#>  8:                     1  3.6098994     1       3.6098994                1
#>  9:                     1  3.6629335     1       3.6629335                1
#> 10:                     1  4.8105031     0      23.5155189                1
#> 11:                     2  6.2697466     2       6.2697466                2
#> 12:                     2  6.5779521     0      69.4200065                2
#> 13:                     1  6.7803246     1       6.7803246                1
#> 14:                     1  9.0773123     1       9.0773123                1
#> 15:                     1 11.2627134     0      45.4914540                1
#> 16:                     1 11.5683349     0      11.7161419                1
#> 17:                     2 14.7330279     0      18.7251927                2
#> 18:                     1 16.0248620     0      38.5529257                1
#> 19:                     1 17.8147570     0      35.5325054                1
#> 20:                     1 17.9405586     0      59.6181982                1
#>     uncensored_status_cvd       time event uncensored_time uncensored_event
#>                     <int>      <num> <int>           <num>            <int>
```
