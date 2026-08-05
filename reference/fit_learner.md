# Fit a single base learner

Pre-processes subject-level time-to-event data into a long Poisson
format on a piecewise-constant time grid, then fits **one** initialized
learner object. For competing risks, a separate model is fit for each
event type (cause) using the standard cause-specific Poisson likelihood
on the long data.

## Usage

``` r
fit_learner(
  data,
  learner,
  id = "id",
  stratified_k_fold = FALSE,
  status = "status",
  event_time = NULL,
  number_of_nodes = NULL,
  nodes = NULL,
  variable_transformation = NULL,
  ...
)
```

## Arguments

- data:

  `data.frame`. Subject-level input data (one row per subject).

- learner:

  Reference-class learner object (e.g. from
  [`Learner_glmnet()`](https://gpitt71.github.io/poissonsuperlearner/reference/Learner_glmnet-class.md),
  [`Learner_hal()`](https://gpitt71.github.io/poissonsuperlearner/reference/Learner_hal-class.md)
  or
  [`Learner_gam()`](https://gpitt71.github.io/poissonsuperlearner/reference/Learner_gam-class.md)).
  Must implement a `$private_fit(dt_long)` method that fits the learner
  on long Poisson data for one cause.

- id:

  `character(1)`. Name of the subject identifier column. If not found in
  `data`, an `id` column is created automatically.

- stratified_k_fold:

  `logical(1)`. Reserved argument for future fold strategy. Currently
  ignored.

- status:

  `character(1)`. Name of the event-status column. Must be coded with
  `0` = censoring and `1,2,...,K` for event types (causes). If there is
  no `0` in `status`, the data are treated as uncensored.

- event_time:

  `character(1)`. Name of the event/censoring time column. Must be
  present in `data`.

- number_of_nodes:

  `numeric(1)` or `NULL`. If not `NULL`, constructs a quantile-based
  node grid with `number_of_nodes + 1` cut points (including endpoints),
  then adds `0` if missing.

- nodes:

  `numeric` or `NULL`. Explicit time-node grid (cut points). If
  supplied, `number_of_nodes` is ignored. `0` is added if missing. Nodes
  beyond `max(event_time)` are dropped.

- variable_transformation:

  `list`/`character`/`formula` or `NULL`. Optional transformations
  applied to the internally created long Poisson data before fitting
  (via `apply_transformations()`).

- ...:

  Additional arguments currently ignored.

## Value

An object of class `base_learner`, i.e. a named `list` with:

- model:

  The **learner object** that was fit (the input `learner`), stored for
  later prediction. This contains the learner specification (e.g.,
  covariates, tuning parameters).

- learner_fit:

  A `list` of fitted model objects, **one per cause**. Its length equals
  `data_info$n_crisks`. The list is created by splitting the internally
  pre-processed long data by cause indicator `k` and calling
  `model$private_fit()` on each split.

  - Names typically correspond to the cause labels `"1"`, `"2"`, ...,
    `"K"`.

  - Each element is **learner-dependent**: e.g. for `Learner_glmnet` it
    may be a `"glmnet"` (often wrapped, e.g. `"fishnet"`) fit; for other
    learners it will be whatever `$private_fit()` returns.

  - Each fitted object is trained on long Poisson data representing the
    piecewise-constant hazard for that cause across the node intervals.

- data_info:

  A `list` of bookkeeping information needed for prediction and
  interpretation:

  id

  :   Identifier column name used.

  status

  :   Status column name used.

  event_time

  :   Event/censoring time column name used.

  nodes

  :   Numeric vector of node cut points used for the piecewise grid
      (includes `0` and is sorted). These are the interval boundaries
      used in the long Poisson representation.

  maximum_followup

  :   `max(data[[event_time]])`.

  n_crisks

  :   Number of event types (causes) detected. If censoring is present
      (`0` in `status`), then `n_crisks = #unique(status) - 1`;
      otherwise `n_crisks = #unique(status)`.

  variable_transformation

  :   The transformation specification passed in
      `variable_transformation` (or `NULL`).

## Examples

``` r
d <- simulateStenoT1(50, competing_risks = TRUE)
lrn <- Learner_glmnet(covariates = c("age", "value_LDL"),
                      lambda = 0,
                      cross_validation = FALSE)
bl <- fit_learner(d,
                  learner = lrn,
                  id = "id",
                  status = "status_cvd",
                  event_time = "time_cvd",
                  number_of_nodes = 2)
```
