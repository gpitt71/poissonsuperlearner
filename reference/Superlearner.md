# Fit a Poisson Super Learner ensemble

Fits an ensemble of cause-specific piecewise-constant hazard models
using a long-format Poisson representation and combines them through a
meta-learner (stacking).

## Usage

``` r
Superlearner(
  data,
  id = "id",
  status = "status",
  event_time = NULL,
  learners,
  number_of_nodes = NULL,
  nodes = NULL,
  variable_transformation = NULL,
  nfold = 3,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  `data.frame`. Subject-level input data, one row per subject.

- id:

  `character(1)`. Name of the subject identifier column. If missing, an
  `id` column is created automatically.

- status:

  `character(1)`. Name of the event-status column. It must be coded with
  `0` for censoring and `1, 2, ..., K` for event types. If there is no
  `0` in `status`, the data are treated as uncensored.

- event_time:

  `character(1)`. Name of the event or censoring time column.

- learners:

  `list`. Either a single learner library used for every cause, or a
  list of cause-specific learner libraries. A learner library is a named
  or unnamed list of initialized learner reference-class objects, for
  example
  [`Learner_glmnet()`](https://gpitt71.github.io/poissonsuperlearner/reference/Learner_glmnet-class.md),
  [`Learner_hal()`](https://gpitt71.github.io/poissonsuperlearner/reference/Learner_hal-class.md),
  or
  [`Learner_gam()`](https://gpitt71.github.io/poissonsuperlearner/reference/Learner_gam-class.md).
  Missing learner names are filled as `"learner_1"`, `"learner_2"`, and
  so on; missing cause-library names are filled as `"cause_1"`,
  `"cause_2"`, and so on. Each learner must implement
  `$private_fit(dt_long)` and `$private_predictor(model, newdata)`.

- number_of_nodes:

  `numeric(1)` or `NULL`. If not `NULL`, constructs a quantile-based
  node grid with `number_of_nodes + 1` cut points. Ignored when `nodes`
  is supplied.

- nodes:

  `numeric` or `NULL`. Explicit time-node grid. If supplied,
  `number_of_nodes` is ignored. `0` is added if missing, and nodes
  larger than `max(event_time)` are dropped.

- variable_transformation:

  Optional transformation specification passed to
  `apply_transformations()` on the internally created long-format data.

- nfold:

  `numeric(1)`. Number of folds for cross-validation stacking.

- verbose:

  logical(1). If TRUE, display progress bars during full-data fitting
  and cross-validation fitting. Defaults to FALSE.

- ...:

  Additional arguments currently ignored.

## Value

An object of class `poisson_superlearner`, stored as a named `list` with
the following components:

`learners`: a cause-specific list of retained base learner libraries.
Thus `learners[[k]][[j]]` is the `j`-th retained learner object for
cause `k`.

`metalearner`: a list describing the internal meta-learner used for
stacking (`engine = "glmnet::glmnet"`, Poisson family, no intercept,
`lambda = 0`, `add_nodes = FALSE`, log-hazard scale). If no stacking is
performed because only one learner remains for every cause,
`metalearner` is `NULL`.

`superlearner`: a `list` of length `data_info$n_crisks`, one entry per
cause. For cause `k`, `superlearner[[k]]` is a `list` with two elements:

- `learners_fit`: the fitted base learner object or objects for cause
  `k`. If more than one learner is retained, this is a `list` with one
  fitted object per retained learner. If only one learner remains, this
  is the single fitted learner object itself.

- `meta_learner_fit`: the fitted cause-specific meta-learner for cause
  `k`. If no stacking is performed, this is `NULL`.

`cross_validation_deviance`: a `data.table` with columns `cause_index`,
`cause`, `learner_index`, `learner`, and `deviance`, giving the
cross-validated Poisson deviance for each retained base learner within
each cause. This component is absent when all causes are fitted directly
with a single retained learner.

`data_info`: a `list` of bookkeeping information used for prediction and
interpretation, containing:

- `id`: identifier column name used.

- `status`: status column name used.

- `event_time`: event-time column name used.

- `nodes`: numeric vector of node cut points used for the piecewise
  grid.

- `nfold`: number of folds used for stacking.

- `maximum_followup`: maximum observed follow-up time.

- `n_crisks`: number of event types detected.

- `learners_labels`: list of character vectors with retained learner
  labels for each cause.

- `variable_transformation`: the transformation specification passed in
  `variable_transformation`, or `NULL`.

## Details

Internally, the function:

1.  builds a time grid (`nodes`) and converts the subject-level data to
    a long Poisson format;

2.  fits each base learner once on the full long data for each cause;

3.  removes learners that already fail on the full data;

4.  uses `nfold` cross-validation to obtain out-of-sample base-learner
    predictions (`Z1`, `Z2`, ...) for stacking;

5.  removes learners whose cross-validated prediction column is entirely
    missing for at least one cause;

6.  fits a cause-specific meta-learner on the retained stacked
    predictions.

If all learners fail on the full data, the function stops with an error.
If only one learner remains after the full-data screening step or after
the cross-validation screening step, no meta-learner is fit. In that
case, `metalearner` is `NULL`, each `superlearner[[k]]$meta_learner_fit`
is `NULL`, and prediction is based directly on the stored fitted base
learner. If some, but not all, causes retain only one learner after
screening, those causes are predicted directly while other causes may
still use a fitted meta-learner. Numeric learner positions always refer
to the learners actually retained for the corresponding cause in the
fitted object.

## Examples

``` r
data <- simulateStenoT1(50, competing_risks = TRUE)

learners <- list(
  glm = Learner_glmnet(
    covariates = c("sex", "value_LDL"),
    lambda = 0,
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
  data = data,
  id = "id",
  status = "status_cvd",
  event_time = "time_cvd",
  learners = learners,
  number_of_nodes = 3,
  nfold = 2
)
```
