# Penalized Poisson learner via `glmnet`

`Learner_glmnet` is a Reference Class implementing the learner interface
used by
[`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md)
and
[`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md).

## Details

**User-facing API:** users are expected to **initialize** the learner
(i.e., call `Learner_glmnet(...)`) and pass the resulting object to
[`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md)
or
[`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md).
The remaining methods documented below are part of the internal learner
interface and are **not meant to be called directly by users**.

**Wrapper role:** this class is a user-friendly wrapper around the
existing `glmnet` implementation. The package-specific contribution is
to provide a piecewise-constant hazard workflow: create the long-format
Poisson data with offsets for time at risk, include interval ("node")
indicators for the baseline hazard, and forward standard `glmnet`
arguments supplied at initialization to the backend fitter.

## Model

Let \\0=t_0 \< t_1 \< \cdots \< t_m\\ denote time knots and define
interval indicators \\I_k(t)=1\\t\in(t_k,t\_{k+1}\]\\\\. The
piecewise-constant hazard model is \$\$ \lambda(t \mid x) =
\sum\_{k=0}^{m} I_k(t)\\\lambda_k(x), \qquad \lambda_k(x) =
\exp(\beta^\top x + \gamma_k). \$\$ Penalization is applied to the
regression coefficients through the `glmnet` elastic-net penalty. Node
(baseline) terms are given zero penalty by default; if this backend call
fails, the learner retries with a fully penalized design.

## Fields

- `covariates` (`character`):

  Names of covariate columns used in the model.

- `cross_validation` (`logical`):

  If `TRUE`, chooses `lambda` by
  [`glmnet::cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html).

- `intercept` (`logical`):

  Backend intercept flag; currently fixed to `TRUE` by the constructor.

- `lambda` (`numeric`):

  If `cross_validation=FALSE`, the `lambda` used in the final fit.

- `formula` (`character`):

  Formula string used to create the design matrix in long format.

- `learner` (`function`):

  Backend fitter
  ([`glmnet::glmnet`](https://glmnet.stanford.edu/reference/glmnet.html)
  or
  [`glmnet::cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html)).

- `fit_arguments` (`list`):

  Additional arguments forwarded to the backend fitter.

## Methods (internal learner interface)

- `initialize(...)`:

  Construct and configure the learner. This is the only method users
  should call.

- `private_fit(data, ...)`:

  Internal. Fits a Poisson model with offset `log(tij)` on long-format
  data.

- `private_fit_all_causes(data, ...)`:

  Internal. Fits cause-specific Poisson models for all requested causes
  using a shared long-format data setup.

- `private_predictor(model, newdata, ...)`:

  Internal. Predicts hazards on the response scale for long-format
  `newdata`.

## Examples

``` r
lrn <- Learner_glmnet(covariates = c("age", "sex"), alpha = 1, cross_validation = TRUE)
```
