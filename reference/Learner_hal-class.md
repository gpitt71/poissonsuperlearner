# HAL learner for piecewise Poisson hazards

`Learner_hal` is a Reference Class implementing the learner interface
used by
[`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md)
and
[`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md).

## Details

**User-facing API:** users should **only initialize** the learner and
pass it to
[`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md)
/
[`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md).
The remaining methods documented below are part of the internal learner
interface and are **not meant to be called directly by users**.

**Wrapper role:** this class provides a piecewise-constant hazard
wrapper around a HAL-style indicator-basis construction, estimated by
L1-penalized Poisson regression using a `glmnet` backend. The
package-specific contribution is to (i) construct the long-format
Poisson representation with offsets for time at risk, (ii) generate
indicator bases compatible with piecewise hazards, and (iii) forward
backend fitting arguments supplied via `...`.

## Model

Let \\0=t_0 \< t_1 \< \cdots \< t_m\\ denote time knots and define
interval indicators \\I_k(t)=1\\t\in(t_k,t\_{k+1}\]\\\\. The HAL
piecewise-constant hazard model is \$\$ \lambda(t \mid x) =
\sum\_{k=0}^{m} I_k(t)\\\exp\\f(t,x)\\, \$\$ where \\f(t,x)\\ is
approximated by a finite linear combination of indicator basis
functions.

## Two-covariate illustration

Let \\x=(x_1,x_2)\\ be two covariates and let \\t_0 \< t_1 \< \cdots \<
t_R\\ be time grid points used to create step functions in time. Choose
covariate cutpoints \\c\_{1,1},\ldots,c\_{1,K_1}\\ for \\x_1\\ and
\\c\_{2,1},\ldots,c\_{2,K_2}\\ for \\x_2\\.

Define indicator bases: \$\$B_r(t) = 1\\t_r \le t\\\$\$ \$\$B\_{1,p}(x)
= 1\\c\_{1,p} \le x_1\\\$\$ \$\$B\_{2,q}(x) = 1\\c\_{2,q} \le x_2\\\$\$

A main-effects HAL approximation on the log-hazard scale can be written
as: \$\$ f\_\beta(t,x) = \beta_0 + \sum\_{r=1}^R \beta_r B_r(t) +
\sum\_{r=1}^R\sum\_{p=1}^{K_1} \beta\_{r,1,p} B_r(t) B\_{1,p}(x) +
\sum\_{r=1}^R\sum\_{q=1}^{K_2} \beta\_{r,2,q} B_r(t) B\_{2,q}(x). \$\$
If `max_degree >= 2`, the learner additionally includes interaction
bases such as \$\$ \sum\_{r=1}^R\sum\_{p=1}^{K_1}\sum\_{q=1}^{K_2}
\beta\_{r,12,pq} B_r(t) B\_{1,p}(x) B\_{2,q}(x). \$\$

## How reference class parameters map to the model

- `covariates`:

  Covariate columns used to build covariate indicator bases.

- `num_knots`:

  Controls the number of cutpoints per covariate used for indicator
  bases.

- `max_degree`:

  Maximum interaction order included in the basis expansion.

- `intercept`:

  Whether the backend penalized regression includes an intercept term.

- `cross_validation`:

  If `TRUE`, selects the penalty level using
  [`glmnet::cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html).

- `maxit_prefit`:

  Optional `maxit` value used for the initial HAL backend fit. Leave as
  `NA` to use the backend default.

- `fit_arguments`:

  Additional arguments forwarded to the `glmnet` backend (e.g.
  `nfolds`).

## Fields

- `covariates` (`character`):

  Names of covariate columns used in the basis.

- `cross_validation` (`logical`):

  Whether to use `cv.glmnet` to select the penalty.

- `intercept` (`logical`):

  Backend intercept flag.

- `max_degree` (`integer`):

  Maximum interaction order.

- `num_knots` (`numeric`):

  Knots used for basis construction.

- `lambda_opt` (`numeric`):

  Selected penalty level when using cross-validation.

- `maxit_prefit` (`numeric`):

  Optional `maxit` value used for the initial HAL backend fit.

- `fit_arguments` (`list`):

  Extra backend arguments forwarded to `glmnet`.

## Methods (internal learner interface)

- `initialize(...)`:

  Construct and configure the learner. This is the only method users
  should call.

- `hal_basis(...)`:

  Internal helper. Constructs HAL basis matrices and metadata for
  fitting.

- `hal_prepare_new(...)`:

  Internal helper. Builds prediction-time HAL basis matrices from fitted
  basis metadata.

- `private_fit(data, ...)`:

  Internal. Builds bases and fits the penalized Poisson model with
  offset `log(tij)`.

- `private_fit_all_causes(data, ...)`:

  Internal. Fits penalized Poisson HAL models for all requested causes
  using a shared basis setup.

- `private_predictor(model, newdata, ...)`:

  Internal. Evaluates the fitted approximation and returns hazards on
  the response scale.

## Examples

``` r
lrn <- Learner_hal(covariates = c("age", "sex"), max_degree = 2L, num_knots = c(10L, 5L))
```
