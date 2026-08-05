# GAM learner via `mgcv::bam`

`Learner_gam` is a Reference Class implementing the learner interface
used by
[`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md)
and
[`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md).

## Arguments

- covariates:

  `character`. Right-hand-side terms, including `mgcv` smooths (e.g.
  `"s(age)"`) and/or linear terms (e.g. `"value_LDL"`).

- cross_validation:

  `logical`. Included for compatibility with the learner interface;
  smoothing selection is controlled by `mgcv` and arguments in `...`.

## Details

**User-facing API:** users should **only initialize** the learner and
pass it to
[`Superlearner()`](https://gpitt71.github.io/poissonsuperlearner/reference/Superlearner.md)
/
[`fit_learner()`](https://gpitt71.github.io/poissonsuperlearner/reference/fit_learner.md).
The remaining methods documented below are part of the internal learner
interface and are **not meant to be called directly by users**.

**Wrapper role:** this class wraps
[`mgcv::bam`](https://rdrr.io/pkg/mgcv/man/bam.html) in a
piecewise-constant hazard workflow. The package-specific contribution is
to provide a convenient interface for the long-format Poisson likelihood
with offsets for time at risk, and optional node terms encoding the
baseline hazard, while forwarding standard
[`mgcv::bam`](https://rdrr.io/pkg/mgcv/man/bam.html) arguments supplied
via `...`.

## Model

Let \\0=t_0 \< t_1 \< \cdots \< t_m\\ denote time knots and define
interval indicators \\I_k(t)=1\\t\in(t_k,t\_{k+1}\]\\\\. The
piecewise-constant hazard model with an additive predictor is \$\$
\lambda(t \mid x) = \sum\_{k=0}^{m} I_k(t)\\\exp\\\eta(x) + \gamma_k\\.
\$\$ The additive predictor \\\eta(x)\\ is constructed from `covariates`
(smooth terms such as `s(age)` and/or linear terms) and estimated by
`mgcv`.

## Fields

- `covariates` (`character`):

  Terms used to build the additive predictor (may include `s()` terms).

- `cross_validation` (`logical`):

  Workflow flag; see Details.

- `intercept` (`logical`):

  Whether to include an intercept.

- `formula` (`character`):

  Formula string passed to
  [`mgcv::bam`](https://rdrr.io/pkg/mgcv/man/bam.html).

- `learner` (`function`):

  Backend fitter ([`mgcv::bam`](https://rdrr.io/pkg/mgcv/man/bam.html)).

- `fit_arguments` (`list`):

  Additional arguments forwarded to
  [`mgcv::bam`](https://rdrr.io/pkg/mgcv/man/bam.html).

## Methods (internal learner interface)

- `initialize(...)`:

  Construct and configure the learner. This is the only method users
  should call.

- `private_fit(data, ...)`:

  Internal. Fits a Poisson GAM with offset `log(tij)` on long-format
  data.

- `private_fit_all_causes(data, ...)`:

  Internal. Fits cause-specific Poisson GAMs for all requested causes
  using a shared long-format data setup.

- `private_predictor(model, newdata, ...)`:

  Internal. Predicts hazards on the response scale.

## Examples

``` r
lrn <- Learner_gam(covariates = c("s(age)", "value_LDL"))
```
