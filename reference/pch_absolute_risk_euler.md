# Absolute risk (Euler approximation) for a cause under piecewise-constant hazards

Computes the cumulative incidence function using the first-order Euler
(discrete) approximation: \$\$F_j(t) \approx \sum S(t\_{k-1})
\lambda\_{j,k} \Delta t_k\$\$ Grouped by `id`, this returns the
cumulative incidence at the end of each interval.

## Usage

``` r
pch_absolute_risk_euler(
  id,
  dt,
  haz,
  cause_idx,
  one_based = TRUE,
  na_is_zero = FALSE
)
```

## Arguments

- id:

  Integer vector. Sorted by `id` then time.

- dt:

  Numeric vector of interval lengths.

- haz:

  Numeric matrix (n x C) of cause-specific hazards per interval.

- cause_idx:

  Integer. Index of the cause of interest (1-based by default).

- one_based:

  Logical. If `TRUE`, `cause_idx` is 1-based. If `FALSE`, 0-based.

- na_is_zero:

  Logical. If `TRUE`, treat NA/Inf hazards as zero.

## Value

Numeric vector of cumulative incidence values (Euler approximation) at
the end of each interval.

## Examples

``` r
id <- c(1L, 1L, 2L, 2L)
dt <- c(1, 1, 1, 1)
haz <- rbind(
  c(0.10, 0.05),
  c(0.20, 0.10),
  c(0.05, 0.02),
  c(0.10, 0.03)
)
pch_absolute_risk_euler(id = id, dt = dt, haz = haz, cause_idx = 1)
#> [1] 0.1000000 0.2721416 0.0500000 0.1432394
```
