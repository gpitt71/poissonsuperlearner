# Absolute risk (cumulative incidence) for a cause under piecewise-constant hazards

Computes, per row, the cumulative incidence function at the end of each
interval, grouped by `id`. The number of causes is inferred from the
number of columns in `haz`.

## Usage

``` r
pch_absolute_risk(id, dt, haz, cause_idx, one_based = TRUE, na_is_zero = FALSE)
```

## Arguments

- id:

  Integer vector. Sorted by `id` then time.

- dt:

  Numeric vector of interval lengths.

- haz:

  Numeric matrix (n x C) of cause-specific hazards per interval. Columns
  correspond to causes 1..C.

- cause_idx:

  Integer. Index of the cause of interest (1-based by default).

- one_based:

  Logical. If `TRUE`, `cause_idx` is 1-based. If `FALSE`, 0-based.

- na_is_zero:

  Logical. If `TRUE`, treat NA/Inf hazards as zero.

## Value

Numeric vector of cumulative incidence values at the end of each
interval.

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
pch_absolute_risk(id = id, dt = dt, haz = haz, cause_idx = 1)
#> [1] 0.09286135 0.24158123 0.04829013 0.13572326

```
