# Piecewise-constant hazards survival function

Computes survival at the end of each interval for competing risks with
piecewise constant hazards.

## Usage

``` r
pch_survival(id, dt, haz, na_is_zero = FALSE)
```

## Arguments

- id:

  Integer vector of subject IDs, sorted by id then time.

- dt:

  Numeric vector of interval lengths.

- haz:

  Numeric matrix (n x C) of cause-specific hazards.

- na_is_zero:

  Logical. If TRUE, treat NA hazards as zero.

## Value

Numeric vector of survival probabilities at the end of each interval.

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
pch_survival(id = id, dt = dt, haz = haz)
#> [1] 0.8607080 0.6376282 0.9323938 0.8187308
```
