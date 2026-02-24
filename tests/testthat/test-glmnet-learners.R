### test-glmnet-learners.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 12 2026 (06:30)
## Version:
## Last-Updated: feb 12 2026 (15:24)
##           By: Thomas Alexander Gerds
##     Update #: 18
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
library(testthat)
library(poissonsuperlearner)
library(data.table)
library(survival)
library(riskRegression)


test_that("data_pre_processing preserves time at risk and covariates", {

  skip_if_not_installed("riskRegression")
  skip_if_not_installed("data.table")

  {set.seed(42)

  # Synthetic data
  Xvars <- paste0("X", 1:10)
  d <- riskRegression::sampleData(
    n = 100,
    formula = ~ f(X1, 2) + f(X2, 0) + f(X3, 0) + f(X6, 0) + f(X7, 0) +
      f(X8, 0) + f(X9, 0) + f(X10, 0)
  )}

  d <- as.data.table(d)
  #make it survival data set
  d[,event:=as.numeric(event==1)]
  d[,id:=1:.N]

  id_col <- "id"
  time_col <- "time"
  status_col <- "event"

  # Sanity: unique id in input
  testthat::expect_equal(nrow(d), uniqueN(d[[id_col]]))

  # Define 10 time knots
  qs <- as.numeric(stats::quantile(d[[time_col]], probs = seq(0.1, 1, by = 0.1), na.rm = TRUE))
  nodes <- sort(unique(c(0,qs)))
  testthat::expect_true(length(nodes) >= 2)  # otherwise intervals are degenerate

  # Run preprocessing
  dt_fit <- data_pre_processing(
    data = d,
    id = id_col,
    status = status_col,
    event_time = time_col,
    nodes = nodes,
    predictions = FALSE,
    uncensored_01 = FALSE
  )

  dt_fit <- as.data.table(dt_fit)

  # ---------- 1) Time preservation ----------
  # For each subject: sum(tij) should equal the observed time
  time_check <- dt_fit[, .(sum_tij = sum(tij)), by = id][
    d[, .(id = get(id_col), obs_time = get(time_col))],
    on = "id"
  ]

  testthat::expect_true(all(is.finite(time_check$sum_tij)))
  testthat::expect_true(all(is.finite(time_check$obs_time)))

  testthat::expect_true(
    all(abs(time_check$sum_tij - time_check$obs_time) < 1e-10)
  )

  # ---------- 2) Covariate preservation ----------
  # Each covariate must be constant within id in the expanded data,
  # and equal to its value in the original data.
  for (v in Xvars) {
    if (!v %in% names(d)) next
    if (!v %in% names(dt_fit)) {
      testthat::fail(paste0("Covariate '", v, "' missing from preprocessed data"))
    }

    # constant within id
    varying <- dt_fit[, data.table::uniqueN(get(v)), by = id][V1 != 1]
    testthat::expect_equal(nrow(varying), 0)

    # matches original
    merged <- unique(dt_fit[, .(id, val_fit = get(v))])[
      d[, .(id = get(id_col), val_orig = get(v))],
      on = "id"
    ]

    # handle factors vs numeric gracefully
    if (is.factor(merged$val_fit) || is.character(merged$val_fit)) {
      testthat::expect_true(all(as.character(merged$val_fit) == as.character(merged$val_orig)))
    } else {
      testthat::expect_true(all(abs(merged$val_fit - merged$val_orig) < 1e-12))
    }
  }

  # ---------- 3) Structural checks ----------
  testthat::expect_true("tij" %in% names(dt_fit))
  testthat::expect_true("deltaij" %in% names(dt_fit))
  testthat::expect_true("node" %in% names(dt_fit))
  testthat::expect_true("k" %in% names(dt_fit))

  # Node levels match what preprocessing decided (predictions=FALSE uses dt_fit$node max)
  testthat::expect_true(is.factor(dt_fit$node))
  testthat::expect_true(is.factor(dt_fit$k))

})

test_that("Poisson piecewise-constant fit reproduces Cox coefficient (X1) - both Superlearner and fit_learner", {

  skip_if_not_installed("riskRegression")
  skip_if_not_installed("survival")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("data.table")

  library(data.table)

  {set.seed(42)
  d <- riskRegression::sampleData(
    n = 50,
    formula = ~ f(X1, 2)
  )
  d <- as.data.table(d)}

  # sampleData sometimes has competing risks in 'event' (0,1,2,...).
  # For this test, reduce to a standard survival endpoint:
  # event of interest = 1, everything else treated as censored.
  if (!("id" %in% names(d))) d[, id := .I]
  if (!("time" %in% names(d))) stop("sampleData did not produce a 'time' column; adjust the test.")
  if (!("event" %in% names(d))) stop("sampleData did not produce an 'event' column; adjust the test.")
  if (!("X1" %in% names(d))) stop("sampleData did not produce 'X1'; adjust the test.")

  d[, status := as.integer(event == 1L)]

  # ---- Cox model ----
  fit_cox <- survival::coxph(
    survival::Surv(time, status) ~ X1,
    data = d,
    ties = "breslow"
  )
  beta_cox <- unname(coef(fit_cox)[["X11"]])

  # ---- Poisson piecewise-constant hazard model via your infrastructure ----
  # Use all observed times as nodes by leaving nodes=NULL and number_of_nodes=NULL
  # Ensure unpenalized fit: lambda=0.
  learner <- Learner_glmnet(
    covariates = "X1",
    cross_validation = FALSE,
    intercept = FALSE,
    add_nodes = TRUE,
    penalise_nodes = FALSE,
    lambda = 0,
    alpha = 1
  )

  olcheck <- Superlearner(d,
                          learner = list(learner),
                          id = "id",
                          status = "status",
                          event_time = "time",
                          nodes = NULL,
                          number_of_nodes = NULL)

  fit_ps <- fit_learner(
    data = d,
    learner = learner,
    id = "id",
    status = "status",
    event_time = "time",
    nodes = NULL,
    number_of_nodes = NULL
  )

  # fit_ps$learner_fit is a list by competing-risk stratum k.
  # For survival (after status recode), k should be length 1.
  testthat::expect_true(inherits(fit_ps, "base_learner"))
  testthat::expect_true(length(fit_ps$learner_fit) >= 1)

  glmnet_fit <- fit_ps$learner_fit[[1]]

  # Extract coefficient for X1 from glmnet fit.
  # coef(glmnet_fit) returns a sparse matrix with rownames including "X1".
  beta_pois <- as.matrix(stats::coef(glmnet_fit))
  testthat::expect_true("X11" %in% rownames(beta_pois))
  beta_pois <- unname(beta_pois["X11", 1])


  ## superlearner w one learner
  beta_ol <- coef(olcheck$superlearner[[1]]$learners_fit)["X11", 1]

  # ---- compare ----
  # Numerical equality depends on ties + discretization at all observed times;
  # this should be close, but allow a small tolerance.
  testthat::expect_equal(beta_pois, beta_cox, tolerance = 1e-3)
  testthat::expect_equal(beta_ol, beta_cox, tolerance = 1e-3)
})




