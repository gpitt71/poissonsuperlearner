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
# tests/testthat/test-glmnet-learners.R

testthat::test_that("data_pre_processing preserves time at risk and covariates", {

  testthat::skip_if_not_installed("riskRegression")
  testthat::skip_if_not_installed("data.table")

  set.seed(42)

  Xvars <- paste0("X", 1:10)

  d <- riskRegression::sampleData(
    n = 100,
    formula = ~ f(X1, 2) + f(X2, 0) + f(X3, 0) + f(X6, 0) + f(X7, 0) +
      f(X8, 0) + f(X9, 0) + f(X10, 0)
  )
  d <- data.table::as.data.table(d)

  d[, event := as.integer(event == 1L)]
  d[, id := .I]

  id_col <- "id"
  time_col <- "time"
  status_col <- "event"

  testthat::expect_equal(nrow(d), data.table::uniqueN(d[[id_col]]))

  qs <- as.numeric(stats::quantile(d[[time_col]], probs = seq(0.1, 1, by = 0.1), na.rm = TRUE))
  nodes <- sort(unique(c(0, qs)))
  testthat::expect_true(length(nodes) >= 2)

  dt_fit <- poissonsuperlearner:::data_pre_processing(
    data = d,
    id = id_col,
    status = status_col,
    event_time = time_col,
    nodes = nodes,
    predictions = FALSE,
    uncensored_01 = FALSE
  )
  dt_fit <- data.table::as.data.table(dt_fit)

  # 1) time preservation
  time_check <- dt_fit[, .(sum_tij = sum(tij)), by = id][
    d[, .(id = get(id_col), obs_time = get(time_col))],
    on = "id"
  ]

  testthat::expect_true(all(is.finite(time_check$sum_tij)))
  testthat::expect_true(all(is.finite(time_check$obs_time)))

  testthat::expect_true(
    all(abs(time_check$sum_tij - time_check$obs_time) < 1e-8)
  )

  # 2) covariate preservation
  for (v in Xvars) {
    if (!v %in% names(d)) next

    if (!v %in% names(dt_fit)) {
      testthat::fail(paste0("Covariate '", v, "' missing from preprocessed data"))
    }

    varying <- dt_fit[, data.table::uniqueN(get(v)), by = id][V1 != 1]
    testthat::expect_equal(nrow(varying), 0)

    merged <- unique(dt_fit[, .(id, val_fit = get(v))])[
      d[, .(id = get(id_col), val_orig = get(v))],
      on = "id"
    ]

    if (is.factor(merged$val_fit) || is.character(merged$val_fit)) {
      testthat::expect_true(all(as.character(merged$val_fit) == as.character(merged$val_orig)))
    } else {
      testthat::expect_true(all(abs(merged$val_fit - merged$val_orig) < 1e-10))
    }
  }

  # 3) structure
  testthat::expect_true(all(c("tij", "deltaij", "node", "k") %in% names(dt_fit)))
  testthat::expect_true(is.factor(dt_fit$node))
  testthat::expect_true(is.factor(dt_fit$k))
})


testthat::test_that("Poisson piecewise-constant fit reproduces Cox coefficient (X1) - both Superlearner and fit_learner", {

  testthat::skip_if_not_installed("riskRegression")
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("glmnet")
  testthat::skip_if_not_installed("data.table")

  set.seed(42)

  d <- riskRegression::sampleData(n = 50, formula = ~ f(X1, 2))
  d <- data.table::as.data.table(d)

  if (!("id" %in% names(d))) d[, id := .I]
  if (!all(c("time", "event", "X1") %in% names(d))) {
    stop("sampleData did not produce expected columns: time, event, X1")
  }

  d[, status := as.integer(event == 1L)]

  fit_cox <- survival::coxph(survival::Surv(time, status) ~ X1, data = d, ties = "breslow")
  beta_cox_vec <- stats::coef(fit_cox)

  # robust name lookup
  if ("X1" %in% names(beta_cox_vec)) {
    beta_cox <- unname(beta_cox_vec[["X1"]])
  } else {
    ix <- grep("^X1", names(beta_cox_vec))
    testthat::expect_true(length(ix) >= 1)
    beta_cox <- unname(beta_cox_vec[ix[1]])
  }

  learner <- poissonsuperlearner::Learner_glmnet(
    covariates = "X1",
    cross_validation = FALSE,
    intercept = FALSE,
    add_nodes = TRUE,
    penalise_nodes = FALSE,
    lambda = 0,
    alpha = 1
  )

  olcheck <- poissonsuperlearner::Superlearner(
    data = d,
    learners = list(learner),
    id = "id",
    status = "status",
    event_time = "time",
    nodes = NULL,
    number_of_nodes = NULL
  )

  fit_ps <- poissonsuperlearner::fit_learner(
    data = d,
    learner = learner,
    id = "id",
    status = "status",
    event_time = "time",
    nodes = NULL,
    number_of_nodes = NULL
  )

  testthat::expect_true(inherits(fit_ps, "base_learner"))
  testthat::expect_true(length(fit_ps$learner_fit) >= 1)

  glmnet_fit <- fit_ps$learner_fit[[1]]
  beta_pois_mat <- as.matrix(stats::coef(glmnet_fit))

  # robust rowname lookup
  if ("X1" %in% rownames(beta_pois_mat)) {
    beta_pois <- unname(beta_pois_mat["X1", 1])
  } else {
    ix <- grep("^X1", rownames(beta_pois_mat))
    testthat::expect_true(length(ix) >= 1)
    beta_pois <- unname(beta_pois_mat[ix[1], 1])
  }

  beta_ol_mat <- as.matrix(stats::coef(olcheck$superlearner[[1]]$learners_fit))
  if ("X1" %in% rownames(beta_ol_mat)) {
    beta_ol <- unname(beta_ol_mat["X1", 1])
  } else {
    ix <- grep("^X1", rownames(beta_ol_mat))
    testthat::expect_true(length(ix) >= 1)
    beta_ol <- unname(beta_ol_mat[ix[1], 1])
  }

  testthat::expect_equal(beta_pois, beta_cox, tolerance = 1e-3)
  testthat::expect_equal(beta_ol,   beta_cox, tolerance = 1e-3)
})
