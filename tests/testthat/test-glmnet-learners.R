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

testthat::test_that("current compressed preprocessing preserves time at risk and covariates", {

  testthat::skip_if_not_installed("riskRegression")
  testthat::skip_if_not_installed("data.table")

  set.seed(42)

  Xvars <- paste0("X", 1:10)

  d <- riskRegression::sampleData(
    n = 100,
    formula = ~ f(X1, 2) + f(X2, 0) + f(X3, 0) + f(X1, 0) + f(X7, 0) +
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

  grid_dt <- data.table::data.table(node = nodes)
  grid_dt[, prev_node_psl := data.table::shift(node)]
  grid_dt[, (time_col) := node]

  dt_terminal <- grid_dt[d, on = time_col, roll = Inf]

  dt_terminal[
    ,
    node_start := data.table::fifelse(
      get(time_col) == node,
      prev_node_psl,
      node
    )
  ]

  dt_terminal <- dt_terminal[
    !is.na(node_start),
    c("node", "tij") := list(node_start, get(time_col) - node_start)
  ]

  data.table::setnames(dt_terminal, old = status_col, new = "deltaij")

  # The compressed Superlearner()/fit_learner() preprocessing keeps the terminal
  # interval only. The learner expands those terminal rows back to the complete
  # piecewise-constant grid after grouping by covariate pattern.
  train_d <- data.table::copy(
    dt_terminal[, .SD, .SDcols = unique(c(Xvars, "node", "tij", "deltaij"))]
  )

  train_d[, gid := .GRP, by = Xvars]

  group_key <- unique(train_d[, c(Xvars, "gid"), with = FALSE])

  train_d <- train_d[!is.na(gid) & !is.na(node) & !is.na(tij)]
  train_d <- train_d[, .(
    tij = sum(tij),
    deltaij = sum(as.numeric(deltaij == 1L)),
    N_terminal = .N
  ), by = .(gid, node)]

  train_d[, node_ix := match(node, nodes)]
  data.table::setorder(train_d, gid, node_ix)

  expanded <- poissonsuperlearner:::expand_terminal_grouped_cpp(
    gid = train_d$gid,
    node_ix = train_d$node_ix,
    tij = train_d$tij,
    deltaij = train_d$deltaij,
    N = train_d$N_terminal,
    widths = c(diff(nodes), 1.0)
  )

  dt_fit <- data.table::as.data.table(expanded)
  dt_fit[, node := nodes[node_ix]]
  dt_fit <- group_key[dt_fit, on = "gid"]

  # 1) time preservation after learner expansion
  d_with_gid <- group_key[d, on = Xvars]
  time_check <- dt_fit[, .(sum_tij = sum(tij)), by = gid][
    d_with_gid[, .(obs_time = sum(get(time_col))), by = gid],
    on = "gid"
  ]

  testthat::expect_true(all(is.finite(time_check$sum_tij)))
  testthat::expect_true(all(is.finite(time_check$obs_time)))
  testthat::expect_true(all(abs(time_check$sum_tij - time_check$obs_time) < 1e-8))

  # 2) covariate preservation after learner expansion
  for (v in Xvars) {
    if (!v %in% names(d)) next

    if (!v %in% names(dt_fit)) {
      testthat::fail(paste0("Covariate '", v, "' missing from preprocessed data"))
    }

    varying <- dt_fit[, data.table::uniqueN(get(v)), by = gid][V1 != 1]
    testthat::expect_equal(nrow(varying), 0)

    merged <- unique(dt_fit[, .(gid, val_fit = get(v))])[
      group_key[, .(gid, val_orig = get(v))],
      on = "gid"
    ]

    if (is.factor(merged$val_fit) || is.character(merged$val_fit)) {
      testthat::expect_true(all(as.character(merged$val_fit) == as.character(merged$val_orig)))
    } else {
      testthat::expect_true(all(abs(merged$val_fit - merged$val_orig) < 1e-10))
    }
  }

  # 3) structure
  testthat::expect_true(all(c("gid", "node", "node_ix", "tij", "deltaij", "N_terminal") %in% names(dt_fit)))
  observed_node_ix <- sort(unique(dt_fit$node_ix))
  testthat::expect_equal(observed_node_ix, seq_len(max(observed_node_ix)))
  testthat::expect_true(max(observed_node_ix) <= length(nodes))
  testthat::expect_true(all(dt_fit$tij >= 0))
})


testthat::test_that("Poisson piecewise-constant fit reproduces Cox coefficient (X1) - both Superlearner and fit_learner", {

  testthat::skip_if_not_installed("riskRegression")
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("glmnet")
  testthat::skip_if_not_installed("data.table")
  {
    set.seed(42)
    d <- riskRegression::sampleData(n = 50, formula = ~ f(X1, 2))
  }
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

  learner <- poissonsuperlearner::Learner_gam(
    covariates = "X1",
    cross_validation = T
  )


  olcheck <- poissonsuperlearner::Superlearner(
    data = d,
    learners = list(list(learner)),
    id = "id",
    status = "status",
    event_time = "time")

  fit_ps <- poissonsuperlearner::fit_learner(
    data = d,
    learner = learner,
    id = "id",
    status = "status",
    event_time = "time"
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
