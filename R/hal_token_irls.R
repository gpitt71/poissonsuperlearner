# Token-only HAL construction. Matches the primitive knot definition and
# column enumeration of Learner_hal for arbitrary interaction degrees.
hal3_tokens <- function(data, vars, max_degree = 2L, num_knots = c(10L, 5L)) {
  started <- as.numeric(Sys.time())
  stage <- numeric()
  mark <- function(name) {
    now <- as.numeric(Sys.time())
    stage[[name]] <<- now - started - sum(stage)
    invisible(NULL)
  }
  stopifnot(length(max_degree) == 1L, !is.na(max_degree), max_degree >= 1L,
            length(num_knots) == max_degree,
            all(vars %in% names(data)), !anyDuplicated(vars))
  n <- nrow(data); d <- length(vars)
  xnum <- matrix(NA_real_, n, d)
  xcat <- matrix(NA_integer_, n, d)
  levels_used <- setNames(vector('list', d), vars)
  primitives <- setNames(vector('list', d), vars)
  for (v in seq_along(vars)) {
    z <- data[[vars[v]]]
    if (is.character(z) || is.logical(z)) z <- factor(z)
    if (is.numeric(z)) {
      xnum[,v] <- z
      zz <- z[is.finite(z)]
      if (!length(zz)) { cuts <- numeric() } else {
        # Equal to mk_main_numeric_cpp: type-1 quantiles, drop the two
        # endpoints, unique, and exclude the minimum/all-ones column.
        sorted <- sort(zz)
        K <- as.integer(num_knots[1])
        ranks <- pmax(1L, ceiling((0:(K+1)) / (K+1) * length(sorted)))
        cuts <- sort(unique(sorted[ranks[-c(1, length(ranks))]]))
        cuts <- cuts[cuts > sorted[1]]
      }
      primitives[[v]] <- lapply(cuts, function(cut)
        list(var = vars[v], kind = 0L, cut = cut, level = -1L))
    } else if (is.factor(z)) {
      levels_used[[v]] <- levels(z)
      xcat[,v] <- as.integer(z)
      primitives[[v]] <- lapply(seq_along(levels(z))[-1L], function(lev)
        list(var = vars[v], kind = 1L, cut = 0, level = as.integer(lev)))
    } else stop('Unsupported variable type: ',vars[v])
  }
  mark("data_and_primitives")
  cols <- unlist(primitives, recursive = FALSE)
  if (max_degree >= 2L && d >= 2L) {
    for (degree in 2:min(max_degree, d)) {
      for (combination in utils::combn(seq_along(vars), degree,
                                       simplify = FALSE)) {
        parts <- lapply(combination, function(v) {
          z <- primitives[[v]]
          if (length(z) && z[[1L]]$kind == 0L)
            z <- head(z, num_knots[degree])
          z
        })
        if (any(lengths(parts) == 0L)) next
        # expand.grid varies its first argument fastest, exactly like the
        # mixed-radix enumeration in inter_cols_batch_cpp().
        selection <- do.call(expand.grid, c(
          lapply(parts, seq_along), KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE))
        for (row in seq_len(nrow(selection))) {
          cols[[length(cols) + 1L]] <- lapply(seq_len(degree), function(k)
            parts[[k]][[selection[row, k]]])
        }
      }
    }
  }
  mark("column_enumeration")
  # Normalize main columns to lists of parts.
  cols <- lapply(cols, function(z) if (!is.null(z$var)) list(z) else z)
  encode <- function(part) c(var = match(part$var, vars)-1L,
                              kind = part$kind, cut = part$cut, level = part$level)
  p <- length(cols)
  if (!p) stop('Empty HAL token basis.')
  degree <- lengths(cols)
  offsets <- c(0L, cumsum(degree))
  encoded <- do.call(rbind, lapply(cols, function(z)
    do.call(rbind, lapply(z, encode))))
  out <- list(numeric_data=xnum, cat_data=xcat,
    token_offsets=as.integer(offsets), token_var=as.integer(encoded[,1L]),
    token_kind=as.integer(encoded[,2L]), token_cut=as.numeric(encoded[,3L]),
    token_level=as.integer(encoded[,4L]),
    vars=vars, levels_used=levels_used, columns=cols)
  mark("token_encoding")
  # Eliminate empty columns in the same manner as hal_basis() without
  # constructing a dgCMatrix. The sums use the token structure itself.
  index <- hal_token_index_cpp(out$numeric_data, out$cat_data,
    out$token_offsets, out$token_var, out$token_kind, out$token_cut,
    out$token_level)
  sums <- hal_token_sums_cached_cpp(index,out$numeric_data,out$cat_data,rep(1,n),rep(0,n))
  keep <- which(sums$curvature > 0)
  mark("index_and_empty_scan")
  if (!length(keep)) stop('All HAL token columns are empty.')
  if (length(keep) != length(out$columns)) {
    out$columns <- out$columns[keep]
    degree <- lengths(out$columns)
    out$token_offsets <- as.integer(c(0L, cumsum(degree)))
    encoded <- do.call(rbind, lapply(out$columns, function(z)
      do.call(rbind, lapply(z, encode))))
    out$token_var <- as.integer(encoded[,1L])
    out$token_kind <- as.integer(encoded[,2L])
    out$token_cut <- as.numeric(encoded[,3L])
    out$token_level <- as.integer(encoded[,4L])
    out$index <- hal_token_index_cpp(out$numeric_data, out$cat_data,
      out$token_offsets, out$token_var, out$token_kind, out$token_cut,
      out$token_level)
  } else {
    # The first index is already final. Avoid repeating every sort and index
    # allocation when the candidate enumeration contains no empty supports.
    out$index <- index
  }
  mark("final_index_rebuild")
  out$native_index_bytes <- hal_token_index_bytes_cpp(out$index)
  out$colnames <- unname(vapply(out$columns, function(z) paste(vapply(z, function(part) {
    if (part$kind == 0L) sprintf('I(%s>=%.10g)',part$var,part$cut)
    else sprintf('I(%s==%s)',part$var,make.names(out$levels_used[[part$var]][part$level]))
  }, character(1)), collapse=':'), character(1)))
  out$default_penalty_factor <- vapply(out$columns, function(parts)
    as.numeric(!(length(parts) == 1L &&
      identical(parts[[1L]]$var, "node") && parts[[1L]]$kind == 1L)),
    numeric(1L))
  out$ncol <- length(out$columns)
  # The nested lists are convenient while enumerating and naming columns but
  # duplicate the compact flattened representation used by every kernel.
  out$columns <- NULL
  mark("column_names")
  out$preprocessing_seconds <- stage
  out
}

hal3_fit_tokens <- function(tokens,y,offset,lambda,penalty_factor=NULL,intercept=TRUE,
                            tol=1e-8,maxit_irls=100L,maxit_cd=10000L,active_set=TRUE,
                            strong_rules=TRUE,rescale_penalty_factor=TRUE,
                            standardize=FALSE) {
  p <- tokens$ncol
  if (is.null(penalty_factor)) penalty_factor <- rep(1,p)
  if (length(standardize) != 1L || is.na(standardize))
    stop("`standardize` must be TRUE or FALSE.")
  penalty_factor <- as.numeric(penalty_factor)
  if (length(penalty_factor) != p || any(!is.finite(penalty_factor)) ||
      any(penalty_factor < 0))
    stop("Invalid `penalty_factor`.")

  # glmnet first rescales the user-supplied penalty factors to sum to the
  # number of columns, and then fits on standardized columns.  Expressed in
  # the original (0/1 indicator) coordinates, this multiplies each L1 penalty
  # by that column's population standard deviation.  Applying that equivalent
  # penalty lets the token solver retain its sparse indicator representation.
  penalty_scale <- if (isTRUE(rescale_penalty_factor)) {
    if (!any(penalty_factor > 0)) stop("All penalties are zero.")
    p / sum(penalty_factor)
  } else 1
  column_scale <- rep(1, p)
  penalty_use <- penalty_factor
  rescale_use <- rescale_penalty_factor
  if (isTRUE(standardize)) {
    prevalence <- hal_token_sums_cached_cpp(
      tokens$index, tokens$numeric_data, tokens$cat_data,
      rep(1, nrow(tokens$numeric_data)), rep(0, nrow(tokens$numeric_data)))
    prevalence <- as.numeric(prevalence$curvature)
    column_scale <- if (isTRUE(intercept)) {
      sqrt(pmax(0, prevalence - prevalence^2))
    } else {
      sqrt(prevalence)
    }
    penalty_use <- penalty_factor * penalty_scale * column_scale
    rescale_use <- FALSE
  }
  z <- hal_token_irls_cpp(tokens$index,tokens$numeric_data,tokens$cat_data,
    as.numeric(y),as.numeric(offset),as.numeric(lambda),penalty_use,
    isTRUE(intercept),as.numeric(tol),as.integer(maxit_irls),as.integer(maxit_cd),
    isTRUE(active_set),isTRUE(rescale_use),isTRUE(strong_rules))
  # Report the scaling of the caller's factors, rather than the internal
  # already-combined factors passed to the native solver.
  z$penalty_factor_scale <- penalty_scale
  z$standardize <- isTRUE(standardize)
  z$column_scale <- column_scale
  z$effective_penalty_factor <- penalty_factor * penalty_scale * column_scale
  rownames(z$beta) <- tokens$colnames
  z$tokens <- tokens
  z$selected_lambda <- tail(z$lambda,1L)
  class(z) <- 'hal_token_irls_fit'
  z
}

#' @export
coef.hal_token_irls_fit <- function(object,s=NULL,...) {
  j <- if(is.null(s)) length(object$lambda) else which.min(abs(object$lambda-s))
  rbind('(Intercept)'=object$a0[j],object$beta[,j,drop=FALSE])
}

hal3_predict_tokens <- function(object,newdata,s=NULL,type=c('response','link'),newoffset=NULL) {
  type <- match.arg(type)
  t <- object$tokens
  if (!all(t$vars %in% names(newdata))) stop('Missing covariates')
  j <- if(is.null(s)) length(object$lambda) else which.min(abs(object$lambda-s))
  eta <- rep(object$a0[j],nrow(newdata))
  for (k in which(object$beta[,j] != 0)) {
    good <- hal3_token_indicator(t, newdata, k)
    eta[which(good)] <- eta[which(good)] + object$beta[k,j]
  }
  if(!is.null(newoffset)) eta <- eta+newoffset
  if(type=='response') exp(eta) else eta
}

#' @export
predict.hal_token_irls_fit <- function(object,newdata,s=NULL,type=c('response','link'),newoffset=NULL,...)
  hal3_predict_tokens(object,newdata,s,type,newoffset)

hal3_token_indicator <- function(tokens, data, column) {
  if (length(column) != 1L || is.na(column) || column < 1L ||
      column > tokens$ncol) stop("Invalid token column.")
  begin <- tokens$token_offsets[column] + 1L
  end <- tokens$token_offsets[column + 1L]
  good <- rep(TRUE, nrow(data))
  for (at in begin:end) {
    variable <- tokens$vars[tokens$token_var[at] + 1L]
    if (tokens$token_kind[at] == 0L) {
      good <- good & is.finite(data[[variable]]) &
        data[[variable]] >= tokens$token_cut[at]
    } else {
      good <- good & !is.na(data[[variable]]) &
        as.character(data[[variable]]) ==
          tokens$levels_used[[variable]][tokens$token_level[at]]
    }
  }
  good
}
