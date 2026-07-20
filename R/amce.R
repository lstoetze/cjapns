#' Direct pairwise difference-in-means estimator (profile 1/2 conditioned)
#'
#' Estimates the pairwise contrast
#' \eqn{E[Y_{i1k} \mid T_{i1kl}=t_q, T_{i2kl}=t_p] - E[Y_{i1k} \mid T_{i1kl}=t_p, T_{i2kl}=t_q]}
#' directly on tasks informative for \eqn{(t_q,t_p)} (Definition 3), by
#' conditioning jointly on both profiles' levels for attribute \code{l}. This
#' is the plug-in estimator for Theorem 1 (APNS/AMCE) and, on a subclassed
#' \code{data}, for Theorem 2/4 (CAPNS/MAPNS) — evaluated directly at
#' whichever pair is needed, with no detour through a base-level regression
#' and no reconstruction of non-base contrasts from base-relative
#' coefficients (which is only valid under an unstated additivity
#' assumption and can diverge materially from the direct estimate whenever
#' neither \code{tq} nor \code{tp} is the reference level).
#'
#' @param data A data.frame in long format (one row per profile).
#' @param outcome_var Character name of the outcome column.
#' @param attribute Character name of the attribute column.
#' @param tq,tp The two levels to compare.
#' @param id_var Character name of the respondent ID column, or `NULL`.
#'   Used for cluster-robust standard errors.
#' @param task_var Character name of the task-number variable, or `NULL`.
#'   Required when `informative = "informative"`.
#' @param informative Whether to restrict to informative tasks
#'   (`"informative"`) or use all tasks (`"all"`, default).
#' @param profile_var Character name of the profile indicator variable, or
#'   `NULL`. See \code{\link{.filter_informative}}.
#' @param weights_col Character name of a numeric weight column in `data`,
#'   or `NULL` (default) for the unweighted pooled estimator (Corollary
#'   "Pooled Estimation under Full Randomization"). When supplied — see
#'   \code{\link{.attach_design_weights}} — implements the general,
#'   context-weighted estimator (Propositions "Nonparametric Estimation of
#'   the APNS/CAPNS"). Passing a column of all 1s reproduces the unweighted
#'   estimate exactly.
#' @return A list with `estimate` (signed difference in means), `se`
#'   (cluster-robust if `id_var` is given, else a simple two-sample SE),
#'   and `n_tq`, `n_tp` (informative-task profile counts for each level).
#' @keywords internal
estimate_pairwise <- function(data, outcome_var, attribute, tq, tp,
                              id_var = NULL, task_var = NULL,
                              informative = c("all", "informative"),
                              profile_var = NULL, weights_col = NULL) {
  informative <- match.arg(informative)
  do_filter <- informative == "informative" && !is.null(task_var) && !is.null(id_var)

  if (do_filter) {
    task_key     <- paste(data[[id_var]], data[[task_var]], sep = ":::")
    profile_vals <- if (!is.null(profile_var)) data[[profile_var]] else NULL
    keep         <- .filter_informative(data[[attribute]], task_key, tq, tp, profile_vals)
    d_pair       <- data[keep, , drop = FALSE]
  } else {
    d_pair <- data
  }

  Y_pair    <- d_pair[[outcome_var]]
  w_pair    <- if (!is.null(weights_col)) d_pair[[weights_col]] else rep(1, nrow(d_pair))
  attr_char <- as.character(d_pair[[attribute]])
  idx_tq    <- which(attr_char == as.character(tq))
  idx_tp    <- which(attr_char == as.character(tp))

  if (length(idx_tq) == 0 || length(idx_tp) == 0)
    return(list(estimate = NA_real_, se = NA_real_, n_tq = length(idx_tq), n_tp = length(idx_tp)))

  est <- stats::weighted.mean(Y_pair[idx_tq], w_pair[idx_tq]) -
         stats::weighted.mean(Y_pair[idx_tp], w_pair[idx_tp])

  is_weighted <- !is.null(weights_col) && any(w_pair != 1)

  se <- if (!is.null(id_var)) {
    .cluster_se_dim(Y_pair, d_pair[[attribute]], tq, tp, d_pair[[id_var]],
                    weights = if (is_weighted) w_pair else NULL)
  } else if (is_weighted) {
    .weighted_two_sample_se(Y_pair[idx_tq], w_pair[idx_tq], Y_pair[idx_tp], w_pair[idx_tp])
  } else {
    sqrt(stats::var(Y_pair[idx_tq], na.rm = TRUE) / length(idx_tq) +
         stats::var(Y_pair[idx_tp], na.rm = TRUE) / length(idx_tp))
  }

  list(estimate = est, se = se, n_tq = length(idx_tq), n_tp = length(idx_tp))
}

#' Horvitz-Thompson-style SE for a weighted two-sample mean difference
#' @keywords internal
.weighted_two_sample_se <- function(y1, w1, y2, w2) {
  .wvar <- function(y, w) {
    ybar <- stats::weighted.mean(y, w)
    sum(w^2 * (y - ybar)^2) / sum(w)^2
  }
  sqrt(.wvar(y1, w1) + .wvar(y2, w2))
}


#' Internal AMCE estimation via difference-in-means
#'
#' Computes Average Marginal Component Effects, relative to a base level, via
#' \code{\link{estimate_pairwise}} — one direct estimate per non-base level.
#' Used for `estimand = "amce"` reporting; pairwise APNS/ACMCE contrasts
#' between two arbitrary (possibly non-base) levels are computed by calling
#' \code{estimate_pairwise} directly rather than reconstructing them from
#' these base-relative coefficients.
#'
#' @param formula A formula: `outcome ~ attribute1 + attribute2 + ...`.
#' @param data A data.frame in long format (one row per profile).
#' @param id A one-sided formula for the respondent ID (e.g., `~ ResponseId`).
#'   Used for cluster-robust standard errors.
#' @param task_var Character name of the task-number variable (e.g., `"time"`).
#'   Required when `informative = "informative"`.
#' @param informative Whether to restrict to informative tasks (`"informative"`)
#'   or use all tasks (`"all"`, default). See `.filter_informative`.
#' @param profile_var Character name of the profile indicator variable
#'   (e.g., `"profile"`). When provided, informative task detection compares
#'   profiles explicitly rather than counting levels.
#' @param design Either `"uniform"` (default; unweighted pooled estimator)
#'   or a `"cj_design"` object from \code{\link{make_design}}, applying the
#'   general context-weighted estimator via \code{\link{.attach_design_weights}}.
#' @return A list with elements `amce` (per-attribute estimates),
#'   `coefficients`, `se`, `attributes`.
#' @keywords internal
estimate_amce <- function(formula, data, id = NULL, task_var = NULL,
                          informative = c("all", "informative"),
                          profile_var = NULL, design = "uniform") {
  informative <- match.arg(informative)

  tt <- stats::terms(formula, data = data)
  outcome_var <- all.vars(formula)[1]
  attr_names <- attr(tt, "term.labels")
  attr_names <- attr_names[!grepl(":", attr_names)]

  attributes_info <- lapply(attr_names, function(a) {
    if (!is.factor(data[[a]])) data[[a]] <- as.factor(data[[a]])
    levels(data[[a]])
  })
  names(attributes_info) <- attr_names

  for (a in attr_names) {
    if (!is.factor(data[[a]])) data[[a]] <- as.factor(data[[a]])
  }

  id_var <- if (!is.null(id)) all.vars(id) else NULL
  use_design <- !identical(design, "uniform")

  amce_list <- list()
  all_beta <- c()
  all_se <- c()

  for (a in attr_names) {
    levs <- attributes_info[[a]]
    base_level <- levs[1]
    a_coefs <- c()
    a_se <- c()

    weights_col <- NULL
    if (use_design) {
      data$.design_weight <- .attach_design_weights(data, attr_names, a, design,
                                                     id_var, task_var)
      weights_col <- ".design_weight"
    }

    for (lev in levs[-1]) {
      pw <- estimate_pairwise(data, outcome_var, a, lev, base_level,
                              id_var = id_var, task_var = task_var,
                              informative = informative, profile_var = profile_var,
                              weights_col = weights_col)
      a_coefs[lev] <- pw$estimate
      a_se[lev]    <- pw$se
      all_beta[paste0(a, lev)] <- pw$estimate
      all_se[paste0(a, lev)]   <- pw$se
    }

    amce_list[[a]] <- list(
      base_level = base_level, levels = levs,
      estimate = a_coefs, se = a_se
    )
  }

  list(coefficients = all_beta, se = all_se, formula = formula,
       attributes = attributes_info, amce = amce_list)
}


#' Filter rows to informative tasks for a level pair (Definition 3 in paper)
#'
#' A task is informative for pair (tq, tp) if exactly one profile shows tq and
#' all J-1 remaining profiles show tp, or vice versa.
#'
#' When `profile_vals` is provided the check compares the two profiles
#' explicitly (profile a vs profile b). Otherwise falls back to counting levels
#' within the task, which generalises to J > 2.
#'
#' @param attr_vals Attribute column from the data (vector).
#' @param task_key Character vector identifying each task (respondent x task).
#' @param tq,tp The two levels being compared.
#' @param profile_vals Optional profile indicator vector (e.g. "a"/"b").
#' @return Logical vector of length `length(attr_vals)`, TRUE for rows in
#'   informative tasks.
#' @keywords internal
.filter_informative <- function(attr_vals, task_key, tq, tp, profile_vals = NULL) {
  a_char <- as.character(attr_vals)
  tq_c   <- as.character(tq)
  tp_c   <- as.character(tp)

  if (!is.null(profile_vals)) {
    prof_char <- as.character(profile_vals)
    prof_levs <- sort(unique(prof_char))
    val_mat   <- tapply(a_char, list(task_key, prof_char), function(x) x[1])
    p1 <- val_mat[, prof_levs[1]]
    p2 <- val_mat[, prof_levs[2]]
    is_inf <- (!is.na(p1) & !is.na(p2)) &
              ((p1 == tq_c & p2 == tp_c) | (p1 == tp_c & p2 == tq_c))
    return(task_key %in% names(is_inf)[is_inf])
  }

  # Count-based fallback: works for any J
  task_n   <- tapply(rep(1L, length(a_char)), task_key, sum)
  task_ntq <- tapply(a_char == tq_c, task_key, sum)
  task_ntp <- tapply(a_char == tp_c, task_key, sum)
  is_inf   <- (task_ntq == 1L & task_ntp == (task_n - 1L)) |
              (task_ntp == 1L & task_ntq == (task_n - 1L))
  task_key %in% names(is_inf)[is_inf]
}


#' Cluster-robust SE for a (possibly weighted) difference-in-means
#'
#' When `weights` is `NULL`, this is an ordinary-least-squares cluster-robust
#' sandwich SE (unchanged from the unweighted estimator). When `weights` is
#' supplied, fits weighted least squares and uses the corresponding weighted
#' cluster-robust sandwich (bread \eqn{(X'WX)^{-1}}, meat built from
#' \eqn{w_i X_i e_i}), which reduces to the unweighted formula exactly when
#' all weights equal 1.
#' @keywords internal
.cluster_se_dim <- function(Y, treatment, level_tq, level_tp, cluster, weights = NULL) {
  keep <- treatment %in% c(level_tq, level_tp)
  Y <- Y[keep]; D <- as.integer(treatment[keep] == level_tq)
  cluster <- cluster[keep]
  w <- if (!is.null(weights)) weights[keep] else rep(1, length(Y))
  fit <- stats::lm(Y ~ D, weights = w)
  X <- stats::model.matrix(fit)
  n <- nrow(X); p <- ncol(X); e <- stats::residuals(fit)
  clusters <- unique(cluster); M <- length(clusters)
  bread <- solve(crossprod(X, w * X))
  meat <- matrix(0, p, p)
  for (g in clusters) {
    idx <- which(cluster == g)
    score_g <- crossprod(X[idx, , drop = FALSE], w[idx] * e[idx])
    meat <- meat + tcrossprod(score_g)
  }
  correction <- (M / (M - 1)) * ((n - 1) / (n - p))
  V <- bread %*% (correction * meat) %*% bread
  sqrt(V[2, 2])
}
