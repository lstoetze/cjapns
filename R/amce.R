#' Internal AMCE estimation via difference-in-means
#'
#' Computes Average Marginal Component Effects using the bivariate
#' subclassification estimator (Proposition 5 in Stoetzer & Magazinnik, 2026).
#' For each attribute, AMCE = mean(Y | level = tq) - mean(Y | level = tp).
#'
#' @param formula A formula: `outcome ~ attribute1 + attribute2 + ...`.
#' @param data A data.frame in long format (one row per profile).
#' @param id A one-sided formula for the respondent ID (e.g., `~ ResponseId`).
#'   Used for cluster-robust standard errors.
#' @param task_var Character name of the task-number variable (e.g., `"time"`).
#'   Required when `informative = "informative"`.
#' @param informative Whether to restrict to informative tasks (`"informative"`)
#'   or use all tasks (`"all"`, default). See `.filter_informative`.
#' @return A list with elements `amce` (per-attribute estimates),
#'   `coefficients`, `se`, `attributes`.
#' @keywords internal
estimate_amce <- function(formula, data, id = NULL, task_var = NULL,
                          informative = c("all", "informative")) {
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
  do_filter <- informative == "informative" && !is.null(task_var) && !is.null(id_var)
  if (do_filter) task_key <- paste(data[[id_var]], data[[task_var]], sep = ":::")

  amce_list <- list()
  all_beta <- c()
  all_se <- c()

  for (a in attr_names) {
    levs <- attributes_info[[a]]
    base_level <- levs[1]
    a_coefs <- c()
    a_se <- c()

    for (lev in levs[-1]) {
      if (do_filter) {
        keep   <- .filter_informative(data[[a]], task_key, lev, base_level)
        d_pair <- data[keep, , drop = FALSE]
      } else {
        d_pair <- data
      }
      Y_pair <- d_pair[[outcome_var]]
      idx_tq <- which(d_pair[[a]] == lev)
      idx_tp <- which(d_pair[[a]] == base_level)

      if (length(idx_tq) == 0 || length(idx_tp) == 0) {
        a_coefs[lev] <- NA_real_
        a_se[lev]    <- NA_real_
        all_beta[paste0(a, lev)] <- NA_real_
        all_se[paste0(a, lev)]   <- NA_real_
        next
      }

      amce_val <- mean(Y_pair[idx_tq], na.rm = TRUE) - mean(Y_pair[idx_tp], na.rm = TRUE)

      if (!is.null(id_var)) {
        se_val <- .cluster_se_dim(Y_pair, d_pair[[a]], lev, base_level, d_pair[[id_var]])
      } else {
        se_val <- sqrt(var(Y_pair[idx_tq], na.rm = TRUE) / length(idx_tq) +
                       var(Y_pair[idx_tp], na.rm = TRUE) / length(idx_tp))
      }

      a_coefs[lev] <- amce_val
      a_se[lev]    <- se_val
      all_beta[paste0(a, lev)] <- amce_val
      all_se[paste0(a, lev)]   <- se_val
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
#' @param attr_vals Attribute column from the data (vector).
#' @param task_key Character vector identifying each task (respondent x task).
#' @param tq,tp The two levels being compared.
#' @return Logical vector of length `length(attr_vals)`, TRUE for rows in
#'   informative tasks.
#' @keywords internal
.filter_informative <- function(attr_vals, task_key, tq, tp) {
  a_char   <- as.character(attr_vals)
  tq_c     <- as.character(tq)
  tp_c     <- as.character(tp)
  task_n   <- tapply(rep(1L, length(a_char)), task_key, sum)
  task_ntq <- tapply(a_char == tq_c, task_key, sum)
  task_ntp <- tapply(a_char == tp_c, task_key, sum)
  is_inf   <- (task_ntq == 1L & task_ntp == (task_n - 1L)) |
              (task_ntp == 1L & task_ntq == (task_n - 1L))
  task_key %in% names(is_inf)[is_inf]
}


#' Cluster-robust SE for a difference-in-means
#' @keywords internal
.cluster_se_dim <- function(Y, treatment, level_tq, level_tp, cluster) {
  keep <- treatment %in% c(level_tq, level_tp)
  Y <- Y[keep]; D <- as.integer(treatment[keep] == level_tq)
  cluster <- cluster[keep]
  fit <- stats::lm(Y ~ D)
  X <- stats::model.matrix(fit)
  n <- nrow(X); p <- ncol(X); e <- stats::residuals(fit)
  clusters <- unique(cluster); M <- length(clusters)
  bread <- solve(crossprod(X))
  meat <- matrix(0, p, p)
  for (g in clusters) {
    idx <- which(cluster == g)
    score_g <- crossprod(X[idx, , drop = FALSE], e[idx])
    meat <- meat + tcrossprod(score_g)
  }
  correction <- (M / (M - 1)) * ((n - 1) / (n - p))
  V <- bread %*% (correction * meat) %*% bread
  sqrt(V[2, 2])
}


#' Get AMCE for an arbitrary pair of levels
#'
#' Reconstructs AMCE(tq, tp) from regression coefficients relative to the
#' base level.
#'
#' @param amce_result Output of \code{estimate_amce}.
#' @param attribute Attribute name (character).
#' @param tq,tp The two levels to compare.
#' @param base The base (reference) level.
#' @return Scalar AMCE estimate.
#' @keywords internal
get_amce_for_pair <- function(amce_result, attribute, tq, tp, base) {
  amce_a <- amce_result$amce[[attribute]]
  if (tq == base && tp == base) return(0)
  if (tq == base) return(-amce_a$estimate[tp])
  if (tp == base) return(amce_a$estimate[tq])
  amce_a$estimate[tq] - amce_a$estimate[tp]
}
