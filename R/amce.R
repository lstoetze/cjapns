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
#' @return A list with elements `amce` (per-attribute estimates),
#'   `coefficients`, `se`, `attributes`.
#' @keywords internal
estimate_amce <- function(formula, data, id = NULL) {

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

  Y <- data[[outcome_var]]
  amce_list <- list()
  all_beta <- c()
  all_se <- c()

  for (a in attr_names) {
    levs <- attributes_info[[a]]
    base_level <- levs[1]
    a_coefs <- c()
    a_se <- c()

    for (lev in levs[-1]) {
      idx_tq <- which(data[[a]] == lev)
      idx_tp <- which(data[[a]] == base_level)
      amce_val <- mean(Y[idx_tq], na.rm = TRUE) - mean(Y[idx_tp], na.rm = TRUE)

      if (!is.null(id)) {
        id_var <- all.vars(id)
        se_val <- .cluster_se_dim(Y, data[[a]], lev, base_level, data[[id_var]])
      } else {
        se_val <- sqrt(var(Y[idx_tq], na.rm = TRUE) / length(idx_tq) +
                       var(Y[idx_tp], na.rm = TRUE) / length(idx_tp))
      }

      a_coefs[lev] <- amce_val
      a_se[lev] <- se_val
      all_beta[paste0(a, lev)] <- amce_val
      all_se[paste0(a, lev)] <- se_val
    }

    amce_list[[a]] <- list(
      base_level = base_level, levels = levs,
      estimate = a_coefs, se = a_se
    )
  }

  list(coefficients = all_beta, se = all_se, formula = formula,
       attributes = attributes_info, amce = amce_list)
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
