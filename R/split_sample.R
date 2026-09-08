#' Point difference-in-means for a level pair, without standard errors
#'
#' Same point-estimate logic as \code{\link{estimate_pairwise}} (informative-
#' task filtering, optional weighting), but skips the cluster-robust SE
#' computation entirely. Used by \code{\link{.split_sample_mapns}} for the
#' selection-half computation, which calls this thousands of times (every
#' candidate pair, on every split) and only ever needs the point estimate
#' there to pick the argmax; the estimation-half computation needs an SE
#' and uses \code{\link{estimate_pairwise}} directly instead.
#' @keywords internal
.point_diff_pairwise <- function(data, outcome_var, attribute, tq, tp,
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

  if (length(idx_tq) == 0 || length(idx_tp) == 0) return(NA_real_)

  stats::weighted.mean(Y_pair[idx_tq], w_pair[idx_tq]) -
    stats::weighted.mean(Y_pair[idx_tp], w_pair[idx_tp])
}

#' Split-sample point estimate (+ SE) of the separability MAPNS, one attribute
#'
#' Internal workhorse shared by \code{\link{mapns_split_test}} (diagnostic:
#' reports the full distribution of split-based estimates against the
#' plug-in max) and \code{cj_apns(split_sample = TRUE)} (uses the median
#' split-based estimate in place of the plug-in max as the reported point
#' estimate, with an SE/CI via the combining rule described below). See
#' \code{\link{mapns_split_test}} for the method and its limitations -- in
#' particular, this is a bias-reduction device, not an unbiased estimator of
#' the true max.
#'
#' The reported SE combines two variance components across the \code{m}
#' valid splits, following Rubin's rules for combining several randomly
#' varying analyses of the same data (as in multiple imputation): the
#' within-split variance \eqn{\bar W} (the mean of each split's own
#' cluster-robust sampling variance for its selected pair, estimated on the
#' held-out estimation half) and the between-split variance \eqn{B} (the
#' variance, across splits, of the split-level point estimates themselves).
#' The combined variance is \eqn{\bar W + (1 + 1/m) B}. The between-split
#' term is what inflates the interval when pair-selection is unstable across
#' splits (near-tied candidate pairs) -- exactly the winner's-curse-relevant
#' case -- and shrinks to near zero when the same pair reliably wins every
#' split. Rubin's rules are derived for combining *means*; applying the
#' resulting SE alongside the reported *median* (chosen for robustness
#' against the occasional degenerate split -- see `cj_apns`'s `split_sample`
#' documentation) is a practical approximation, not an exact match, and can
#' understate/overstate the true sampling variability of the median when
#' `split_vals` is heavily skewed.
#'
#' @param data A data.frame already restricted/prepared for this attribute
#'   (e.g. with a `.design_weight` column already attached if `weights_col`
#'   is supplied).
#' @param levs Character vector of this attribute's levels.
#' @param weights_col Character name of a weight column in `data`, or `NULL`.
#' @param n_splits Number of independent random respondent-level splits.
#' @param prop Share of respondents assigned to the selection half. Default
#'   0.3, i.e. a 30/70 selection/estimation split -- shifting more data to
#'   the estimation half tightens the confidence interval on the selected
#'   pair somewhat, at the cost of a noisier (but still consistent)
#'   selection step.
#' @param alpha Significance level for the combined-SE confidence interval.
#' @return A list with `estimate` (median split-based |contrast| at the
#'   argmax pair, across valid splits), `se`, `lower`, `upper` (`NA` if fewer
#'   than two splits yield a valid within-split SE), `n_valid_splits`
#'   (splits with both a valid point estimate and a valid, finite
#'   cluster-robust SE), `n_splits`, and the raw per-split vectors
#'   `split_vals`, `split_ses`, `selected_pair` (used by
#'   \code{\link{mapns_split_test}} to additionally report `split_mean`,
#'   selection entropy, and per-split detail).
#' @keywords internal
.split_sample_mapns <- function(data, outcome_var, attribute, levs,
                                 id_var, task_var, informative, profile_var,
                                 weights_col = NULL, n_splits = 200, prop = 0.3,
                                 alpha = 0.05) {
  pairs <- list()
  for (q in seq_along(levs)) for (p in seq_along(levs)) {
    if (q >= p) next
    pairs[[length(pairs) + 1]] <- c(levs[q], levs[p])
  }
  pair_labels <- vapply(pairs, function(pr) paste0(pr[1], " vs ", pr[2]), character(1))

  unique_ids <- unique(data[[id_var]])
  n_a        <- round(prop * length(unique_ids))

  split_vals    <- rep(NA_real_, n_splits)
  split_ses     <- rep(NA_real_, n_splits)
  selected_pair <- rep(NA_character_, n_splits)

  for (s in seq_len(n_splits)) {
    ids_a <- sample(unique_ids, size = n_a)
    in_a  <- data[[id_var]] %in% ids_a
    dA <- data[in_a, , drop = FALSE]
    dB <- data[!in_a, , drop = FALSE]

    sel_vals <- vapply(pairs, function(pr) {
      abs(.point_diff_pairwise(dA, outcome_var, attribute, pr[1], pr[2],
                                id_var = id_var, task_var = task_var,
                                informative = informative, profile_var = profile_var,
                                weights_col = weights_col))
    }, numeric(1))
    if (all(is.na(sel_vals))) next
    winner <- which.max(sel_vals)
    pr     <- pairs[[winner]]
    selected_pair[s] <- pair_labels[winner]

    # Full estimate_pairwise() (not the lean .point_diff_pairwise()) for this
    # one estimation-half call only, since we need its cluster-robust SE for
    # the within-split variance component below; the selection step above
    # never needs an SE, only the magnitude used to pick the argmax.
    est_b <- estimate_pairwise(dB, outcome_var, attribute, pr[1], pr[2],
                               id_var = id_var, task_var = task_var,
                               informative = informative,
                               profile_var = profile_var,
                               weights_col = weights_col)
    split_vals[s] <- abs(est_b$estimate)
    split_ses[s]  <- est_b$se
  }

  estimate <- stats::median(split_vals, na.rm = TRUE)
  valid <- !is.na(split_vals) & !is.na(split_ses) & is.finite(split_ses)
  m <- sum(valid)

  se <- NA_real_; lower <- NA_real_; upper <- NA_real_
  if (m >= 2) {
    w_bar     <- mean(split_ses[valid]^2)
    b_var     <- stats::var(split_vals[valid])
    total_var <- w_bar + (1 + 1 / m) * b_var
    se        <- sqrt(total_var)
    z         <- stats::qnorm(1 - alpha / 2)
    lower     <- max(0, estimate - z * se)
    upper     <- estimate + z * se
  }

  list(
    estimate       = estimate,
    se             = se,
    lower          = lower,
    upper          = upper,
    n_valid_splits = m,
    n_splits       = n_splits,
    split_vals     = split_vals,
    split_ses      = split_ses,
    selected_pair  = selected_pair
  )
}

#' Split-sample diagnostic for max-selection ("winner's curse") bias in MAPNS
#'
#' Under separable monotonicity, \eqn{\mathrm{MAPNS}_l} is the maximum of
#' several noisy pairwise APNS estimates. Because \code{max(.)} is convex,
#' the plug-in sample max is upward-biased for the true max (Jensen's
#' inequality): the same estimate that "wins" the max is also the one used
#' to report its value.
#'
#' This function repeatedly splits respondents into two independent halves:
#' one half (the selection sample) picks which pair attains the max; the
#' other half (the estimation sample) estimates that already-fixed pair's
#' APNS. Because selection and estimation use disjoint respondents, the
#' reported value is no longer conditioned on its own selection, which
#' removes the winner's-curse channel specifically.
#'
#' This is a diagnostic, not a replacement point estimator: it does not
#' recover an unbiased estimate of the true max. When the selection sample
#' picks the wrong pair (more likely when the selection sample is noisy, or
#' when the best and second-best pairs are a near-tie), the split estimate
#' instead targets that generically smaller pair's true APNS, biasing the
#' split-based figure downward rather than upward. Comparing the ordinary
#' plug-in max to the distribution of split-based estimates across many
#' random splits indicates whether winner's-curse bias is large enough to
#' matter in a given application; it is not a corrected estimate to report
#' in place of the plug-in MAPNS.
#'
#' @param formula A formula: `outcome ~ attribute1 + attribute2 + ...`.
#' @param data A data.frame in long format (one row per profile).
#' @param id A one-sided formula for the respondent ID (e.g., `~ ResponseId`).
#'   Splitting is done at the respondent level, since informative tasks from
#'   the same respondent are correlated and the two halves must be
#'   independent.
#' @param tasks A one-sided formula for the task-number variable, or `NULL`.
#' @param profile A one-sided formula for the profile indicator, or `NULL`.
#' @param informative Whether to restrict to informative tasks
#'   (`"informative"`, default) or use all tasks (`"all"`).
#' @param design Either `"uniform"` (default) or a `"cj_design"` object from
#'   \code{\link{make_design}}; see \code{\link{cj_apns}}.
#' @param n_splits Number of independent random respondent-level splits.
#'   Default 200.
#' @param prop Share of respondents assigned to the selection half. Default
#'   0.3, i.e. a 30/70 selection/estimation split.
#' @param alpha Significance level for the `lower`/`upper` confidence
#'   interval on the split-based estimate. Default 0.05.
#' @param seed Optional integer seed for reproducibility.
#' @return A data.frame with one row per attribute in `formula`:
#'   \describe{
#'     \item{`item`}{Attribute name.}
#'     \item{`plugin_mapns`}{The ordinary full-sample plug-in max estimate
#'       (Theorem "Nonparametric Estimation of the MAPNS").}
#'     \item{`split_median`, `split_mean`}{Median/mean of the split-based
#'       estimate (selection on one random half, estimation on the other)
#'       across valid splits. Comparable in magnitude to `plugin_mapns`; a
#'       systematic gap between them is the winner's-curse signature.}
#'     \item{`se`, `lower`, `upper`}{SE and confidence interval for
#'       `split_median`, via the same Rubin's-rules combination of within-
#'       and between-split variance used internally by
#'       \code{cj_apns(split_sample = TRUE)} -- see
#'       \code{\link{.split_sample_mapns}} for the formula and its
#'       approximations (in particular, the combining rule is derived for
#'       combining *means*; applying it to the reported *median* is a
#'       practical approximation). `NA` if fewer than two splits yield a
#'       usable within-split SE.}
#'     \item{`n_valid_splits`}{Number of splits (out of `n_splits`) with
#'       both a valid point estimate and a valid, finite within-split SE.
#'       Low values mean the check itself is underpowered (thin attribute),
#'       not that there is no bias.}
#'     \item{`selection_entropy`}{Shannon entropy, in bits, of how often
#'       each pair was selected as the argmax across splits. 0 means the
#'       same pair always wins (stable, low winner's-curse risk); higher
#'       values indicate a near-tie among several pairs, where selection
#'       noise -- and hence winner's-curse risk -- is largest.}
#'   }
#'   The per-split raw values (`split_estimate`, `split_se`, `selected_pair`
#'   for every split and attribute) are attached as `attr(x, "detail")` for
#'   plotting.
#' @export
mapns_split_test <- function(formula, data, id, tasks = NULL, profile = NULL,
                              informative = c("informative", "all"),
                              design = "uniform", n_splits = 200,
                              prop = 0.3, alpha = 0.05, seed = NULL) {
  informative <- match.arg(informative)
  if (!is.null(seed)) set.seed(seed)

  if (missing(id)) stop("'id' must be specified (e.g., id = ~ ResponseId).")
  id_var <- all.vars(id)
  stopifnot(length(id_var) == 1, id_var %in% names(data))

  if (informative == "informative" && is.null(tasks)) {
    warning("'informative = \"informative\"' requires 'tasks' to be specified. ",
            "Falling back to 'informative = \"all\"'.")
    informative <- "all"
  }
  task_var <- if (!is.null(tasks)) all.vars(tasks)[1] else NULL
  if (!is.null(task_var))
    stopifnot("'tasks' variable not found in data" = task_var %in% names(data))

  profile_var <- if (!is.null(profile)) all.vars(profile)[1] else NULL
  if (!is.null(profile_var))
    stopifnot("'profile' variable not found in data" = profile_var %in% names(data))

  tt <- stats::terms(formula, data = data)
  attr_names <- attr(tt, "term.labels")
  attr_names <- attr_names[!grepl(":", attr_names)]
  outcome_var <- all.vars(formula)[1]
  for (a in attr_names)
    if (!is.factor(data[[a]])) data[[a]] <- as.factor(data[[a]])

  detail <- list()

  rows <- lapply(attr_names, function(a) {
    levs <- levels(data[[a]])
    pairs <- list()
    for (q in seq_along(levs)) for (p in seq_along(levs)) {
      if (q >= p) next
      pairs[[length(pairs) + 1]] <- c(levs[q], levs[p])
    }

    d_a <- data
    weights_col <- NULL
    if (!identical(design, "uniform")) {
      d_a$.design_weight <- .attach_design_weights(d_a, attr_names, a, design,
                                                     id_var, task_var)
      weights_col <- ".design_weight"
    }

    full_vals <- vapply(pairs, function(pr) {
      abs(.point_diff_pairwise(d_a, outcome_var, a, pr[1], pr[2],
                                id_var = id_var, task_var = task_var,
                                informative = informative, profile_var = profile_var,
                                weights_col = weights_col))
    }, numeric(1))
    plugin_mapns <- max(full_vals, na.rm = TRUE)

    # Shared workhorse with cj_apns(split_sample = TRUE): same splits, same
    # point estimate (median) + SE/CI (Rubin's-rules combination) logic, so
    # the two can never drift apart again.
    ss <- .split_sample_mapns(d_a, outcome_var, a, levs, id_var, task_var,
                               informative, profile_var, weights_col,
                               n_splits = n_splits, prop = prop, alpha = alpha)

    sel_tab <- table(ss$selected_pair[!is.na(ss$selected_pair)])
    sel_p   <- sel_tab / sum(sel_tab)
    entropy <- if (length(sel_p) > 0) -sum(sel_p * log2(sel_p)) else NA_real_

    detail[[a]] <<- data.frame(item = a, split_estimate = ss$split_vals,
                                split_se = ss$split_ses,
                                selected_pair = ss$selected_pair,
                                stringsAsFactors = FALSE)

    data.frame(
      item               = a,
      plugin_mapns       = plugin_mapns,
      split_median       = ss$estimate,
      split_mean         = mean(ss$split_vals, na.rm = TRUE),
      se                 = ss$se,
      lower              = ss$lower,
      upper              = ss$upper,
      n_valid_splits     = ss$n_valid_splits,
      n_splits           = ss$n_splits,
      selection_entropy  = entropy,
      stringsAsFactors   = FALSE
    )
  })

  out <- do.call(rbind, rows)
  attr(out, "detail") <- do.call(rbind, detail)
  out
}
