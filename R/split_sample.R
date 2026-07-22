#' Point difference-in-means for a level pair, without standard errors
#'
#' Same point-estimate logic as \code{\link{estimate_pairwise}} (informative-
#' task filtering, optional weighting), but skips the cluster-robust SE
#' computation entirely. Used by \code{\link{mapns_split_test}}, which calls
#' this thousands of times (every candidate pair, on every split) and only
#' ever needs the point estimate.
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
#'   0.5.
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
#'     \item{`n_valid_splits`}{Number of splits (out of `n_splits`) where
#'       the selected pair had informative-task coverage in the estimation
#'       half. Low values mean the check itself is underpowered (thin
#'       attribute), not that there is no bias.}
#'     \item{`selection_entropy`}{Shannon entropy, in bits, of how often
#'       each pair was selected as the argmax across splits. 0 means the
#'       same pair always wins (stable, low winner's-curse risk); higher
#'       values indicate a near-tie among several pairs, where selection
#'       noise -- and hence winner's-curse risk -- is largest.}
#'   }
#'   The per-split raw values (`split_estimate`, `selected_pair` for every
#'   split and attribute) are attached as `attr(x, "detail")` for plotting.
#' @export
mapns_split_test <- function(formula, data, id, tasks = NULL, profile = NULL,
                              informative = c("informative", "all"),
                              design = "uniform", n_splits = 200,
                              prop = 0.5, seed = NULL) {
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

  unique_ids <- unique(data[[id_var]])
  n_ids      <- length(unique_ids)
  n_a        <- round(prop * n_ids)

  detail <- list()

  rows <- lapply(attr_names, function(a) {
    levs  <- levels(data[[a]])
    pairs <- list()
    for (q in seq_along(levs)) for (p in seq_along(levs)) {
      if (q >= p) next
      pairs[[length(pairs) + 1]] <- c(levs[q], levs[p])
    }
    pair_labels <- vapply(pairs, function(pr) paste0(pr[1], " vs ", pr[2]), character(1))

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

    split_vals    <- rep(NA_real_, n_splits)
    selected_pair <- rep(NA_character_, n_splits)

    for (s in seq_len(n_splits)) {
      ids_a <- sample(unique_ids, size = n_a)
      in_a  <- d_a[[id_var]] %in% ids_a
      dA <- d_a[in_a, , drop = FALSE]
      dB <- d_a[!in_a, , drop = FALSE]

      sel_vals <- vapply(pairs, function(pr) {
        abs(.point_diff_pairwise(dA, outcome_var, a, pr[1], pr[2],
                                  id_var = id_var, task_var = task_var,
                                  informative = informative, profile_var = profile_var,
                                  weights_col = weights_col))
      }, numeric(1))
      if (all(is.na(sel_vals))) next
      winner <- which.max(sel_vals)
      pr     <- pairs[[winner]]
      selected_pair[s] <- pair_labels[winner]

      split_vals[s] <- abs(.point_diff_pairwise(dB, outcome_var, a, pr[1], pr[2],
                                                 id_var = id_var, task_var = task_var,
                                                 informative = informative,
                                                 profile_var = profile_var,
                                                 weights_col = weights_col))
    }

    sel_tab <- table(selected_pair[!is.na(selected_pair)])
    sel_p   <- sel_tab / sum(sel_tab)
    entropy <- if (length(sel_p) > 0) -sum(sel_p * log2(sel_p)) else NA_real_

    detail[[a]] <<- data.frame(item = a, split_estimate = split_vals,
                                selected_pair = selected_pair,
                                stringsAsFactors = FALSE)

    data.frame(
      item               = a,
      plugin_mapns       = plugin_mapns,
      split_median       = stats::median(split_vals, na.rm = TRUE),
      split_mean         = mean(split_vals, na.rm = TRUE),
      n_valid_splits     = sum(!is.na(split_vals)),
      n_splits           = n_splits,
      selection_entropy  = entropy,
      stringsAsFactors   = FALSE
    )
  })

  out <- do.call(rbind, rows)
  attr(out, "detail") <- do.call(rbind, detail)
  out
}
