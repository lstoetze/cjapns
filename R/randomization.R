# ══════════════════════════════════════════════════════════════════════════════
# Randomization test for attribute irrelevance (Algorithm 2 in the paper)
# ══════════════════════════════════════════════════════════════════════════════
#
# The parametric bootstrap behind cj_apns(se = "parametric") propagates
# uncertainty through |.| and max(.), but because those maps are non-negative
# its intervals structurally exclude zero and cannot test a null. This file
# supplies the separate test: re-draw each informative task's orientation --
# which of the two profiles carries which level of attribute l -- from a fair
# coin, hold everything else fixed (the unordered pair shown, all other
# attributes, preference groups, and the outcomes), and recompute the same
# plug-in estimator. The resulting reference distribution is exact in finite
# samples for any statistic, since it re-draws under the true assignment
# mechanism; respondent clustering is preserved automatically because only
# task-level orientation moves.


#' Build the task-level orientation table for one attribute
#'
#' Collapses the long (one row per profile) data to one row per
#' \eqn{(t_q,t_p)}-informative task, recording the unordered level pair shown,
#' the realized orientation \eqn{D_{ik}}, and both profiles' outcomes and
#' design weights. Everything the randomization test needs is here: a draw
#' re-randomizes only the \code{D} column, and the estimator is recomputed
#' from the fixed \code{y}/\code{w} columns.
#'
#' Restricted to two-profile designs (\eqn{J = 2}), where "swap the levels of
#' attribute \eqn{l} between the two profiles" is a single binary orientation.
#' A task is informative for attribute \eqn{l} exactly when its two profiles
#' show different levels of \eqn{l}, which matches
#' \code{\link{.filter_informative}} at \eqn{J = 2}.
#'
#' The level pair is canonicalized by the attribute's factor level order:
#' \eqn{t_q} is whichever of the two levels comes first, so the \code{pair}
#' label matches the \code{"<tq> vs <tp>"} labels used by \code{\link{cj_apns}}
#' and \eqn{D = 1} means profile 1 shows \eqn{t_q}.
#'
#' @param data A data.frame in long format (one row per profile).
#' @param outcome_var,attribute Character column names.
#' @param id_var,task_var Character names of the respondent-ID and task
#'   columns. Both required: tasks are what gets re-randomized.
#' @param profile_var Character name of the profile indicator, or `NULL`.
#'   When supplied it fixes which row of a task is profile 1; otherwise the
#'   data's own row order within a task is used.
#' @param levs Character vector of the attribute's levels, in factor order.
#' @param weights_col Character name of a design-weight column, or `NULL`.
#' @param context_vars Character names of the remaining attributes, stored per
#'   profile for the infeasible-swap check.
#' @return A data.frame with one row per informative task: `task_key`, `id`,
#'   `pair`, `D`, `y1`, `y2`, `w1`, `w2`, and `ctx1`, `ctx2` (each profile's
#'   remaining-attribute context, used by the infeasible-swap check).
#' @keywords internal
.ri_task_table <- function(data, outcome_var, attribute, id_var, task_var,
                           profile_var, levs, weights_col = NULL,
                           context_vars = character(0)) {
  task_key <- paste(data[[id_var]], data[[task_var]], sep = ":::")

  ord <- if (!is.null(profile_var))
    order(task_key, as.character(data[[profile_var]])) else order(task_key)
  d  <- data[ord, , drop = FALSE]
  tk <- task_key[ord]

  sizes <- rle(tk)$lengths
  if (any(sizes != 2L))
    stop("The randomization test is implemented for two-profile designs only, ",
         "but ", sum(sizes != 2L), " task(s) do not have exactly 2 rows. ",
         "Check 'id', 'tasks', and 'profile'.", call. = FALSE)

  i1 <- seq.int(1L, nrow(d), by = 2L)
  i2 <- i1 + 1L

  lev1 <- as.character(d[[attribute]])[i1]
  lev2 <- as.character(d[[attribute]])[i2]
  li1  <- match(lev1, levs)
  li2  <- match(lev2, levs)

  keep <- !is.na(li1) & !is.na(li2) & li1 != li2
  if (!any(keep)) return(NULL)

  i1 <- i1[keep]; i2 <- i2[keep]
  li1 <- li1[keep]; li2 <- li2[keep]

  w <- if (!is.null(weights_col)) d[[weights_col]] else rep(1, nrow(d))

  ctx <- if (length(context_vars) > 0)
    do.call(paste, c(lapply(context_vars, function(v) as.character(d[[v]])), sep = "\r"))
  else rep("", nrow(d))

  data.frame(
    task_key = tk[i1],
    id       = as.character(d[[id_var]])[i1],
    # canonical pair label, matching cj_apns(): tq is the earlier factor level
    pair     = paste0(levs[pmin(li1, li2)], " vs ", levs[pmax(li1, li2)]),
    D        = as.integer(li1 < li2),
    y1       = as.numeric(d[[outcome_var]])[i1],
    y2       = as.numeric(d[[outcome_var]])[i2],
    w1       = w[i1],
    w2       = w[i2],
    ctx1     = ctx[i1],
    ctx2     = ctx[i2],
    stringsAsFactors = FALSE
  )
}


#' Recompute the signed pairwise contrast for every re-randomization draw
#'
#' The long-format plug-in contrast of \code{\link{estimate_pairwise}} is the
#' weighted mean outcome over profiles showing \eqn{t_q} minus the weighted
#' mean over profiles showing \eqn{t_p}, both taken over informative tasks.
#' Flipping a task's orientation moves that task's two profiles between the
#' two means, which is all this computes -- no data are rewritten and no
#' re-filtering happens per draw.
#'
#' @param tt Task table from \code{\link{.ri_task_table}}.
#' @param idx Integer row indices of \code{tt} entering this contrast.
#' @param Dm Orientation matrix, \code{nrow(tt)} x \code{B}.
#' @return Numeric vector of length \code{ncol(Dm)}; \code{NA} where the cell
#'   has no informative task (or zero total weight).
#' @keywords internal
.ri_delta_cols <- function(tt, idx, Dm) {
  B <- ncol(Dm)
  if (length(idx) == 0) return(rep(NA_real_, B))

  D  <- Dm[idx, , drop = FALSE]
  y1 <- tt$y1[idx]; y2 <- tt$y2[idx]
  w1 <- tt$w1[idx]; w2 <- tt$w2[idx]

  # D = 1: profile 1 carries t_q, so the t_q-mean takes (y1, w1); D = 0 swaps.
  Ytq <- D * y1 + (1 - D) * y2
  Ytp <- D * y2 + (1 - D) * y1
  Wtq <- D * w1 + (1 - D) * w2
  Wtp <- D * w2 + (1 - D) * w1

  num_q <- colSums(Wtq * Ytq); den_q <- colSums(Wtq)
  num_p <- colSums(Wtp * Ytp); den_p <- colSums(Wtp)

  out <- num_q / den_q - num_p / den_p
  out[!is.finite(out)] <- NA_real_
  out
}


#' Recompute a full test statistic for every re-randomization draw
#'
#' @param tt Task table from \code{\link{.ri_task_table}}.
#' @param cells List of cells, each \code{list(idx = , weight = )}. A cell is
#'   one pairwise contrast (separability) or one preference group's contrast
#'   at its own pair (conditional/heterogeneous).
#' @param agg \code{"max"} for the homogeneous MAPNS (the maximum is retaken
#'   on every draw, which is what makes the test account for pair selection)
#'   or \code{"wsum"} for the share-weighted conditional estimators.
#' @param Dm Orientation matrix, \code{nrow(tt)} x \code{B}.
#' @param ok Logical vector over \code{cells} marking the cells with
#'   informative-task coverage. Coverage is a property of which tasks exist,
#'   not of their orientation, so it is fixed across draws; it is computed
#'   once on the observed statistic and reused, which keeps the `"wsum"`
#'   renormalization (see \code{\link{.renormalized_mapns}}) identical in the
#'   observed statistic and in every draw.
#' @return Numeric vector of length \code{ncol(Dm)}.
#' @keywords internal
.ri_statistic <- function(tt, cells, agg, Dm, ok = NULL) {
  B <- ncol(Dm)
  M <- matrix(NA_real_, nrow = B, ncol = length(cells))
  for (j in seq_along(cells)) M[, j] <- .ri_delta_cols(tt, cells[[j]]$idx, Dm)
  A <- abs(M)

  if (is.null(ok)) ok <- !is.na(A[1, ])
  if (!any(ok)) return(rep(NA_real_, B))
  A <- A[, ok, drop = FALSE]

  if (agg == "max") {
    apply(A, 1, max)
  } else {
    wts <- vapply(cells, function(ce) ce$weight, numeric(1))[ok]
    total <- sum(wts)
    if (!is.finite(total) || total <= 0) return(rep(NA_real_, B))
    as.numeric(A %*% (wts / total))
  }
}


#' Run one randomization test
#'
#' @param tt,cells,agg As in \code{\link{.ri_statistic}}.
#' @param B Number of re-randomization draws.
#' @return A list with `statistic`, `p_value`, `draws`, `n_informative`.
#' @keywords internal
.ri_run <- function(tt, cells, agg, B) {
  used <- sort(unique(unlist(lapply(cells, `[[`, "idx"))))
  n    <- nrow(tt)

  obs_all <- {
    Dobs <- matrix(tt$D, ncol = 1)
    M <- matrix(NA_real_, nrow = 1, ncol = length(cells))
    for (j in seq_along(cells)) M[, j] <- .ri_delta_cols(tt, cells[[j]]$idx, Dobs)
    abs(M)
  }
  ok  <- !is.na(obs_all[1, ])
  obs <- .ri_statistic(tt, cells, agg, matrix(tt$D, ncol = 1), ok = ok)

  if (is.na(obs))
    return(list(statistic = NA_real_, p_value = NA_real_, draws = NULL,
                n_informative = length(used)))

  draws <- numeric(B)
  # Chunk the draws so the orientation matrix stays a few million cells at
  # most, whatever B and the number of informative tasks are.
  chunk <- max(1L, as.integer(floor(2e6 / max(n, 1))))
  done  <- 0L
  while (done < B) {
    b  <- min(chunk, B - done)
    Dm <- matrix(stats::rbinom(n * b, 1L, 0.5), nrow = n, ncol = b)
    draws[(done + 1L):(done + b)] <- .ri_statistic(tt, cells, agg, Dm, ok = ok)
    done <- done + b
  }

  # Ties with the observed statistic are real (a flipped configuration can
  # reproduce the observed partition) and must count, so compare with a
  # tolerance rather than letting floating-point noise drop them, which would
  # make the test anti-conservative.
  tol <- 1e-10 * max(1, abs(obs))
  list(statistic = obs,
       p_value = (1 + sum(draws >= obs - tol, na.rm = TRUE)) / (B + 1),
       draws = draws,
       n_informative = length(used))
}


#' Map respondents to preference groups and group shares
#'
#' Mirrors the grouping that \code{\link{cj_apns}} applies internally, so the
#' observed test statistic equals the estimate \code{cj_apns()} reports.
#' Returns `NULL` when the attribute has no usable preference information.
#' @keywords internal
.ri_pref_groups <- function(data, id_var, attribute, levs, preferences) {
  if (!attribute %in% preferences$attributes) return(NULL)
  pref <- preferences$data
  pref <- pref[!duplicated(pref[[preferences$id_var]]), , drop = FALSE]
  pref_ids <- as.character(pref[[preferences$id_var]])

  if (preferences$type == "ranking") {
    top <- as.character(pref[[paste0(attribute, "_top")]])
    bot <- as.character(pref[[paste0(attribute, "_bot")]])
    lab <- paste0(top, " vs ", bot)
    # pi_g over respondents in the preferences object, as in .estimate_point()
    groups <- unique(lab[!is.na(top) & !is.na(bot)])
    if (length(groups) == 0) return(NULL)
    info <- lapply(groups, function(g) {
      j  <- which(lab == g)[1]
      ti <- match(top[j], levs); bi <- match(bot[j], levs)
      list(label = g, pi = mean(!is.na(lab) & lab == g),
           # the group's extreme pair, relabelled in canonical factor order so
           # it matches the task table's `pair` column
           pair = if (is.na(ti) || is.na(bi)) NA_character_
                  else paste0(levs[min(ti, bi)], " vs ", levs[max(ti, bi)]))
    })
    names(info) <- groups
    return(list(type = "ranking", info = info,
                map = stats::setNames(lab, pref_ids)))
  }

  vals <- pref[[attribute]]
  map  <- stats::setNames(vals, pref_ids)
  # pi_g over profile rows with non-missing preference, as in .estimate_point()
  row_g <- map[as.character(data[[id_var]])]
  has   <- !is.na(row_g)
  if (!any(has)) return(NULL)

  if (is.character(vals)) {
    groups <- sort(unique(vals[!is.na(vals)]))
    info <- lapply(groups, function(g)
      list(label = g, pi = mean(row_g[has] == g), pair = NA_character_))
  } else {
    pi_val <- mean(as.numeric(row_g[has]))
    groups <- c("1", "0")
    info <- list(list(label = "1", pi = pi_val, pair = NA_character_),
                 list(label = "0", pi = 1 - pi_val, pair = NA_character_))
  }
  names(info) <- groups
  list(type = "subclass", info = info,
       map = stats::setNames(as.character(vals), pref_ids))
}


# ══════════════════════════════════════════════════════════════════════════════
# Design checks
# ══════════════════════════════════════════════════════════════════════════════

#' Design checks for the e = 1/2 orientation assumption
#'
#' The test re-draws orientation from a fair coin, which is the design
#' probability only when attribute \eqn{l} is assigned independently of the
#' other attributes and the two profile slots are drawn from the same
#' distribution. These three checks look for the ways that can fail. They
#' warn rather than fail: a flagged design is not necessarily invalid, but
#' its \eqn{e_{ik}} should be computed from the design rather than fixed at
#' 1/2 (the package implements only the 1/2 case).
#'
#' \describe{
#'   \item{`balance`}{Per pair, an exact binomial test of \eqn{\#\{D=1\}}
#'     against \eqn{n/2}. A strong imbalance suggests position-dependent
#'     randomization.}
#'   \item{`independence`}{Chi-square test of attribute \eqn{l} against each
#'     other attribute at the profile level. Dependence implies restrictions
#'     or weighted joint draws, under which \eqn{e_{ik} \neq 1/2}.}
#'   \item{`swaps`}{Share of informative tasks whose swapped configuration --
#'     the same two remaining-attribute vectors with the levels of \eqn{l}
#'     exchanged -- is never observed. Under restrictions such tasks carry no
#'     randomization information. A raw share is not interpretable on its own,
#'     because a configuration seen once will usually not have its swap seen
#'     either, however unrestricted the design: with many other attributes
#'     almost every configuration is unique and the raw share approaches 1.
#'     The share is therefore calibrated against the same orientation re-draw
#'     the test itself uses: flipping which member of each configuration pair
#'     was realized gives the share expected under no restriction, and the
#'     reported `p_value` is the share's position in that reference
#'     distribution. Pure sparsity moves the observed and simulated shares
#'     together and leaves `p_value` large; a restriction pins every task to
#'     one member of its pair and makes the observed share stand out.
#'     Past a point sparsity leaves nothing to detect: once the expected share
#'     approaches 1, every configuration is unique whatever the design, the
#'     reference distribution collapses onto a point, and a difference of no
#'     practical size can still land in its tail. `informative` is `FALSE`
#'     there --- the check abstains rather than reporting a significant
#'     `p_value` it cannot back --- and no warning is raised.}
#' }
#'
#' @param tt Task table from \code{\link{.ri_task_table}}.
#' @param data,attr_names,attribute As in \code{\link{cj_ri_test}}.
#' @param n_sim Re-draws used to calibrate the infeasible-swap share.
#' @return A list with `balance`, `independence`, `swaps`.
#' @keywords internal
.ri_checks <- function(tt, data, attr_names, attribute, n_sim = 200) {
  # ── orientation balance ────────────────────────────────────────────────
  balance <- do.call(rbind, lapply(sort(unique(tt$pair)), function(pr) {
    d1 <- tt$D[tt$pair == pr]
    n  <- length(d1)
    pv <- if (n > 0) stats::binom.test(sum(d1), n, 0.5)$p.value else NA_real_
    data.frame(item = attribute, pair = pr, n = n, n_D1 = sum(d1),
               p_value = pv, stringsAsFactors = FALSE)
  }))
  row.names(balance) <- NULL

  # ── independence of attribute l from the other attributes ──────────────
  others <- setdiff(attr_names, attribute)
  independence <- if (length(others) == 0) NULL else
    do.call(rbind, lapply(others, function(o) {
      tab <- table(as.character(data[[attribute]]), as.character(data[[o]]))
      pv <- NA_real_
      if (nrow(tab) > 1 && ncol(tab) > 1)
        pv <- suppressWarnings(stats::chisq.test(tab)$p.value)
      data.frame(item = attribute, other = o, p_value = pv,
                 stringsAsFactors = FALSE)
    }))
  if (!is.null(independence)) row.names(independence) <- NULL

  # ── infeasible swaps ───────────────────────────────────────────────────
  # A configuration is the unordered set of (level of l, remaining attributes)
  # over the two profiles; its swap exchanges the levels of l between the two
  # remaining-attribute vectors. Keys are built so that the unordered set maps
  # to one string.
  tq <- sub(" vs .*$", "", tt$pair)
  tp <- sub("^.* vs ", "", tt$pair)
  lv1 <- ifelse(tt$D == 1L, tq, tp)   # level carried by profile 1
  lv2 <- ifelse(tt$D == 1L, tp, tq)
  key <- function(a, b) {
    lo <- pmin(a, b); hi <- pmax(a, b)
    paste(lo, hi, sep = "\n")
  }
  cfg  <- key(paste(lv1, tt$ctx1, sep = "\t"), paste(lv2, tt$ctx2, sep = "\t"))
  swp  <- key(paste(lv2, tt$ctx1, sep = "\t"), paste(lv1, tt$ctx2, sep = "\t"))

  # A raw share of unobserved swaps is uninterpretable: a configuration seen
  # once will usually not have its swap seen either, whatever the design. So
  # calibrate it against the null of no restriction, under which the realized
  # member of each {configuration, swap} pair is a coin flip -- the same
  # re-draw the test itself performs, applied to the full configuration key.
  obs_share <- mean(!(swp %in% cfg))
  sim <- vapply(seq_len(n_sim), function(b) {
    coin     <- stats::rbinom(length(cfg), 1L, 0.5) == 1L
    realized <- ifelse(coin, cfg, swp)
    other    <- ifelse(coin, swp, cfg)
    mean(!(other %in% realized))
  }, numeric(1))

  exp_share <- mean(sim)
  swaps <- list(
    share_infeasible = obs_share,
    share_expected   = exp_share,
    p_value          = (1 + sum(sim >= obs_share - 1e-12)) / (n_sim + 1),
    # With little room left between the expected share and 1, there is nothing
    # the check could detect: every configuration is unique whatever the
    # design. Abstain rather than report a p-value the data cannot support.
    informative      = exp_share < 0.9,
    n_distinct_configs = length(unique(cfg)),
    n_tasks          = nrow(tt))

  list(balance = balance, independence = independence, swaps = swaps)
}


# ══════════════════════════════════════════════════════════════════════════════
# Main entry point
# ══════════════════════════════════════════════════════════════════════════════

#' Randomization Test for Attribute Irrelevance
#'
#' Tests the sharp null that swapping the levels of an attribute between the
#' two profiles of a task changes no respondent's choice — the null of
#' attribute irrelevance (Algorithm 2 in Stoetzer and Magazinnik 2026). For
#' the APNS the null concerns one level pair; for the MAPNS, any pair. Under
#' the conditional and heterogeneous assumptions the estimands are weighted
#' averages of non-negative group terms, so a zero aggregate implies a zero
#' term in every group and the same null applies.
#'
#' This complements, rather than replaces, the intervals from
#' \code{cj_apns(se = "parametric")}. Those propagate sampling uncertainty
#' through \eqn{|\cdot|} and \eqn{\max(\cdot)}, but because both maps are
#' non-negative the resulting intervals structurally exclude zero and cannot
#' be read as tests of a null. Use the intervals for uncertainty and this
#' function for testing.
#'
#' @details
#' The procedure keeps the informative tasks — those whose two profiles show
#' different levels of the attribute — and records each one's orientation
#' \eqn{D_{ik} = 1} if profile 1 carries \eqn{t_q} (the pair's first level in
#' the attribute's factor order). It then re-draws every \eqn{D_{ik}} from a
#' fair coin, holding fixed the unordered pair shown, all other attributes,
#' preference-group memberships, and the outcomes, and recomputes the same
#' plug-in estimator from scratch. The \eqn{p}-value is
#' \eqn{(1 + \#\{S^{(b)} \ge S^{\mathrm{obs}}\})/(B+1)}; the test is
#' one-sided because every statistic is non-negative.
#'
#' The statistic recomputed per draw is the full estimator:
#' \describe{
#'   \item{APNS, separability}{\eqn{|\hat\Delta|} for the pair.}
#'   \item{APNS, conditional}{Group-share-weighted sum of the within-group
#'     absolute contrasts, with shares fixed at their observed values.}
#'   \item{MAPNS, separability}{The maximum over pairs. The maximum is
#'     retaken on every draw, which is what makes the test account for the
#'     selection of the largest pairwise contrast.}
#'   \item{MAPNS, conditional}{The share-weighted sum over extreme-pair
#'     groups of each group's absolute contrast at its own extreme pair
#'     (`preferences` of `type = "ranking"`), or the conditional APNS at the
#'     single pair when the attribute is binary.}
#' }
#' Because only the orientation vector is re-drawn and outcomes are held
#' fixed, the clustering of tasks within respondents is preserved and the
#' test is exact in finite samples.
#'
#' The orientation probability is fixed at \eqn{e_{ik} = 1/2}. This is the
#' design probability whenever the attribute is assigned independently of the
#' other attributes and the two profile slots are drawn from the same
#' distribution — including under non-uniform level probabilities, which
#' change how many informative tasks each pair receives but not the
#' orientation balance within them. It departs from 1/2 only when the
#' attribute is correlated with the others within a profile (restrictions or
#' weighted joint draws) or the two slots are drawn differently. `checks`
#' looks for both; see \code{\link{.ri_checks}}.
#'
#' Target (context) weights affect the statistic but not the reference
#' distribution, so `design` changes the estimator being tested and the
#' test's power, never its size.
#'
#' Implemented for two-profile designs only.
#'
#' @param formula A formula: `outcome ~ attr1 + attr2 + ...`, as in
#'   \code{\link{cj_apns}}.
#' @param data A data.frame in long format (one row per profile).
#' @param id A one-sided formula for the respondent ID (e.g. `~ ResponseId`).
#' @param tasks A one-sided formula for the task-number variable (e.g.
#'   `~ time`). Required: tasks are the unit being re-randomized.
#' @param profile A one-sided formula for the profile indicator (e.g.
#'   `~ profile`). Recommended; without it the data's row order within a task
#'   determines which profile is profile 1.
#' @param estimand `"mapns"` (default) or `"apns"`.
#' @param assumption `"separability"` (default), `"conditional"`, or
#'   `"both"`.
#' @param preferences A `"cj_preferences"` object from
#'   \code{\link{make_preferences}}. Required when `assumption` is
#'   `"conditional"` or `"both"`.
#' @param attributes Character vector of attributes to test. Default `NULL`
#'   tests every attribute on the right-hand side of `formula`.
#' @param pair For `estimand = "apns"`, a length-2 character vector
#'   `c(tq, tp)` naming the level pair to test. Default `NULL` tests every
#'   observed pair. Ignored for `estimand = "mapns"`.
#' @param design `"uniform"` (default) or a `"cj_design"` object from
#'   \code{\link{make_design}}, applying the same context weights
#'   \code{\link{cj_apns}} would.
#' @param B Number of re-randomization draws. Default 2000.
#' @param seed Integer seed, default 123. The caller's RNG stream is saved
#'   and restored on exit, so the call is reproducible without fixing the
#'   randomness of anything run after it. `NULL` draws from the ambient
#'   stream.
#' @param checks Logical, default `TRUE`: run the design checks described in
#'   \code{\link{.ri_checks}} and warn on violations.
#' @param min_informative Warn when a preference group contributes fewer than
#'   this many informative tasks to a conditional statistic. Default 10.
#' @param return_draws Logical, default `FALSE`. When `TRUE`, the \eqn{B}
#'   re-randomized statistics are returned for each test (for plotting the
#'   reference distribution).
#'
#' @return An object of class `"cj_ri_test"`:
#'   \describe{
#'     \item{`results`}{data.frame with one row per test: `item`,
#'       `estimand`, `assumption`, `pair` (`NA` for MAPNS), `statistic`
#'       (\eqn{S^{\mathrm{obs}}}), `p_value`, `n_informative`, `B`.}
#'     \item{`draws`}{Named list of length-`B` numeric vectors, or `NULL`
#'       unless `return_draws = TRUE`.}
#'     \item{`checks`}{Per-attribute list of `balance`, `independence`, and
#'       `swaps`, or `NULL` when `checks = FALSE`.}
#'     \item{`skipped`}{data.frame of tests that could not be run, with the
#'       reason. `NULL` if none.}
#'   }
#'
#' @references
#' Stoetzer, L.F. and Magazinnik, A. (2026). Measuring Attribute Relevance
#' in Conjoint Analysis.
#'
#' @seealso \code{\link{cj_apns}} for estimation and confidence intervals,
#'   \code{\link{mapns_split_test}} for the winner's-curse diagnostic.
#'
#' @examples
#' \dontrun{
#' data(cnj_cand)
#'
#' # MAPNS under separable monotonicity, every attribute
#' rt <- cj_ri_test(vote ~ borders + eurobonds + immucard + schools,
#'                  data = cnj_cand, id = ~ ResponseId,
#'                  tasks = ~ time, profile = ~ profile, B = 2000)
#' rt
#'
#' # One pairwise APNS
#' cj_ri_test(vote ~ borders, data = cnj_cand, id = ~ ResponseId,
#'            tasks = ~ time, profile = ~ profile,
#'            estimand = "apns", pair = c("open", "closed"))
#' }
#'
#' @export
cj_ri_test <- function(formula, data, id, tasks, profile = NULL,
                       estimand = c("mapns", "apns"),
                       assumption = c("separability", "conditional", "both"),
                       preferences = NULL,
                       attributes = NULL,
                       pair = NULL,
                       design = "uniform",
                       B = 2000, seed = 123,
                       checks = TRUE, min_informative = 10,
                       return_draws = FALSE) {

  cl <- match.call()
  estimand   <- match.arg(estimand)
  assumption <- match.arg(assumption)
  stopifnot(is.numeric(B), length(B) == 1, B >= 1)
  B <- as.integer(B)

  if (missing(id))    stop("'id' must be specified (e.g., id = ~ ResponseId).")
  if (missing(tasks) || is.null(tasks))
    stop("'tasks' must be specified (e.g., tasks = ~ time): the randomization ",
         "test re-draws the orientation of informative tasks, which cannot be ",
         "identified without a task variable.", call. = FALSE)

  id_var <- all.vars(id)
  stopifnot(length(id_var) == 1, id_var %in% names(data))
  task_var <- all.vars(tasks)[1]
  stopifnot("'tasks' variable not found in data" = task_var %in% names(data))
  profile_var <- if (!is.null(profile)) all.vars(profile)[1] else NULL
  if (!is.null(profile_var))
    stopifnot("'profile' variable not found in data" = profile_var %in% names(data))

  if (assumption %in% c("conditional", "both") && is.null(preferences))
    stop("'preferences' must be provided for assumption = \"", assumption,
         "\". Use make_preferences().", call. = FALSE)

  # Same RNG etiquette as cj_apns(): reproducible by default, and the
  # caller's stream is left exactly as we found it.
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1, !is.na(seed))
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
    }
    set.seed(seed)
  }

  tt_terms   <- stats::terms(formula, data = data)
  attr_names <- attr(tt_terms, "term.labels")
  attr_names <- attr_names[!grepl(":", attr_names)]
  outcome_var <- all.vars(formula)[1]
  for (a in attr_names)
    if (!is.factor(data[[a]])) data[[a]] <- as.factor(data[[a]])

  if (is.null(attributes)) attributes <- attr_names
  unknown <- setdiff(attributes, attr_names)
  if (length(unknown) > 0)
    stop("Attribute(s) not in 'formula': ", paste(unknown, collapse = ", "),
         call. = FALSE)

  assumptions <- if (assumption == "both") c("separability", "conditional")
                 else assumption

  rows <- list(); draw_store <- list(); check_store <- list(); skipped <- list()
  thin_groups <- character(0)

  for (a in attributes) {
    levs <- levels(data[[a]])

    weights_col <- NULL
    d_a <- data
    if (!identical(design, "uniform")) {
      d_a$.design_weight <- .attach_design_weights(d_a, attr_names, a, design,
                                                   id_var, task_var)
      weights_col <- ".design_weight"
    }

    tt <- .ri_task_table(d_a, outcome_var, a, id_var, task_var, profile_var,
                         levs, weights_col,
                         context_vars = setdiff(attr_names, a))
    if (is.null(tt)) {
      skipped[[length(skipped) + 1L]] <- data.frame(
        item = a, estimand = estimand, assumption = assumption,
        reason = "no informative tasks", stringsAsFactors = FALSE)
      next
    }

    observed_pairs <- sort(unique(tt$pair))
    target_pairs <- observed_pairs
    if (estimand == "apns" && !is.null(pair)) {
      stopifnot("'pair' must be a length-2 character vector" = length(pair) == 2)
      pi_ <- match(as.character(pair), levs)
      if (anyNA(pi_))
        stop("'pair' levels not found in attribute '", a, "': ",
             paste(pair, collapse = ", "), call. = FALSE)
      target_pairs <- paste0(levs[min(pi_)], " vs ", levs[max(pi_)])
      if (!target_pairs %in% observed_pairs)
        stop("No informative tasks for pair '", target_pairs, "' on attribute '",
             a, "'.", call. = FALSE)
    }

    for (asm in assumptions) {
      specs <- list()

      if (asm == "separability") {
        if (estimand == "mapns") {
          specs[[1]] <- list(
            pair = NA_character_, agg = "max",
            cells = lapply(observed_pairs, function(pr)
              list(idx = which(tt$pair == pr), weight = 1)))
        } else {
          for (pr in target_pairs)
            specs[[length(specs) + 1L]] <- list(
              pair = pr, agg = "max",
              cells = list(list(idx = which(tt$pair == pr), weight = 1)))
        }

      } else {
        pg <- .ri_pref_groups(d_a, id_var, a, levs, preferences)
        if (is.null(pg)) {
          skipped[[length(skipped) + 1L]] <- data.frame(
            item = a, estimand = estimand, assumption = asm,
            reason = paste0("no usable preference information for this ",
                            "attribute in 'preferences'"),
            stringsAsFactors = FALSE)
          next
        }
        grp <- pg$map[tt$id]

        if (pg$type == "ranking") {
          if (estimand != "mapns") {
            skipped[[length(skipped) + 1L]] <- data.frame(
              item = a, estimand = estimand, assumption = asm,
              reason = paste0("preferences of type \"ranking\" identify the ",
                              "conditional estimator only at each group's own ",
                              "extreme pair; use estimand = \"mapns\""),
              stringsAsFactors = FALSE)
            next
          }
          cells <- lapply(pg$info, function(gi) list(
            idx = if (is.na(gi$pair)) integer(0)
                  else which(!is.na(grp) & grp == gi$label & tt$pair == gi$pair),
            weight = gi$pi, label = gi$label))
          specs[[1]] <- list(pair = NA_character_, agg = "wsum", cells = cells)

        } else {
          if (estimand == "mapns" && length(levs) > 2) {
            skipped[[length(skipped) + 1L]] <- data.frame(
              item = a, estimand = estimand, assumption = asm,
              reason = paste0("conditional MAPNS is not identified for an ",
                              "attribute with more than two levels unless ",
                              "preferences are of type \"ranking\" (each ",
                              "group's own extreme pair is needed)"),
              stringsAsFactors = FALSE)
            next
          }
          # Binary attribute: the conditional MAPNS is the conditional APNS
          # at the attribute's single pair (Proposition 6).
          prs <- if (estimand == "mapns") observed_pairs else target_pairs
          for (pr in prs) {
            cells <- lapply(pg$info, function(gi) list(
              idx = which(!is.na(grp) & grp == gi$label & tt$pair == pr),
              weight = gi$pi, label = gi$label))
            specs[[length(specs) + 1L]] <- list(
              pair = if (estimand == "mapns") NA_character_ else pr,
              agg = "wsum", cells = cells)
          }
        }
      }

      for (sp in specs) {
        if (asm == "conditional") {
          thin <- vapply(sp$cells, function(ce)
            length(ce$idx) > 0 && length(ce$idx) < min_informative, logical(1))
          if (any(thin))
            thin_groups <- c(thin_groups, paste0(
              a, " [", paste(vapply(sp$cells[thin], function(ce)
                as.character(ce$label), character(1)), collapse = ", "), "]"))
        }

        res <- .ri_run(tt, sp$cells, sp$agg, B)
        key <- paste(a, estimand, asm,
                     if (is.na(sp$pair)) "" else sp$pair, sep = " | ")
        if (return_draws && !is.null(res$draws)) draw_store[[key]] <- res$draws

        rows[[length(rows) + 1L]] <- data.frame(
          item = a, estimand = estimand, assumption = asm, pair = sp$pair,
          statistic = res$statistic, p_value = res$p_value,
          n_informative = res$n_informative, B = B,
          stringsAsFactors = FALSE)
      }
    }

    # After the tests, so that turning `checks` on or off cannot shift the
    # reported p-values by consuming from the RNG stream ahead of them.
    if (checks) check_store[[a]] <- .ri_checks(tt, d_a, attr_names, a)
  }

  results <- if (length(rows) > 0) do.call(rbind, rows) else
    data.frame(item = character(0), estimand = character(0),
               assumption = character(0), pair = character(0),
               statistic = numeric(0), p_value = numeric(0),
               n_informative = integer(0), B = integer(0),
               stringsAsFactors = FALSE)
  row.names(results) <- NULL

  skipped_df <- if (length(skipped) > 0) {
    s <- do.call(rbind, skipped); row.names(s) <- NULL; s
  } else NULL

  if (length(thin_groups) > 0)
    warning("Preference group(s) with fewer than ", min_informative,
            " informative tasks contribute to a conditional statistic: ",
            paste(unique(thin_groups), collapse = "; "),
            ". Their contrasts are estimated from very little data.",
            call. = FALSE)

  if (checks && length(check_store) > 0) .ri_warn_checks(check_store)

  structure(list(results = results,
                 draws = if (return_draws) draw_store else NULL,
                 checks = if (checks) check_store else NULL,
                 skipped = skipped_df,
                 B = B, call = cl),
            class = "cj_ri_test")
}


#' Emit consolidated warnings from the design checks
#' @keywords internal
.ri_warn_checks <- function(check_store) {
  bal <- do.call(rbind, lapply(check_store, `[[`, "balance"))
  if (!is.null(bal)) {
    bad <- bal[!is.na(bal$p_value) & bal$p_value < 0.01, , drop = FALSE]
    if (nrow(bad) > 0)
      warning("Orientation imbalance (binomial test against 1/2, p < 0.01) for: ",
              paste0(bad$item, " [", bad$pair, "]", collapse = "; "),
              ". The test fixes e_ik = 1/2; a position-dependent randomizer ",
              "would need e_ik computed from the design.", call. = FALSE)
  }

  ind <- do.call(rbind, lapply(check_store, `[[`, "independence"))
  if (!is.null(ind)) {
    bad <- ind[!is.na(ind$p_value) & ind$p_value < 0.01, , drop = FALSE]
    if (nrow(bad) > 0)
      warning("Attribute levels are not independent across attributes ",
              "(chi-square, p < 0.01): ",
              paste0(bad$item, " ~ ", bad$other, collapse = "; "),
              ". Restrictions or weighted joint draws make e_ik != 1/2.",
              call. = FALSE)
  }

  swp <- vapply(check_store, function(x) x$swaps$p_value, numeric(1))
  obs <- vapply(check_store, function(x) x$swaps$share_infeasible, numeric(1))
  exp_ <- vapply(check_store, function(x) x$swaps$share_expected, numeric(1))
  inf_ <- vapply(check_store, function(x) isTRUE(x$swaps$informative), logical(1))
  bad <- names(swp)[inf_ & !is.na(swp) & swp < 0.01]
  if (length(bad) > 0)
    warning("More informative tasks lack an observed swapped configuration ",
            "than sparsity alone explains (observed vs expected share): ",
            paste0(bad, " = ", round(obs[bad], 3), " vs ",
                   round(exp_[bad], 3), collapse = "; "),
            ". This is the signature of a restricted design, under which ",
            "those tasks carry no randomization information and e_ik is not ",
            "1/2.", call. = FALSE)
  invisible(NULL)
}


#' @export
print.cj_ri_test <- function(x, ...) {
  cat("\nRandomization Test for Attribute Irrelevance\n")
  cat("============================================\n")
  cat("Null: swapping the attribute's levels between the two profiles\n")
  cat("      changes no respondent's choice.\n")
  cat("Draws (B):", x$B, " | orientation probability e = 1/2\n\n")

  if (nrow(x$results) == 0) {
    cat("No tests could be run.\n")
  } else {
    out <- x$results[, c("item", "estimand", "assumption", "pair",
                         "statistic", "p_value", "n_informative")]
    out$pair      <- ifelse(is.na(out$pair), "-", out$pair)
    out$statistic <- round(out$statistic, 4)
    out$p_value   <- format.pval(out$p_value, digits = 3, eps = 1 / (x$B + 1))
    print(out, row.names = FALSE)
  }

  if (!is.null(x$skipped) && nrow(x$skipped) > 0) {
    cat("\nSkipped:\n")
    for (i in seq_len(nrow(x$skipped)))
      cat("  ", x$skipped$item[i], " (", x$skipped$assumption[i], "): ",
          x$skipped$reason[i], "\n", sep = "")
  }
  cat("\n")
  invisible(x)
}
