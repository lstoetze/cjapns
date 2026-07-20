# ══════════════════════════════════════════════════════════════════════════════
# Standard error methods for cj_apns
# ══════════════════════════════════════════════════════════════════════════════


# ── Parametric bootstrap ─────────────────────────────────────────────────────
#
# Each direct pairwise contrast (Theorem 1/2/4's plug-in estimator) is
# asymptotically normal. We draw B samples of each pairwise contrast actually
# needed directly from N(estimate, se^2) — evaluated at the exact pair/group,
# never reconstructed from base-level coefficients — apply the absolute-value
# / weighted plugin transformation, and take the SD of the resulting MAPNS
# distribution. This correctly propagates uncertainty through the |.| mapping.

#' Draw one simulated value from a pairwise estimate's sampling distribution
#' @keywords internal
.sim_pairwise <- function(info) {
  if (is.null(info) || is.na(info$estimate) || is.na(info$se)) return(NA_real_)
  stats::rnorm(1, info$estimate, info$se)
}

#' Precompute direct pairwise (point, se) for every contrast inference needs
#'
#' Evaluated directly at the exact pair/group required by Theorems 1-4
#' (never reconstructed from base-level coefficients). Shared by
#' `.se_parametric` and `.se_folded_normal` so the two SE methods read the
#' subclassing logic from one place and cannot drift out of sync on which
#' subsample each pairwise estimate is computed on.
#'
#' @return A list with `sep_info` (named by attribute, then by pair) and
#'   `cond_info` (named by attribute; each entry has `type` and `groups`,
#'   where `groups` is keyed by pair for ranking-type — direct
#'   `estimate_pairwise` output — or by pair then group for binary/multilevel).
#' @keywords internal
.precompute_pairwise_info <- function(data, outcome_var, attr_names, attributes_info,
                                      preferences, id_var, task_var, informative,
                                      profile_var, do_sep, do_cond, design = "uniform") {
  sep_info <- list()
  if (do_sep) {
    for (a in attr_names) {
      levs <- attributes_info[[a]]
      weights_col <- NULL
      if (!identical(design, "uniform")) {
        data$.design_weight <- .attach_design_weights(data, attr_names, a, design,
                                                       id_var, task_var)
        weights_col <- ".design_weight"
      }
      pairs <- list()
      for (q in seq_along(levs)) for (p in seq_along(levs)) {
        if (q >= p) next
        pair <- paste0(levs[q], " vs ", levs[p])
        pairs[[pair]] <- estimate_pairwise(
          data, outcome_var, a, levs[q], levs[p],
          id_var = id_var, task_var = task_var,
          informative = informative, profile_var = profile_var,
          weights_col = weights_col)
      }
      sep_info[[a]] <- pairs
    }
  }

  cond_info <- list()
  if (do_cond && !is.null(preferences)) {
    for (a in attr_names) {
      if (!a %in% preferences$attributes) next
      levs <- attributes_info[[a]]

      weights_col <- NULL
      if (!identical(design, "uniform")) {
        data$.design_weight <- .attach_design_weights(data, attr_names, a, design,
                                                       id_var, task_var)
        weights_col <- ".design_weight"
      }

      if (preferences$type == "ranking") {
        pref_data <- preferences$data[
          !duplicated(preferences$data[[preferences$id_var]]), , drop = FALSE]
        pref_ids  <- as.character(pref_data[[preferences$id_var]])
        row_ids   <- as.character(data[[id_var]])
        top_vals  <- as.character(pref_data[[paste0(a, "_top")]])
        bot_vals  <- as.character(pref_data[[paste0(a, "_bot")]])
        ep_label  <- paste0(top_vals, " vs ", bot_vals)
        ep_row    <- ep_label[match(row_ids, pref_ids)]

        ep_info <- list()
        for (ep in unique(ep_label[!is.na(ep_label)])) {
          parts <- strsplit(ep, " vs ")[[1]]
          dm_g  <- data[!is.na(ep_row) & ep_row == ep, , drop = FALSE]
          if (nrow(dm_g) == 0) next
          ep_info[[ep]] <- estimate_pairwise(
            dm_g, outcome_var, a, parts[1], parts[2],
            id_var = id_var, task_var = task_var,
            informative = informative, profile_var = profile_var,
            weights_col = weights_col)
        }
        cond_info[[a]] <- list(type = "ranking", groups = ep_info)

      } else {
        # binary / multilevel: attribute-level groups
        pref_df  <- preferences$data[!duplicated(preferences$data[[preferences$id_var]]),
                                     c(preferences$id_var, a), drop = FALSE]
        pref_map <- stats::setNames(pref_df[[a]], pref_df[[preferences$id_var]])
        data$.pg <- pref_map[as.character(data[[id_var]])]
        has <- !is.na(data$.pg)

        pair_info <- list()
        if (sum(has) > 0) {
          dm          <- data[has, , drop = FALSE]
          pg_sample   <- dm$.pg[!is.na(dm$.pg)]
          pref_groups <- if (is.character(pg_sample)) sort(unique(pg_sample)) else c("1", "0")
          dm$.pg_chr  <- as.character(dm$.pg)
          sub_dm_list <- stats::setNames(
            lapply(pref_groups, function(g) dm[dm$.pg_chr == g, , drop = FALSE]),
            pref_groups)

          for (q in seq_along(levs)) for (p in seq_along(levs)) {
            if (q >= p) next
            pair <- paste0(levs[q], " vs ", levs[p])
            grp_info <- list()
            for (g in pref_groups) {
              sub_dm <- sub_dm_list[[g]]
              if (nrow(sub_dm) == 0) next
              grp_info[[g]] <- estimate_pairwise(
                sub_dm, outcome_var, a, levs[q], levs[p],
                id_var = id_var, task_var = task_var,
                informative = informative, profile_var = profile_var,
                weights_col = weights_col)
            }
            pair_info[[pair]] <- grp_info
          }
        }
        data$.pg <- NULL
        cond_info[[a]] <- list(type = preferences$type, groups = pair_info)
      }
    }
  }

  list(sep_info = sep_info, cond_info = cond_info)
}

#' @keywords internal
.se_parametric <- function(formula, data, id, id_var, attr_names,
                           estimand, assumption, preferences, B = 500, alpha = 0.05,
                           task_var = NULL, informative = "all", profile_var = NULL,
                           design = "uniform") {

  pt <- .estimate_point(formula, data, id, id_var, attr_names,
                        estimand, assumption, preferences, task_var, informative,
                        profile_var, design)
  point_vec <- .flatten_estimates(pt, estimand, assumption, attr_names)

  attributes_info <- pt$attributes
  outcome_var <- all.vars(formula)[1]
  do_sep  <- assumption %in% c("separability", "both")
  do_cond <- assumption %in% c("conditional", "both")

  info <- .precompute_pairwise_info(data, outcome_var, attr_names, attributes_info,
                                    preferences, id_var, task_var, informative,
                                    profile_var, do_sep, do_cond, design)
  sep_info  <- info$sep_info
  cond_info <- info$cond_info

  # ── Draw B parametric bootstrap replicates ────────────────────────────────
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(point_vec))
  colnames(boot_mat) <- names(point_vec)

  for (b in seq_len(B)) {

    # AMCE (vs base level): resample directly from pt$amce's own (point, se),
    # which estimate_amce() already computes via estimate_pairwise().
    for (a in attr_names) {
      amce_a <- pt$amce[[a]]
      for (lev in names(amce_a$estimate)) {
        nm <- paste0("amce.", a, ".", lev)
        if (nm %in% colnames(boot_mat))
          boot_mat[b, nm] <- .sim_pairwise(list(estimate = amce_a$estimate[lev],
                                                se = amce_a$se[lev]))
      }
    }

    for (a in attr_names) {
      levs <- attributes_info[[a]]
      Dl <- length(levs)

      # Separability: simulate each pair directly, MAPNS = max (Prop. 3/7).
      if (do_sep) {
        sim_vals <- vapply(sep_info[[a]], .sim_pairwise, numeric(1))
        mapns_sep_b <- max(abs(sim_vals), na.rm = TRUE)
        nm <- paste0("mapns.separability.", a)
        if (nm %in% colnames(boot_mat)) boot_mat[b, nm] <- mapns_sep_b
      }

      # Conditional/heterogeneous: simulate each (pair, group) directly.
      if (do_cond && a %in% names(cond_info)) {
        info_a <- cond_info[[a]]
        mapns_cond_b <- NA_real_

        if (info_a$type == "ranking") {
          # MAPNS = sum_g pi_g |CAPNS_g| at each group's own extreme pair (Prop. 4/8).
          mapns_cond_b <- 0
          for (ep in names(info_a$groups)) {
            v_ep  <- pt$acmce[[a]][[ep]]
            sim_v <- .sim_pairwise(info_a$groups[[ep]])
            nm_g  <- paste0("acmce.", make.names(ep), ".", a)
            if (nm_g %in% colnames(boot_mat)) boot_mat[b, nm_g] <- sim_v
            if (!is.na(sim_v)) mapns_cond_b <- mapns_cond_b + v_ep$pi * abs(sim_v)
          }
        } else if (length(info_a$groups) > 0) {
          pi_info <- pt$pi_hat[[a]]
          is_multilevel_a <- length(pi_info) > 1
          vals_c <- c()

          for (pair in names(info_a$groups)) {
            grp_info <- info_a$groups[[pair]]

            if (!is_multilevel_a) {
              # binary: single scalar pi for all pairs
              sim_pro <- .sim_pairwise(grp_info[["1"]])
              sim_con <- .sim_pairwise(grp_info[["0"]])
              if (is.na(sim_pro)) sim_pro <- 0
              if (is.na(sim_con)) sim_con <- 0
              vals_c <- c(vals_c, pi_info * abs(sim_pro) + (1 - pi_info) * abs(sim_con))
              nm_pro <- paste0("acmce.pro.", a, ".", pair)
              nm_con <- paste0("acmce.con.", a, ".", pair)
              if (nm_pro %in% colnames(boot_mat)) boot_mat[b, nm_pro] <- sim_pro
              if (nm_con %in% colnames(boot_mat)) boot_mat[b, nm_con] <- sim_con
            } else {
              # multilevel: K groups, attribute-level pi
              sim_grp <- vapply(names(pi_info), function(g) {
                v <- .sim_pairwise(grp_info[[g]])
                if (is.na(v)) 0 else v
              }, numeric(1))
              vals_c <- c(vals_c, sum(pi_info * abs(sim_grp)))
              for (g in names(pi_info)) {
                nm_g <- paste0("acmce.", g, ".", a, ".", pair)
                if (nm_g %in% colnames(boot_mat)) boot_mat[b, nm_g] <- sim_grp[g]
              }
            }
          }
          # Proposition 8: MAPNS is the extreme-pair, group-share-weighted sum,
          # not an average over all pairs. For Dl == 2 the single pair is
          # every group's extreme pair, so this coincides with sum(vals_c).
          # For Dl > 2, non-ranking preferences do not identify each group's
          # own extreme pair, so MAPNS is not estimable here (see
          # .mapns_from_pairwise_cond).
          mapns_cond_b <- if (Dl == 2) sum(vals_c) else NA_real_
        }

        if (!is.na(mapns_cond_b)) {
          nm <- paste0("mapns.conditional.", a)
          if (nm %in% colnames(boot_mat)) boot_mat[b, nm] <- mapns_cond_b
        }
      }
    }
  }

  se_vec  <- apply(boot_mat, 2, stats::sd, na.rm = TRUE)
  # Algorithm 1 (paper): percentile-based CIs from the empirical bootstrap quantiles
  lower_vec <- apply(boot_mat, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE)
  upper_vec <- apply(boot_mat, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  list(se = se_vec, point = point_vec, lower = lower_vec, upper = upper_vec)
}


# ── Nonparametric bootstrap ──────────────────────────────────────────────────

#' @keywords internal
.se_bootstrap <- function(formula, data, id, id_var, attr_names,
                          estimand, assumption, preferences, B = 500,
                          task_var = NULL, informative = "all", profile_var = NULL,
                          design = "uniform") {

  pt <- .estimate_point(formula, data, id, id_var, attr_names,
                        estimand, assumption, preferences, task_var, informative,
                        profile_var, design)
  point_vec <- .flatten_estimates(pt, estimand, assumption, attr_names)

  unique_ids <- unique(data[[id_var]])
  n_ids <- length(unique_ids)

  boot_ests <- vector("list", B)
  for (b in seq_len(B)) {
    sampled <- sample(unique_ids, n_ids, replace = TRUE)
    bd <- do.call(rbind, lapply(seq_along(sampled), function(i) {
      rows <- data[data[[id_var]] == sampled[i], , drop = FALSE]
      rows[[id_var]] <- i; rows
    }))
    bp <- preferences
    if (!is.null(preferences)) {
      pd <- preferences$data
      bpd <- do.call(rbind, lapply(seq_along(sampled), function(i) {
        rows <- pd[pd[[preferences$id_var]] == sampled[i], , drop = FALSE]
        rows[[preferences$id_var]] <- i; rows
      }))
      bp <- preferences; bp$data <- bpd
    }
    bid <- stats::reformulate(id_var, response = NULL)
    tryCatch({
      est <- .estimate_point(formula, bd, bid, id_var, attr_names,
                             estimand, assumption, bp, task_var, informative,
                             profile_var, design)
      boot_ests[[b]] <- .flatten_estimates(est, estimand, assumption, attr_names)
    }, error = function(e) { boot_ests[[b]] <<- NULL })
  }

  boot_ests <- Filter(Negate(is.null), boot_ests)
  mat <- do.call(rbind, lapply(boot_ests, function(x) x[names(point_vec)]))
  se_vec <- apply(mat, 2, stats::sd, na.rm = TRUE)
  list(se = se_vec, point = point_vec)
}


# ── Folded normal approximation ──────────────────────────────────────────────
#
# For X ~ N(mu, sigma^2), |X| follows a folded normal with:
#   E[|X|] = sigma * sqrt(2/pi) * exp(-mu^2/(2*sigma^2)) + mu * (1 - 2*Phi(-mu/sigma))
#   Var(|X|) = mu^2 + sigma^2 - E[|X|]^2
# We use the AMCE point estimate and SE as (mu, sigma).

#' Direct pairwise SE for the folded-normal method, or 0 if unavailable
#' @keywords internal
.se_or_zero <- function(info) {
  if (is.null(info) || is.na(info$se)) 0 else info$se
}

#' @keywords internal
.se_folded_normal <- function(pt, assumption, formula, data, id_var,
                              preferences, task_var, informative, profile_var,
                              design = "uniform") {
  attr_names <- names(pt$attributes)
  point_vec <- .flatten_estimates(pt, "mapns", assumption, attr_names)
  se_vec <- rep(NA_real_, length(point_vec))
  names(se_vec) <- names(point_vec)

  do_sep  <- assumption %in% c("separability", "both")
  do_cond <- assumption %in% c("conditional", "both")

  outcome_var <- all.vars(formula)[1]
  # Only the conditional/heterogeneous SE needs subclassed pairwise info;
  # MAPNS-separability has no closed form here regardless (see below).
  info <- .precompute_pairwise_info(data, outcome_var, attr_names, pt$attributes,
                                    preferences, id_var, task_var, informative,
                                    profile_var, do_sep = FALSE, do_cond = do_cond,
                                    design = design)
  cond_info <- info$cond_info

  for (a in attr_names) {
    amce_a <- pt$amce[[a]]
    levs <- pt$attributes[[a]]
    Dl <- length(levs)

    # SE for AMCEs (these are already known)
    for (lev in names(amce_a$estimate)) {
      nm <- paste0("amce.", a, ".", lev)
      if (nm %in% names(se_vec)) se_vec[nm] <- amce_a$se[lev]
    }

    # Folded normal SE for MAPNS under separability is not supported:
    # MAPNS = max|AMCE| (Proposition 3) and the variance of a max of correlated
    # normals has no closed form. Use se = "parametric" instead.
    if (do_sep) {
      nm <- paste0("mapns.separability.", a)
      if (nm %in% names(se_vec)) se_vec[nm] <- NA_real_
    }

    # For conditional: delta method on the weighted sum, using each group's
    # own direct pairwise SE (from .precompute_pairwise_info) rather than
    # approximating it from base-level AMCE SEs.
    if (do_cond && a %in% names(cond_info)) {
      info_a <- cond_info[[a]]
      is_ranking_a <- info_a$type == "ranking"
      var_mapns_c <- 0

      if (is_ranking_a) {
        # ranking: one contribution per extreme-pair group
        for (ep in names(info_a$groups)) {
          v_ep  <- pt$acmce[[a]][[ep]]
          pi_g  <- v_ep$pi
          sig_g <- .se_or_zero(info_a$groups[[ep]])
          e_g   <- .folded_normal_mean(v_ep$estimate, sig_g)
          var_apns_pair <- pi_g^2 * (v_ep$estimate^2 + sig_g^2 - e_g^2)
          nm_g <- paste0("acmce.", make.names(ep), ".", a)
          if (nm_g %in% names(se_vec)) se_vec[nm_g] <- if (sig_g > 0) sig_g else NA_real_
          var_mapns_c <- var_mapns_c + var_apns_pair
        }
      } else {
        for (pair in names(info_a$groups)) {
          v <- pt$acmce[[a]][[pair]]
          grp_info <- info_a$groups[[pair]]

          if (!is.null(v$pro)) {
            # binary: two groups (pro / con)
            pi_val  <- v$pi
            sig_pro <- .se_or_zero(grp_info[["1"]])
            sig_con <- .se_or_zero(grp_info[["0"]])
            e_pro <- .folded_normal_mean(v$pro, sig_pro)
            e_con <- .folded_normal_mean(v$con, sig_con)
            var_pro <- v$pro^2 + sig_pro^2 - e_pro^2
            var_con <- v$con^2 + sig_con^2 - e_con^2
            var_apns_pair <- pi_val^2 * var_pro + (1 - pi_val)^2 * var_con
            nm_pro <- paste0("acmce.pro.", a, ".", pair)
            nm_con <- paste0("acmce.con.", a, ".", pair)
            if (nm_pro %in% names(se_vec)) se_vec[nm_pro] <- if (sig_pro > 0) sig_pro else NA_real_
            if (nm_con %in% names(se_vec)) se_vec[nm_con] <- if (sig_con > 0) sig_con else NA_real_
          } else {
            # multilevel: K groups
            pi_info <- v$pi
            var_apns_pair <- sum(sapply(names(pi_info), function(g) {
              pi_g   <- pi_info[g]
              amce_g <- v$groups[g]
              sig_g  <- .se_or_zero(grp_info[[g]])
              e_g    <- .folded_normal_mean(amce_g, sig_g)
              pi_g^2 * (amce_g^2 + sig_g^2 - e_g^2)
            }))
            for (g in names(pi_info)) {
              sig_g <- .se_or_zero(grp_info[[g]])
              nm_g <- paste0("acmce.", g, ".", a, ".", pair)
              if (nm_g %in% names(se_vec)) se_vec[nm_g] <- if (sig_g > 0) sig_g else NA_real_
            }
          }
          var_mapns_c <- var_mapns_c + var_apns_pair
        }
      }

      # ranking: MAPNS = sum_g pi_g*|ACMCE_g| (Proposition 8), no averaging
      # over pairs. Non-ranking preferences only identify MAPNS when Dl == 2
      # (the single pair is trivially every group's extreme pair); for Dl > 2
      # MAPNS is not estimable (see .mapns_from_pairwise_cond) and its SE is
      # left unset.
      se_mapns_c <- if (is_ranking_a || Dl == 2) sqrt(var_mapns_c) else NA_real_
      nm <- paste0("mapns.conditional.", a)
      if (nm %in% names(se_vec)) se_vec[nm] <- se_mapns_c
    }
  }

  list(se = se_vec[!is.na(se_vec)], point = point_vec)
}

#' Mean of a folded normal distribution
#' @keywords internal
.folded_normal_mean <- function(mu, sigma) {
  sigma * sqrt(2 / pi) * exp(-mu^2 / (2 * sigma^2)) +
    mu * (1 - 2 * stats::pnorm(-mu / sigma))
}


# ── Jackknife ────────────────────────────────────────────────────────────────

#' @keywords internal
.se_jackknife <- function(formula, data, id, id_var, attr_names,
                          estimand, assumption, preferences,
                          task_var = NULL, informative = "all", profile_var = NULL,
                          design = "uniform") {

  pt <- .estimate_point(formula, data, id, id_var, attr_names,
                        estimand, assumption, preferences, task_var, informative,
                        profile_var, design)
  point_vec <- .flatten_estimates(pt, estimand, assumption, attr_names)

  unique_ids <- unique(data[[id_var]])
  n <- length(unique_ids)

  jack_ests <- vector("list", n)
  for (i in seq_len(n)) {
    d_i <- data[data[[id_var]] != unique_ids[i], , drop = FALSE]
    p_i <- preferences
    if (!is.null(preferences)) {
      p_i <- preferences
      p_i$data <- preferences$data[preferences$data[[preferences$id_var]] != unique_ids[i], ]
    }
    tryCatch({
      est <- .estimate_point(formula, d_i, id, id_var, attr_names,
                             estimand, assumption, p_i, task_var, informative,
                             profile_var, design)
      jack_ests[[i]] <- .flatten_estimates(est, estimand, assumption, attr_names)
    }, error = function(e) { jack_ests[[i]] <<- NULL })
  }

  jack_ests <- Filter(Negate(is.null), jack_ests)
  n_ok <- length(jack_ests)
  mat <- do.call(rbind, lapply(jack_ests, function(x) x[names(point_vec)]))
  theta_bar <- colMeans(mat, na.rm = TRUE)
  se_vec <- sqrt(((n_ok - 1) / n_ok) * colSums(sweep(mat, 2, theta_bar)^2, na.rm = TRUE))
  list(se = se_vec, point = point_vec)
}
