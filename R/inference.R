# ══════════════════════════════════════════════════════════════════════════════
# Standard error methods for cj_apns
# ══════════════════════════════════════════════════════════════════════════════


# ── Parametric bootstrap ─────────────────────────────────────────────────────
#
# The AMCE is asymptotically normal. We draw B samples of AMCE from
# N(amce_hat, se_hat^2), apply the absolute-value / weighted plugin
# transformation, and take the SD of the resulting MAPNS distribution.
# This correctly propagates uncertainty through the |.| mapping.

#' @keywords internal
.se_parametric <- function(formula, data, id, id_var, attr_names,
                           estimand, assumption, preferences, B = 500, alpha = 0.05,
                           task_var = NULL, informative = "all", profile_var = NULL) {

  pt <- .estimate_point(formula, data, id, id_var, attr_names,
                        estimand, assumption, preferences, task_var, informative,
                        profile_var)
  point_vec <- .flatten_estimates(pt, estimand, assumption, attr_names)

  attributes_info <- pt$attributes
  do_sep  <- assumption %in% c("separability", "both")
  do_cond <- assumption %in% c("conditional", "both")

  # Collect AMCE point estimates and SEs
  amce_info <- list()
  for (a in attr_names) {
    amce_a <- pt$amce[[a]]
    for (lev in names(amce_a$estimate))
      amce_info[[paste0(a, ".", lev)]] <- list(
        est = amce_a$estimate[lev], se = amce_a$se[lev])
  }

  # Conditional AMCEs: re-estimate to get SEs per group
  cond_amce_info <- list()
  if (do_cond && !is.null(preferences)) {
    for (a in attr_names) {
      if (!a %in% preferences$attributes) next

      if (preferences$type == "ranking") {
        # ranking: store the direct ACMCE estimate and SE for each extreme-pair
        # group using the same informative-task filter applied to (ep_top, ep_bot)
        # directly — matching the point estimator in .estimate_point.
        pref_data <- preferences$data[
          !duplicated(preferences$data[[preferences$id_var]]), , drop = FALSE]
        pref_ids  <- as.character(pref_data[[preferences$id_var]])
        row_ids   <- as.character(data[[id_var]])
        top_vals  <- as.character(pref_data[[paste0(a, "_top")]])
        bot_vals  <- as.character(pref_data[[paste0(a, "_bot")]])
        ep_label  <- paste0(top_vals, " vs ", bot_vals)
        ep_row    <- ep_label[match(row_ids, pref_ids)]
        outcome_var <- all.vars(formula)[1]

        unique_eps <- unique(ep_label[!is.na(ep_label)])
        for (ep in unique_eps) {
          sub <- data[!is.na(ep_row) & ep_row == ep, , drop = FALSE]
          if (nrow(sub) == 0) next
          parts    <- strsplit(ep, " vs ")[[1]]
          ep_top_i <- parts[1]
          ep_bot_i <- parts[2]

          if (informative == "informative" && !is.null(task_var)) {
            task_key_s     <- paste(sub[[id_var]], sub[[task_var]], sep = ":::")
            profile_vals_s <- if (!is.null(profile_var)) sub[[profile_var]] else NULL
            keep_inf       <- .filter_informative(sub[[a]], task_key_s,
                                                  ep_top_i, ep_bot_i, profile_vals_s)
            sub_inf <- sub[keep_inf, , drop = FALSE]
          } else {
            sub_inf <- sub
          }
          if (nrow(sub_inf) == 0) next

          Y_inf   <- sub_inf[[outcome_var]]
          trt_inf <- as.character(sub_inf[[a]])
          idx_tq  <- which(trt_inf == ep_top_i)
          idx_tp  <- which(trt_inf == ep_bot_i)
          if (length(idx_tq) == 0 || length(idx_tp) == 0) next

          est_direct <- mean(Y_inf[idx_tq], na.rm = TRUE) -
                        mean(Y_inf[idx_tp], na.rm = TRUE)
          se_direct  <- .cluster_se_dim(Y_inf, trt_inf, ep_top_i, ep_bot_i,
                                        sub_inf[[id_var]])

          key <- .make_ep_key(a, "direct", ep)
          cond_amce_info[[key]] <- list(est = est_direct, se = se_direct)
        }
      } else {
        # binary / multilevel: attribute-level groups
        pref_df  <- preferences$data[!duplicated(preferences$data[[preferences$id_var]]),
                                     c(preferences$id_var, a), drop = FALSE]
        pref_map <- stats::setNames(pref_df[[a]], pref_df[[preferences$id_var]])
        data$.pg <- pref_map[as.character(data[[id_var]])]
        has <- !is.na(data$.pg)
        if (sum(has) == 0) { data$.pg <- NULL; next }
        dm <- data[has, ]
        pg_sample   <- dm$.pg[!is.na(dm$.pg)]
        pref_groups <- if (is.character(pg_sample)) sort(unique(pg_sample)) else c(1, 0)
        for (grp in pref_groups) {
          sub <- dm[dm$.pg == grp, ]
          if (nrow(sub) == 0) next
          res <- estimate_amce(formula, sub, id = id,
                               task_var = task_var, informative = informative,
                               profile_var = profile_var)
          for (lev in names(res$amce[[a]]$estimate)) {
            key <- paste0(a, ".", lev, ".grp.", grp)
            cond_amce_info[[key]] <- list(
              est = res$amce[[a]]$estimate[lev], se = res$amce[[a]]$se[lev])
          }
        }
        data$.pg <- NULL
      }
    }
  }

  # Draw B parametric bootstrap replicates
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(point_vec))
  colnames(boot_mat) <- names(point_vec)

  for (b in seq_len(B)) {
    # Resample AMCEs from their asymptotic distribution
    sim_amce <- list()
    for (a in attr_names) {
      levs <- attributes_info[[a]]
      sim_amce[[a]] <- list(base_level = levs[1], levels = levs,
                            estimate = c(), se = c())
      for (lev in levs[-1]) {
        info <- amce_info[[paste0(a, ".", lev)]]
        sim_val <- stats::rnorm(1, info$est, info$se)
        sim_amce[[a]]$estimate[lev] <- sim_val
        sim_amce[[a]]$se[lev] <- info$se
      }
    }

    sim_result <- list(amce = sim_amce, attributes = attributes_info)

    # Compute MAPNS from simulated AMCEs
    sim_mapns <- list()
    sim_acmce <- list()

    for (a in attr_names) {
      levs <- attributes_info[[a]]
      Dl <- length(levs); base <- levs[1]

      if (do_sep) {
        vals <- c()
        for (q in seq_along(levs)) for (p in seq_along(levs)) {
          if (q >= p) next
          vals <- c(vals, abs(get_amce_for_pair(
            sim_result, a, levs[q], levs[p], base)))
        }
        mapns_sep_b <- max(vals)
      }

      mapns_cond_b <- NA_real_
      if (do_cond && a %in% names(pt$acmce) && !is.null(pt$acmce[[a]])) {
        is_ranking_a    <- !is.null(preferences) && preferences$type == "ranking"
        pi_info         <- pt$pi_hat[[a]]
        is_multilevel_a <- !is_ranking_a && length(pi_info) > 1
        vals_c <- c()

        if (is_ranking_a) {
          # ranking: one contribution per extreme-pair group
          for (ep in names(pt$acmce[[a]])) {
            v_ep  <- pt$acmce[[a]][[ep]]
            pi_g  <- v_ep$pi
            v_g   <- .sim_cond_amce_ep(cond_amce_info, a, v_ep$tq, v_ep$tp, base, ep)
            vals_c <- c(vals_c, pi_g * abs(v_g))
            nm_g <- paste0("acmce.", make.names(ep), ".", a)
            if (nm_g %in% names(point_vec)) boot_mat[b, nm_g] <- v_g
          }
          mapns_cond_b <- sum(vals_c)
        } else {
          for (q in seq_along(levs)) for (p in seq_along(levs)) {
            if (q >= p) next
            lev_q <- levs[q]; lev_p <- levs[p]
            pair  <- paste0(lev_q, " vs ", lev_p)
            if (!pair %in% names(pt$acmce[[a]])) next

            if (!is_multilevel_a) {
              # binary: single scalar π for all pairs
              v_pro <- .sim_cond_amce(cond_amce_info, a, lev_q, lev_p, base, 1)
              v_con <- .sim_cond_amce(cond_amce_info, a, lev_q, lev_p, base, 0)
              vals_c <- c(vals_c, pi_info * abs(v_pro) + (1 - pi_info) * abs(v_con))
              nm_pro <- paste0("acmce.pro.", a, ".", pair)
              nm_con <- paste0("acmce.con.", a, ".", pair)
              if (nm_pro %in% names(point_vec)) boot_mat[b, nm_pro] <- v_pro
              if (nm_con %in% names(point_vec)) boot_mat[b, nm_con] <- v_con
            } else {
              # multilevel: K groups, attribute-level π
              grp_amces <- sapply(names(pi_info), function(g)
                .sim_cond_amce(cond_amce_info, a, lev_q, lev_p, base, g))
              vals_c <- c(vals_c, sum(pi_info * abs(grp_amces)))
              for (g in names(pi_info)) {
                nm_g <- paste0("acmce.", g, ".", a, ".", pair)
                if (nm_g %in% names(point_vec)) boot_mat[b, nm_g] <- grp_amces[g]
              }
            }
          }
          mapns_cond_b <- sum(vals_c) / (Dl * (Dl - 1) / 2)
        }
      }

      # Store in boot_mat
      if (do_sep) {
        nm <- paste0("mapns.separability.", a)
        if (nm %in% colnames(boot_mat)) boot_mat[b, nm] <- mapns_sep_b
      }
      if (do_cond && !is.na(mapns_cond_b)) {
        nm <- paste0("mapns.conditional.", a)
        if (nm %in% colnames(boot_mat)) boot_mat[b, nm] <- mapns_cond_b
      }
    }

    # AMCEs
    for (a in attr_names) {
      for (lev in names(sim_amce[[a]]$estimate)) {
        nm <- paste0("amce.", a, ".", lev)
        if (nm %in% colnames(boot_mat))
          boot_mat[b, nm] <- sim_amce[[a]]$estimate[lev]
      }
    }
  }

  se_vec  <- apply(boot_mat, 2, stats::sd, na.rm = TRUE)
  # Algorithm 1 (paper): percentile-based CIs from the empirical bootstrap quantiles
  lower_vec <- apply(boot_mat, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE)
  upper_vec <- apply(boot_mat, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  list(se = se_vec, point = point_vec, lower = lower_vec, upper = upper_vec)
}

#' Simulate a conditional AMCE from the parametric bootstrap (binary/multilevel)
#' @keywords internal
.sim_cond_amce <- function(cond_amce_info, a, tq, tp, base, grp) {
  .get_sim <- function(lev) {
    key <- paste0(a, ".", lev, ".grp.", grp)
    info <- cond_amce_info[[key]]
    if (is.null(info) || is.na(info$se)) return(0)
    stats::rnorm(1, info$est, info$se)
  }
  if (tq == base && tp == base) return(0)
  if (tq == base) return(-.get_sim(tp))
  if (tp == base) return(.get_sim(tq))
  .get_sim(tq) - .get_sim(tp)
}

#' Build the cond_amce_info key for ranking-type extreme-pair groups
#' @keywords internal
.make_ep_key <- function(a, lev, ep) {
  paste0(a, ".", lev, ".ep.", make.names(ep))
}

#' Simulate a conditional AMCE for a ranking-type extreme-pair group
#' @keywords internal
.sim_cond_amce_ep <- function(cond_amce_info, a, tq, tp, base, ep) {
  key  <- .make_ep_key(a, "direct", ep)
  info <- cond_amce_info[[key]]
  if (is.null(info) || is.na(info$se)) return(0)
  stats::rnorm(1, info$est, info$se)
}


# ── Nonparametric bootstrap ──────────────────────────────────────────────────

#' @keywords internal
.se_bootstrap <- function(formula, data, id, id_var, attr_names,
                          estimand, assumption, preferences, B = 500,
                          task_var = NULL, informative = "all", profile_var = NULL) {

  pt <- .estimate_point(formula, data, id, id_var, attr_names,
                        estimand, assumption, preferences, task_var, informative,
                        profile_var)
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
                             profile_var)
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

#' @keywords internal
.se_folded_normal <- function(pt, assumption) {
  attr_names <- names(pt$attributes)
  point_vec <- .flatten_estimates(pt, "mapns", assumption, attr_names)
  se_vec <- rep(NA_real_, length(point_vec))
  names(se_vec) <- names(point_vec)

  do_sep  <- assumption %in% c("separability", "both")
  do_cond <- assumption %in% c("conditional", "both")

  for (a in attr_names) {
    amce_a <- pt$amce[[a]]
    levs <- pt$attributes[[a]]
    Dl <- length(levs); base <- levs[1]

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

    # For conditional: more complex, use delta method on weighted sum
    if (do_cond && !is.null(pt$acmce) && !is.null(pt$acmce[[a]])) {
      # Detect ranking type by checking if acmce entries carry $estimate (vs $pro/$groups)
      is_ranking_a <- length(pt$acmce[[a]]) > 0 &&
                      !is.null(pt$acmce[[a]][[1]]$estimate)
      var_mapns_c <- 0
      for (ep_nm in names(pt$acmce[[a]])) {
        v  <- pt$acmce[[a]][[ep_nm]]
        tq <- strsplit(ep_nm, " vs ")[[1]][1]
        tp <- strsplit(ep_nm, " vs ")[[1]][2]
        if (!is.na(tq) && !is.na(tp)) {
          if (tq == base) sig_approx <- amce_a$se[tp]
          else if (tp == base) sig_approx <- amce_a$se[tq]
          else sig_approx <- sqrt(amce_a$se[tq]^2 + amce_a$se[tp]^2)
          sig_approx <- unname(sig_approx)
          if (is.na(sig_approx)) sig_approx <- 0
        } else {
          sig_approx <- 0
        }

        if (!is.null(v$estimate)) {
          # ranking type: single ACMCE per extreme-pair group
          pi_g  <- v$pi
          sig_g <- if (pi_g > 0) sig_approx / sqrt(pi_g) else 0
          e_g   <- .folded_normal_mean(v$estimate, sig_g)
          var_apns_pair <- pi_g^2 * (v$estimate^2 + sig_g^2 - e_g^2)
          nm_g <- paste0("acmce.", make.names(ep_nm), ".", a)
          if (nm_g %in% names(se_vec)) se_vec[nm_g] <- NA_real_
        } else if (!is.null(v$pro)) {
          # binary: two groups (pro / con)
          pi_val <- v$pi
          sig_pro <- if (pi_val > 0) sig_approx / sqrt(pi_val) else 0
          sig_con <- if ((1 - pi_val) > 0) sig_approx / sqrt(1 - pi_val) else 0
          e_pro <- .folded_normal_mean(v$pro, sig_pro)
          e_con <- .folded_normal_mean(v$con, sig_con)
          var_pro <- v$pro^2 + sig_pro^2 - e_pro^2
          var_con <- v$con^2 + sig_con^2 - e_con^2
          var_apns_pair <- pi_val^2 * var_pro + (1 - pi_val)^2 * var_con
          nm_pro <- paste0("acmce.pro.", a, ".", ep_nm)
          nm_con <- paste0("acmce.con.", a, ".", ep_nm)
          if (nm_pro %in% names(se_vec)) se_vec[nm_pro] <- NA_real_
          if (nm_con %in% names(se_vec)) se_vec[nm_con] <- NA_real_
        } else {
          # multilevel: K groups
          pi_info <- v$pi
          var_apns_pair <- sum(sapply(names(pi_info), function(g) {
            pi_g   <- pi_info[g]
            amce_g <- v$groups[g]
            sig_g  <- if (pi_g > 0) sig_approx / sqrt(pi_g) else 0
            e_g    <- .folded_normal_mean(amce_g, sig_g)
            pi_g^2 * (amce_g^2 + sig_g^2 - e_g^2)
          }))
          for (g in names(pi_info)) {
            nm_g <- paste0("acmce.", g, ".", a, ".", ep_nm)
            if (nm_g %in% names(se_vec)) se_vec[nm_g] <- NA_real_
          }
        }
        var_mapns_c <- var_mapns_c + var_apns_pair
      }
      Dl <- length(levs)
      # ranking: MAPNS = sum_g pi_g*|ACMCE_g|, no averaging over pairs
      se_mapns_c <- if (is_ranking_a) sqrt(var_mapns_c)
                    else sqrt(var_mapns_c) / (Dl * (Dl - 1) / 2)
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
                          task_var = NULL, informative = "all", profile_var = NULL) {

  pt <- .estimate_point(formula, data, id, id_var, attr_names,
                        estimand, assumption, preferences, task_var, informative,
                        profile_var)
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
                             profile_var)
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
