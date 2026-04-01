# ══════════════════════════════════════════════════════════════════════════════
# Standard error methods for cj_apns
# ══════════════════════════════════════════════════════════════════════════════


# ── Parametric bootstrap ─────────────────────────────────────────────────────
#
# The AMCE is asymptotically normal. We draw B samples of AMCE from
# N(amce_hat, se_hat^2), apply the absolute-value / weighted plugin
# transformation, and take the SD of the resulting EAPNS distribution.
# This correctly propagates uncertainty through the |.| mapping.

#' @keywords internal
.se_parametric <- function(formula, data, id, id_var, attr_names,
                           estimand, assumption, preferences, B = 500) {

  pt <- .estimate_point(formula, data, id, id_var, attr_names,
                        estimand, assumption, preferences)
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
        # ranking: pair-specific binary groups → keys include pair identity
        rank_mat <- preferences$rank_data[[a]]
        rank_ids <- rank_mat[[preferences$id_var]]
        row_ids  <- as.character(data[[id_var]])
        idx      <- match(row_ids, rank_ids)
        levs_a   <- attributes_info[[a]]

        for (q in seq_along(levs_a)) for (p in seq_along(levs_a)) {
          if (q >= p) next
          tq <- levs_a[q]; tp <- levs_a[p]
          if (!tq %in% names(rank_mat) || !tp %in% names(rank_mat)) next
          pg <- as.integer(rank_mat[[tq]][idx] < rank_mat[[tp]][idx])
          for (grp in c(1L, 0L)) {
            sub <- data[!is.na(pg) & pg == grp, , drop = FALSE]
            if (nrow(sub) == 0) next
            res <- estimate_amce(formula, sub, id = id)
            if (!a %in% names(res$amce)) next
            for (lev in names(res$amce[[a]]$estimate)) {
              key <- .make_rank_key(a, lev, tq, tp, grp)
              cond_amce_info[[key]] <- list(
                est = res$amce[[a]]$estimate[lev],
                se  = res$amce[[a]]$se[lev])
            }
          }
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
          res <- estimate_amce(formula, sub, id = id)
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

    # Compute EAPNS from simulated AMCEs
    sim_eapns <- list()
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
        eapns_sep_b <- sum(vals) / (Dl - 1)
      }

      eapns_cond_b <- NA_real_
      if (do_cond && a %in% names(pt$acmce) && !is.null(pt$acmce[[a]])) {
        is_ranking_a  <- !is.null(preferences) && preferences$type == "ranking"
        pi_info       <- pt$pi_hat[[a]]
        is_multilevel_a <- !is_ranking_a && length(pi_info) > 1
        vals_c <- c()
        for (q in seq_along(levs)) for (p in seq_along(levs)) {
          if (q >= p) next
          lev_q <- levs[q]; lev_p <- levs[p]
          pair  <- paste0(lev_q, " vs ", lev_p)
          if (!pair %in% names(pt$acmce[[a]])) next

          if (is_ranking_a) {
            # ranking: pair-specific binary groups
            pi_pair <- pt$acmce[[a]][[pair]]$pi
            v_pro <- .sim_cond_amce_rank(cond_amce_info, a, lev_q, lev_p, base, lev_q, lev_p, 1L)
            v_con <- .sim_cond_amce_rank(cond_amce_info, a, lev_q, lev_p, base, lev_q, lev_p, 0L)
            vals_c <- c(vals_c, pi_pair * abs(v_pro) + (1 - pi_pair) * abs(v_con))
            nm_pro <- paste0("acmce.pro.", a, ".", pair)
            nm_con <- paste0("acmce.con.", a, ".", pair)
            if (nm_pro %in% names(point_vec)) boot_mat[b, nm_pro] <- v_pro
            if (nm_con %in% names(point_vec)) boot_mat[b, nm_con] <- v_con
          } else if (!is_multilevel_a) {
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
        eapns_cond_b <- sum(vals_c) / (Dl - 1)
      }

      # Store in boot_mat
      if (do_sep) {
        nm <- paste0("eapns.separability.", a)
        if (nm %in% colnames(boot_mat)) boot_mat[b, nm] <- eapns_sep_b
      }
      if (do_cond && !is.na(eapns_cond_b)) {
        nm <- paste0("eapns.conditional.", a)
        if (nm %in% colnames(boot_mat)) boot_mat[b, nm] <- eapns_cond_b
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

  se_vec <- apply(boot_mat, 2, stats::sd, na.rm = TRUE)
  list(se = se_vec, point = point_vec)
}

#' Simulate a conditional AMCE from the parametric bootstrap (binary/multilevel)
#' @keywords internal
.sim_cond_amce <- function(cond_amce_info, a, tq, tp, base, grp) {
  .get_sim <- function(lev) {
    key <- paste0(a, ".", lev, ".grp.", grp)
    info <- cond_amce_info[[key]]
    if (is.null(info)) return(0)
    stats::rnorm(1, info$est, info$se)
  }
  if (tq == base && tp == base) return(0)
  if (tq == base) return(-.get_sim(tp))
  if (tp == base) return(.get_sim(tq))
  .get_sim(tq) - .get_sim(tp)
}

#' Build the cond_amce_info key for ranking-type preferences
#' @keywords internal
.make_rank_key <- function(a, lev, tq, tp, grp) {
  paste0(a, ".", lev, ".pair.", make.names(tq), ".", make.names(tp), ".grp.", grp)
}

#' Simulate a conditional AMCE for ranking-type preferences (pair-specific keys)
#' @keywords internal
.sim_cond_amce_rank <- function(cond_amce_info, a, tq, tp, base, pair_tq, pair_tp, grp) {
  .get_sim <- function(lev) {
    key  <- .make_rank_key(a, lev, pair_tq, pair_tp, grp)
    info <- cond_amce_info[[key]]
    if (is.null(info)) return(0)
    stats::rnorm(1, info$est, info$se)
  }
  if (tq == base && tp == base) return(0)
  if (tq == base) return(-.get_sim(tp))
  if (tp == base) return(.get_sim(tq))
  .get_sim(tq) - .get_sim(tp)
}


# ── Nonparametric bootstrap ──────────────────────────────────────────────────

#' @keywords internal
.se_bootstrap <- function(formula, data, id, id_var, attr_names,
                          estimand, assumption, preferences, B = 500) {

  pt <- .estimate_point(formula, data, id, id_var, attr_names,
                        estimand, assumption, preferences)
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
                             estimand, assumption, bp)
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
  point_vec <- .flatten_estimates(pt, "eapns", assumption, attr_names)
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

    # Folded normal SE for EAPNS under separability
    if (do_sep) {
      var_eapns <- 0
      n_pairs <- 0
      for (q in seq_along(levs)) for (p in seq_along(levs)) {
        if (q >= p) next
        n_pairs <- n_pairs + 1
        mu <- get_amce_for_pair(list(amce = pt$amce), a, levs[q], levs[p], base)
        # SE for this pair: approximate from base-level SEs
        if (levs[q] == base) {
          sig <- amce_a$se[levs[p]]
        } else if (levs[p] == base) {
          sig <- amce_a$se[levs[q]]
        } else {
          sig <- sqrt(amce_a$se[levs[q]]^2 + amce_a$se[levs[p]]^2)
        }
        e_abs <- .folded_normal_mean(mu, sig)
        var_abs <- mu^2 + sig^2 - e_abs^2
        var_eapns <- var_eapns + var_abs
      }
      # EAPNS = sum / (Dl - 1), so Var(EAPNS) = sum(Var) / (Dl-1)^2
      se_eapns <- sqrt(var_eapns) / (Dl - 1)
      nm <- paste0("eapns.separability.", a)
      if (nm %in% names(se_vec)) se_vec[nm] <- se_eapns
    }

    # For conditional: more complex, use delta method on weighted sum
    if (do_cond && !is.null(pt$acmce) && !is.null(pt$acmce[[a]])) {
      var_eapns_c <- 0
      for (pair_nm in names(pt$acmce[[a]])) {
        v  <- pt$acmce[[a]][[pair_nm]]
        tq <- strsplit(pair_nm, " vs ")[[1]][1]
        tp <- strsplit(pair_nm, " vs ")[[1]][2]
        if (tq == base) sig_approx <- amce_a$se[tp]
        else if (tp == base) sig_approx <- amce_a$se[tq]
        else sig_approx <- sqrt(amce_a$se[tq]^2 + amce_a$se[tp]^2)

        if (!is.null(v$pro)) {
          # binary or ranking: v$pi is the pair-specific (or attribute-level) π
          pi_val <- v$pi
          sig_pro <- sig_approx / sqrt(pi_val)
          sig_con <- sig_approx / sqrt(1 - pi_val)
          e_pro <- .folded_normal_mean(v$pro, sig_pro)
          e_con <- .folded_normal_mean(v$con, sig_con)
          var_pro <- v$pro^2 + sig_pro^2 - e_pro^2
          var_con <- v$con^2 + sig_con^2 - e_con^2
          var_epns_pair <- pi_val^2 * var_pro + (1 - pi_val)^2 * var_con
          nm_pro <- paste0("acmce.pro.", a, ".", pair_nm)
          nm_con <- paste0("acmce.con.", a, ".", pair_nm)
          if (nm_pro %in% names(se_vec)) se_vec[nm_pro] <- NA_real_
          if (nm_con %in% names(se_vec)) se_vec[nm_con] <- NA_real_
        } else {
          # multilevel: K groups
          pi_info <- v$pi
          var_epns_pair <- sum(sapply(names(pi_info), function(g) {
            pi_g   <- pi_info[g]
            amce_g <- v$groups[g]
            sig_g  <- sig_approx / sqrt(pi_g)
            e_g    <- .folded_normal_mean(amce_g, sig_g)
            pi_g^2 * (amce_g^2 + sig_g^2 - e_g^2)
          }))
          for (g in names(pi_info)) {
            nm_g <- paste0("acmce.", g, ".", a, ".", pair_nm)
            if (nm_g %in% names(se_vec)) se_vec[nm_g] <- NA_real_
          }
        }
        var_eapns_c <- var_eapns_c + var_epns_pair
      }
      Dl <- length(levs)
      se_eapns_c <- sqrt(var_eapns_c) / (Dl - 1)
      nm <- paste0("eapns.conditional.", a)
      if (nm %in% names(se_vec)) se_vec[nm] <- se_eapns_c
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
                          estimand, assumption, preferences) {

  pt <- .estimate_point(formula, data, id, id_var, attr_names,
                        estimand, assumption, preferences)
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
                             estimand, assumption, p_i)
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
