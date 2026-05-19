#' Estimate Causal Attribution in Conjoint Experiments
#'
#' Main function for estimating the Average Probability of Necessary and
#' Sufficient conditions (APNS) in conjoint experiments. Computes AMCEs,
#' APNS, and/or MAPNS under separable or conditional separable monotonicity,
#' with multiple options for inference.
#'
#' @param formula A formula: `outcome ~ attr1 + attr2 + ...`. The LHS is the
#'   binary choice outcome; the RHS lists conjoint attributes (as factors or
#'   0/1 integers).
#' @param data A data.frame in long format (one row per profile). Must
#'   contain the outcome, attribute columns, and a respondent identifier.
#' @param id A one-sided formula for the respondent ID variable
#'   (e.g., `~ ResponseId`). Required.
#' @param estimand Which estimand to compute:
#'   * `"mapns"` (default): Expected Average Probability of Necessary and
#'     Sufficient conditions, summarising overall attribute relevance.
#'   * `"apns"`: Pairwise Expected Probability of Necessary and Sufficient.
#'   * `"amce"`: Standard Average Marginal Component Effects only.
#' @param assumption Identifying assumption:
#'   * `"separability"` (default): Separable monotonicity. APNS = |AMCE|.
#'   * `"conditional"`: Conditional separable monotonicity. Requires
#'     `preferences`.
#'   * `"both"`: Estimate under both assumptions for comparison.
#' @param preferences An object of class \code{"cj_preferences"} created by
#'   \code{make_preferences}. Required when \code{assumption} is
#'   \code{"conditional"} or \code{"both"}.
#' @param se Method for standard error estimation:
#'   * `"none"`: No standard errors (fastest).
#'   * `"parametric"` (default): Parametric bootstrap — resamples AMCEs
#'     from their asymptotic normal distribution and applies the plugin
#'     transformation. Fast and appropriate for the absolute-value mapping.
#'   * `"bootstrap"`: Nonparametric bootstrap — resamples respondents with
#'     replacement and re-estimates the full pipeline.
#'   * `"folded_normal"`: Analytical SEs using the folded normal
#'     approximation for |AMCE|.
#'   * `"jackknife"`: Leave-one-respondent-out jackknife.
#' @param B Number of replications for `se = "parametric"` or
#'   `se = "bootstrap"`. Default 500.
#' @param alpha Significance level for confidence intervals. Default 0.05.
#' @param tasks A one-sided formula for the task-number variable
#'   (e.g., `~ time`). Required when `informative = "informative"`.
#' @param informative Whether to restrict AMCE estimation to informative tasks
#'   only (`"informative"`, default) or use all tasks (`"all"`). Informative
#'   tasks are those where exactly one profile shows level tq and all J-1
#'   remaining profiles show tp (or vice versa), per Definition 3 in the paper.
#'   Requires `tasks` to be specified; falls back to `"all"` with a warning if
#'   `tasks` is not provided.
#' @param design Either \code{"uniform"} (default) or an object from
#'   \code{make_design}.
#'
#' @return An object of class `"cj_apns"` with components:
#'   * `estimand`, `assumption`, `se_method`: as requested.
#'   * `amce`: AMCE estimates per attribute.
#'   * `apns`, `mapns`: causal attribution estimates (if requested).
#'   * `acmce`, `pi_hat`: conditional AMCEs and group shares (if conditional).
#'   * `ci`: confidence intervals (if SEs requested).
#'   * `se_detail`: named vector of standard errors.
#'
#' @references
#' Stoetzer, L.F. and Magazinnik, A. (2026). Measuring Attribute Relevance
#' in Conjoint Analysis.
#'
#' Hainmueller, J., Hopkins, D.J. and Yamamoto, T. (2014). Causal Inference
#' in Conjoint Analysis. *Political Analysis*, 22(1), 1--30.
#'
#' @examples
#' \dontrun{
#' data(cnj_cand)
#'
#' # --- MAPNS on informative tasks (default), no SEs ---
#' res <- cj_apns(vote ~ borders + eurobonds + immucard + schools,
#'                data = cnj_cand, id = ~ ResponseId, tasks = ~ time,
#'                se = "none")
#' res
#'
#' # --- MAPNS on all tasks ---
#' res_all <- cj_apns(vote ~ borders + eurobonds + immucard + schools,
#'                    data = cnj_cand, id = ~ ResponseId, tasks = ~ time,
#'                    informative = "all", se = "none")
#'
#' # --- With parametric bootstrap SEs (default) ---
#' res_pb <- cj_apns(vote ~ borders + eurobonds + immucard + schools,
#'                   data = cnj_cand, id = ~ ResponseId, tasks = ~ time)
#' summary(res_pb)
#'
#' # --- Conditional separable monotonicity ---
#' pref_dat <- cnj_cand[!duplicated(cnj_cand$ResponseId),
#'                       c("ResponseId", "Att_borders", "Att_eurobonds",
#'                         "Att_immucard", "Att_schools", "Att_tracingapp")]
#' names(pref_dat) <- gsub("^Att_", "", names(pref_dat))
#' prefs <- make_preferences(pref_dat, id = ~ ResponseId, type = "binary")
#'
#' res_cond <- cj_apns(
#'   vote ~ borders + eurobonds + immucard + schools,
#'   data = cnj_cand, id = ~ ResponseId,
#'   assumption = "both", preferences = prefs,
#'   se = "parametric", B = 500
#' )
#' summary(res_cond)
#' plot(res_cond)
#' }
#'
#' @export
cj_apns <- function(formula, data, id,
                     estimand = c("mapns", "apns", "amce"),
                     assumption = c("separability", "conditional", "both"),
                     preferences = NULL,
                     se = c("parametric", "none", "bootstrap",
                            "folded_normal", "jackknife"),
                     B = 500, alpha = 0.05,
                     tasks = NULL,
                     informative = c("informative", "all"),
                     design = "uniform") {

  cl <- match.call()
  estimand <- match.arg(estimand)
  assumption <- match.arg(assumption)
  se <- match.arg(se)
  informative <- match.arg(informative)

  # ---- validation -------------------------------------------------------
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

  if (assumption %in% c("conditional", "both") && is.null(preferences))
    stop("'preferences' must be provided for assumption = \"", assumption,
         "\". Use make_preferences().")

  tt <- stats::terms(formula, data = data)
  attr_names <- attr(tt, "term.labels")
  attr_names <- attr_names[!grepl(":", attr_names)]
  for (a in attr_names)
    if (!is.factor(data[[a]])) data[[a]] <- as.factor(data[[a]])

  # ---- point estimates --------------------------------------------------
  pt <- .estimate_point(formula, data, id, id_var, attr_names,
                        estimand, assumption, preferences, task_var, informative)

  # ---- standard errors --------------------------------------------------
  se_detail <- NULL; ci <- NULL

  if (se != "none") {
    se_result <- switch(se,
      parametric   = .se_parametric(formula, data, id, id_var, attr_names,
                                     estimand, assumption, preferences, B, alpha,
                                     task_var, informative),
      bootstrap    = .se_bootstrap(formula, data, id, id_var, attr_names,
                                    estimand, assumption, preferences, B,
                                    task_var, informative),
      folded_normal = .se_folded_normal(pt, assumption),
      jackknife    = .se_jackknife(formula, data, id, id_var, attr_names,
                                    estimand, assumption, preferences,
                                    task_var, informative)
    )
    se_detail <- se_result$se
    # Parametric bootstrap uses percentile CIs (Algorithm 1 in paper); all
    # other methods fall back to normal-based CIs from the estimated SE.
    if (se == "parametric" && !is.null(se_result$lower)) {
      ci <- data.frame(
        parameter = names(se_detail),
        estimate = se_result$point[names(se_detail)],
        se = se_detail,
        lower = se_result$lower[names(se_detail)],
        upper = se_result$upper[names(se_detail)],
        row.names = NULL, stringsAsFactors = FALSE
      )
    } else {
      ci <- data.frame(
        parameter = names(se_detail),
        estimate = se_result$point[names(se_detail)],
        se = se_detail,
        lower = se_result$point[names(se_detail)] - stats::qnorm(1 - alpha/2) * se_detail,
        upper = se_result$point[names(se_detail)] + stats::qnorm(1 - alpha/2) * se_detail,
        row.names = NULL, stringsAsFactors = FALSE
      )
    }
  }

  structure(
    list(estimand = estimand, assumption = assumption, se_method = se,
         attributes = pt$attributes, amce = pt$amce,
         apns = pt$apns, mapns = pt$mapns,
         pi_hat = pt$pi_hat, acmce = pt$acmce,
         ci = ci, se_detail = se_detail,
         alpha = alpha, B = B, call = cl,
         task_var = task_var, informative = informative,
         pref_type = if (!is.null(preferences)) preferences$type else NULL),
    class = "cj_apns"
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal point estimation
# ══════════════════════════════════════════════════════════════════════════════

#' @keywords internal
.estimate_point <- function(formula, data, id, id_var, attr_names,
                            estimand, assumption, preferences,
                            task_var = NULL, informative = "all") {

  attributes_info <- lapply(attr_names, function(a) levels(data[[a]]))
  names(attributes_info) <- attr_names

  amce_result <- estimate_amce(formula, data, id = id,
                               task_var = task_var, informative = informative)

  if (estimand == "amce")
    return(list(attributes = attributes_info, amce = amce_result$amce,
                apns = NULL, mapns = NULL, pi_hat = NULL, acmce = NULL))

  apns <- list(); mapns <- list()
  pi_hat <- list(); acmce <- list()
  do_sep  <- assumption %in% c("separability", "both")
  do_cond <- assumption %in% c("conditional", "both")

  for (a in attr_names) {
    levs <- attributes_info[[a]]
    Dl <- length(levs); base <- levs[1]

    # ── separable monotonicity ──────────────────────────────────────────
    if (do_sep) {
      apns_sep <- list()
      for (q in seq_along(levs)) for (p in seq_along(levs)) {
        if (q >= p) next
        pair <- paste0(levs[q], " vs ", levs[p])
        val <- abs(get_amce_for_pair(amce_result, a, levs[q], levs[p], base))
        apns_sep[[pair]] <- list(tq = levs[q], tp = levs[p],
                                 estimate = val, assumption = "separability")
      }
      # Paper Eq. (10): MAPNS = [2/(Dl*(Dl-1))] * sum_{unique pairs} APNS = avg over unique pairs
      mapns_sep <- sum(sapply(apns_sep, `[[`, "estimate")) / (Dl * (Dl - 1) / 2)
    }

    # ── conditional separable monotonicity ──────────────────────────────
    do_cond_a <- FALSE
    if (do_cond && a %in% preferences$attributes) {
      do_cond_a <- TRUE

      # ── ranking type: pair-specific p^{qp}_{il} (Proposition 8) ─────
      if (preferences$type == "ranking") {
        rank_mat  <- preferences$rank_data[[a]]
        rank_ids  <- rank_mat[[preferences$id_var]]
        row_ids   <- as.character(data[[id_var]])
        idx       <- match(row_ids, rank_ids)

        apns_cond <- list(); acmce_a <- list()
        for (q in seq_along(levs)) for (p in seq_along(levs)) {
          if (q >= p) next
          tq <- levs[q]; tp <- levs[p]
          pair <- paste0(tq, " vs ", tp)
          if (!tq %in% names(rank_mat) || !tp %in% names(rank_mat)) next

          # π computed at respondent level (avoids row-duplication bias)
          pi_val <- mean(rank_mat[[tq]] < rank_mat[[tp]], na.rm = TRUE)

          # pair-specific binary indicator mapped to conjoint rows
          data$.pg <- as.integer(rank_mat[[tq]][idx] < rank_mat[[tp]][idx])
          has <- !is.na(data$.pg)
          if (sum(has) == 0) next

          dm       <- data[has, , drop = FALSE]
          has_pro  <- any(dm$.pg == 1L)
          has_con  <- any(dm$.pg == 0L)
          # Need at least one group to estimate anything
          if (!has_pro && !has_con) next

          v_pro <- if (has_pro) get_amce_for_pair(
            estimate_amce(formula, dm[dm$.pg == 1L, ], id = id,
                          task_var = task_var, informative = informative),
            a, tq, tp, base) else 0
          v_con <- if (has_con) get_amce_for_pair(
            estimate_amce(formula, dm[dm$.pg == 0L, ], id = id,
                          task_var = task_var, informative = informative),
            a, tq, tp, base) else 0

          apns_cond[[pair]] <- list(tq = tq, tp = tp,
            estimate = pi_val * abs(v_pro) + (1 - pi_val) * abs(v_con),
            assumption = "conditional")
          acmce_a[[pair]] <- list(pro = v_pro, con = v_con, pi = pi_val)
        }
        data$.pg <- NULL

        if (length(apns_cond) == 0) {
          do_cond_a <- FALSE
        } else {
          mapns_cond <- sum(sapply(apns_cond, `[[`, "estimate")) / (Dl * (Dl - 1) / 2)
          # store pair-specific π as a named vector
          pi_hat[[a]] <- stats::setNames(sapply(acmce_a, `[[`, "pi"), names(acmce_a))
          acmce[[a]]  <- acmce_a
        }

      } else {
        # ── binary / multilevel: attribute-level grouping ──────────────
        pref_df  <- preferences$data[!duplicated(preferences$data[[preferences$id_var]]),
                                     c(preferences$id_var, a), drop = FALSE]
        pref_map <- stats::setNames(pref_df[[a]], pref_df[[preferences$id_var]])
        data$.pg <- pref_map[as.character(data[[id_var]])]
        has      <- !is.na(data$.pg)

        if (sum(has) == 0) { do_cond_a <- FALSE
        } else {
          dm       <- data[has, , drop = FALSE]
          pg_vals  <- dm$.pg[!is.na(dm$.pg)]
          is_multilevel <- is.character(pg_vals)

          if (is_multilevel) {
            pref_groups <- sort(unique(pg_vals))
            pi_vals <- stats::setNames(
              sapply(pref_groups, function(g) mean(dm$.pg == g, na.rm = TRUE)),
              pref_groups
            )
            amce_list <- stats::setNames(
              lapply(pref_groups, function(g) {
                sub_dm <- dm[!is.na(dm$.pg) & dm$.pg == g, , drop = FALSE]
                if (nrow(sub_dm) == 0) return(NULL)
                estimate_amce(formula, sub_dm, id = id,
                              task_var = task_var, informative = informative)
              }),
              pref_groups
            )

            apns_cond <- list(); acmce_a <- list()
            for (q in seq_along(levs)) for (p in seq_along(levs)) {
              if (q >= p) next
              pair <- paste0(levs[q], " vs ", levs[p])
              grp_amces <- sapply(pref_groups, function(g) {
                if (is.null(amce_list[[g]])) return(0)
                get_amce_for_pair(amce_list[[g]], a, levs[q], levs[p], base)
              })
              names(grp_amces) <- pref_groups
              apns_cond[[pair]] <- list(tq = levs[q], tp = levs[p],
                estimate = sum(pi_vals * abs(grp_amces)), assumption = "conditional")
              acmce_a[[pair]] <- list(groups = grp_amces, pi = pi_vals)
            }
            mapns_cond <- sum(sapply(apns_cond, `[[`, "estimate")) / (Dl * (Dl - 1) / 2)
            pi_hat[[a]] <- pi_vals; acmce[[a]] <- acmce_a
          } else {
            pi_val   <- mean(dm$.pg, na.rm = TRUE)
            amce_pro <- estimate_amce(formula, dm[dm$.pg == 1, ], id = id,
                                      task_var = task_var, informative = informative)
            amce_con <- estimate_amce(formula, dm[dm$.pg == 0, ], id = id,
                                      task_var = task_var, informative = informative)

            apns_cond <- list(); acmce_a <- list()
            for (q in seq_along(levs)) for (p in seq_along(levs)) {
              if (q >= p) next
              pair <- paste0(levs[q], " vs ", levs[p])
              v_pro <- get_amce_for_pair(amce_pro, a, levs[q], levs[p], base)
              v_con <- get_amce_for_pair(amce_con, a, levs[q], levs[p], base)
              apns_cond[[pair]] <- list(tq = levs[q], tp = levs[p],
                estimate = pi_val * abs(v_pro) + (1 - pi_val) * abs(v_con),
                assumption = "conditional")
              acmce_a[[pair]] <- list(pro = v_pro, con = v_con, pi = pi_val)
            }
            mapns_cond <- sum(sapply(apns_cond, `[[`, "estimate")) / (Dl * (Dl - 1) / 2)
            pi_hat[[a]] <- pi_val; acmce[[a]] <- acmce_a
          }
        }
        data$.pg <- NULL
      }
    }

    # ── store ───────────────────────────────────────────────────────────
    if (assumption == "separability") {
      apns[[a]] <- apns_sep; mapns[[a]] <- mapns_sep
    } else if (assumption == "conditional" && do_cond_a) {
      apns[[a]] <- apns_cond; mapns[[a]] <- mapns_cond
    } else if (assumption == "both") {
      apns[[a]]  <- list(separability = apns_sep,
                         conditional = if (do_cond_a) apns_cond else NULL)
      mapns[[a]] <- list(separability = mapns_sep,
                         conditional = if (do_cond_a) mapns_cond else NULL)
    }
  }

  list(attributes = attributes_info, amce = amce_result$amce,
       apns = if (estimand %in% c("apns", "mapns")) apns else NULL,
       mapns = if (estimand == "mapns") mapns else NULL,
       pi_hat = if (do_cond) pi_hat else NULL,
       acmce = if (do_cond) acmce else NULL)
}


#' Flatten point estimates into a named vector
#' @keywords internal
.flatten_estimates <- function(pt, estimand, assumption, attr_names) {
  out <- c()
  for (a in attr_names) {
    est <- pt$amce[[a]]$estimate
    names(est) <- paste0("amce.", a, ".", names(est))
    out <- c(out, est)
  }
  if (!is.null(pt$mapns)) {
    for (a in attr_names) {
      val <- pt$mapns[[a]]
      if (assumption == "both" && is.list(val)) {
        for (asn in c("separability", "conditional"))
          if (!is.null(val[[asn]]))
            out <- c(out, stats::setNames(val[[asn]], paste0("mapns.", asn, ".", a)))
      } else if (is.numeric(val))
        out <- c(out, stats::setNames(val, paste0("mapns.", assumption, ".", a)))
    }
  }
  if (!is.null(pt$acmce)) {
    for (a in attr_names) if (!is.null(pt$acmce[[a]])) {
      for (pair in names(pt$acmce[[a]])) {
        v <- pt$acmce[[a]][[pair]]
        if (!is.null(v$pro)) {
          out <- c(out, stats::setNames(v$pro, paste0("acmce.pro.", a, ".", pair)))
          out <- c(out, stats::setNames(v$con, paste0("acmce.con.", a, ".", pair)))
        } else {
          for (g in names(v$groups))
            out <- c(out, stats::setNames(v$groups[[g]], paste0("acmce.", g, ".", a, ".", pair)))
        }
      }
    }
  }
  out
}
