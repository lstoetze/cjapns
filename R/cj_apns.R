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
#'   * `"mapns"` (default): Maximum Average Probability of Necessary and
#'     Sufficient conditions, summarising overall attribute relevance. Under
#'     separability this is the *maximum* |AMCE| across pairwise level
#'     comparisons (Proposition 3/7); under the conditional/heterogeneous
#'     assumption it is the group-share-weighted sum of each preference
#'     group's |ACMCE| at that group's own extreme levels (Proposition 4/8),
#'     which for attributes with more than two levels requires
#'     `preferences` of `type = "ranking"` — see `make_preferences()`.
#'   * `"apns"`: Pairwise Expected Probability of Necessary and Sufficient.
#'   * `"amce"`: Standard Average Marginal Component Effects only.
#' @param assumption Identifying assumption:
#'   * `"separability"` (default): Separable monotonicity. APNS = |AMCE|.
#'   * `"conditional"`: Conditional/heterogeneous separable (transitive)
#'     monotonicity. Requires `preferences`.
#'   * `"both"`: Estimate under both assumptions for comparison.
#' @param preferences An object of class \code{"cj_preferences"} created by
#'   \code{make_preferences}. Required when \code{assumption} is
#'   \code{"conditional"} or \code{"both"}. For MAPNS estimation on
#'   attributes with more than two levels, must be of `type = "ranking"`;
#'   other types only identify MAPNS for two-level attributes (pairwise APNS
#'   is unaffected). See \code{\link{make_preferences}}.
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
#' @param profile A one-sided formula for the profile indicator variable
#'   (e.g., `~ profile`). When provided, informative task detection compares
#'   profile "a" vs profile "b" attribute values explicitly rather than
#'   counting levels within a task.
#' @param design Either \code{"uniform"} (default) or a \code{"cj_design"}
#'   object from \code{\link{make_design}}. Under \code{"uniform"}, every
#'   estimator is the simple pooled difference-in-means (Corollary "Pooled
#'   Estimation under Full Randomization"), unchanged from context-free
#'   estimation. When a design is supplied, `level_probs` reweight (via
#'   inverse-probability/raking weights) the other attributes' realized
#'   levels toward the design's declared marginal probabilities, and
#'   `constraints` exclude tasks whose other-attribute combination is
#'   infeasible — implementing the general context-weighted estimator
#'   (Propositions "Nonparametric Estimation of the APNS/CAPNS"). See
#'   \code{\link{make_design}} and \code{\link{.attach_design_weights}}.
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
#'                data = cnj_cand, id = ~ ResponseId,
#'                tasks = ~ time, profile = ~ profile, se = "none")
#' res
#'
#' # --- MAPNS on all tasks ---
#' res_all <- cj_apns(vote ~ borders + eurobonds + immucard + schools,
#'                    data = cnj_cand, id = ~ ResponseId, tasks = ~ time,
#'                    profile = ~ profile, informative = "all", se = "none")
#'
#' # --- With parametric bootstrap SEs (default) ---
#' res_pb <- cj_apns(vote ~ borders + eurobonds + immucard + schools,
#'                   data = cnj_cand, id = ~ ResponseId,
#'                   tasks = ~ time, profile = ~ profile)
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
#'   tasks = ~ time, profile = ~ profile,
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
                     profile = NULL,
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

  profile_var <- if (!is.null(profile)) all.vars(profile)[1] else NULL
  if (!is.null(profile_var))
    stopifnot("'profile' variable not found in data" = profile_var %in% names(data))

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
                        estimand, assumption, preferences, task_var, informative,
                        profile_var, design)
  .warn_mapns_not_identified(pt, assumption, preferences)

  # ---- standard errors --------------------------------------------------
  se_detail <- NULL; ci <- NULL

  if (se != "none") {
    se_result <- switch(se,
      parametric   = .se_parametric(formula, data, id, id_var, attr_names,
                                     estimand, assumption, preferences, B, alpha,
                                     task_var, informative, profile_var, design),
      bootstrap    = .se_bootstrap(formula, data, id, id_var, attr_names,
                                    estimand, assumption, preferences, B,
                                    task_var, informative, profile_var, design),
      folded_normal = .se_folded_normal(pt, assumption, formula, data, id_var,
                                        preferences, task_var, informative, profile_var,
                                        design),
      jackknife    = .se_jackknife(formula, data, id, id_var, attr_names,
                                    estimand, assumption, preferences,
                                    task_var, informative, profile_var, design)
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
         profile_var = profile_var,
         pref_type = if (!is.null(preferences)) preferences$type else NULL),
    class = "cj_apns"
  )
}



# ══════════════════════════════════════════════════════════════════════════════
# Internal point estimation
# ══════════════════════════════════════════════════════════════════════════════

#' Aggregate pairwise conditional APNS estimates into a MAPNS
#'
#' Proposition 8 identifies MAPNS under heterogeneous separable transitive
#' monotonicity as a share-weighted sum of each preference group's CAPNS at
#' that group's *own* extreme (most-/least-favored) pair, not an average over
#' all pairwise comparisons. For a binary attribute (Dl == 2) there is only
#' one pair, which is trivially every group's extreme pair, so the two
#' formulations coincide. For Dl > 2, computing this requires knowing each
#' group's own extreme pair, which is only available from `type = "ranking"`
#' preferences (see `make_preferences`); `binary`/`scale`/`multilevel`
#' preference types do not carry that information, so MAPNS is not
#' identified and NA is returned (pairwise APNS/ACMCE estimates remain valid
#' and are still reported).
#' @keywords internal
.mapns_from_pairwise_cond <- function(apns_cond, Dl) {
  if (Dl == 2) return(apns_cond[[1]]$estimate)
  NA_real_
}

#' Warn once when conditional MAPNS could not be identified
#'
#' Emits a single consolidated warning listing attributes for which
#' `.mapns_from_pairwise_cond` returned NA (Dl > 2 with non-ranking
#' preferences). Called only from `cj_apns()` on the point estimate, not from
#' inside `.estimate_point` itself, so that bootstrap/jackknife/parametric SE
#' loops (which call `.estimate_point` hundreds of times) do not spam the
#' same warning on every replicate.
#' @keywords internal
.warn_mapns_not_identified <- function(pt, assumption, preferences) {
  if (is.null(preferences) || preferences$type == "ranking") return(invisible(NULL))
  if (assumption == "separability" || is.null(pt$mapns)) return(invisible(NULL))

  bad <- character(0)
  for (a in names(pt$mapns)) {
    val <- pt$mapns[[a]]
    cond_val <- if (assumption == "both" && is.list(val)) val$conditional else val
    if (!is.null(cond_val) && length(cond_val) == 1 && is.na(cond_val) &&
        a %in% names(pt$attributes) && length(pt$attributes[[a]]) > 2) {
      bad <- c(bad, a)
    }
  }

  if (length(bad) > 0)
    warning(
      "MAPNS under conditional/heterogeneous separable monotonicity is not ",
      "identified for attribute(s) ", paste(bad, collapse = ", "),
      " (more than two levels) with preferences$type = \"", preferences$type,
      "\": this requires each preference group's own extreme (most- and ",
      "least-favored) levels, which this preference type does not provide. ",
      "Use type = \"ranking\" (see tidy_ranking_data() / make_preferences()) ",
      "to estimate MAPNS for these attributes. Pairwise APNS/ACMCE estimates ",
      "are still reported.",
      call. = FALSE
    )
  invisible(NULL)
}

#' @keywords internal
.estimate_point <- function(formula, data, id, id_var, attr_names,
                            estimand, assumption, preferences,
                            task_var = NULL, informative = "all",
                            profile_var = NULL, design = "uniform") {

  attributes_info <- lapply(attr_names, function(a) levels(data[[a]]))
  names(attributes_info) <- attr_names
  outcome_var <- all.vars(formula)[1]

  amce_result <- estimate_amce(formula, data, id = id,
                               task_var = task_var, informative = informative,
                               profile_var = profile_var, design = design)

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

    # ── design weighting (Propositions "Nonparametric Estimation of the
    # APNS/CAPNS"): NULL under design = "uniform" reproduces the pooled
    # Corollary exactly; otherwise attach per-task IPW weights for this
    # attribute's context (all other attributes). ────────────────────────
    weights_col <- NULL
    if (!identical(design, "uniform")) {
      data$.design_weight <- .attach_design_weights(data, attr_names, a, design,
                                                     id_var, task_var)
      weights_col <- ".design_weight"
    }

    # ── separable monotonicity ──────────────────────────────────────────
    if (do_sep) {
      apns_sep <- list()
      for (q in seq_along(levs)) for (p in seq_along(levs)) {
        if (q >= p) next
        pair <- paste0(levs[q], " vs ", levs[p])
        # Direct profile-1/2-conditioned estimate for this exact pair
        # (Theorem 1), not reconstructed from base-level coefficients.
        val <- abs(estimate_pairwise(data, outcome_var, a, levs[q], levs[p],
                                     id_var = id_var, task_var = task_var,
                                     informative = informative,
                                     profile_var = profile_var,
                                     weights_col = weights_col)$estimate)
        apns_sep[[pair]] <- list(tq = levs[q], tp = levs[p],
                                 estimate = val, assumption = "separability")
      }
      # Proposition 3: MAPNS = max over unique pairs of APNS = max |AMCE(tq, tp)|
      mapns_sep <- max(sapply(apns_sep, `[[`, "estimate"))
    }

    # ── conditional separable monotonicity ──────────────────────────────
    do_cond_a <- FALSE
    if (do_cond && a %in% preferences$attributes) {
      do_cond_a <- TRUE

      # ── ranking type: look up pre-stored (top, bot) per respondent ─────────
      # preferences$data has <attr>_top and <attr>_bot columns (one row per
      # respondent). Group by extreme pair and compute one ACMCE per group.
      # MAPNS = sum_g pi_g * |ACMCE_g(top_g, bottom_g)|.
      if (preferences$type == "ranking") {
        pref_data <- preferences$data[
          !duplicated(preferences$data[[preferences$id_var]]), , drop = FALSE]
        pref_ids  <- as.character(pref_data[[preferences$id_var]])
        row_ids   <- as.character(data[[id_var]])
        top_vals  <- as.character(pref_data[[paste0(a, "_top")]])
        bot_vals  <- as.character(pref_data[[paste0(a, "_bot")]])
        ep_label  <- paste0(top_vals, " vs ", bot_vals)
        ep_row    <- ep_label[match(row_ids, pref_ids)]

        unique_eps <- unique(ep_label[!is.na(ep_label)])
        apns_cond <- list(); acmce_a <- list()

        for (ep in unique_eps) {
          in_pref  <- !is.na(ep_label) & ep_label == ep
          pi_g     <- mean(in_pref, na.rm = TRUE)
          ep_top   <- top_vals[which(in_pref)[1]]
          ep_bot   <- bot_vals[which(in_pref)[1]]

          data$.pg <- as.integer(!is.na(ep_row) & ep_row == ep)
          dm_g     <- data[data$.pg == 1L, , drop = FALSE]
          if (nrow(dm_g) == 0) next

          # Direct pairwise estimate at this group's own extreme pair
          # (ep_top, ep_bot), restricted to informative tasks (Proposition 8).
          pw_g <- estimate_pairwise(dm_g, outcome_var, a, ep_top, ep_bot,
                                    id_var = id_var, task_var = task_var,
                                    informative = informative,
                                    profile_var = profile_var,
                                    weights_col = weights_col)
          if (is.na(pw_g$estimate)) next
          v_g <- pw_g$estimate

          apns_cond[[ep]] <- list(tq = ep_top, tp = ep_bot,
            estimate = pi_g * abs(v_g), assumption = "conditional")
          acmce_a[[ep]] <- list(estimate = v_g, pi = pi_g,
                                tq = ep_top, tp = ep_bot)
        }
        data$.pg <- NULL

        if (length(apns_cond) == 0) {
          do_cond_a <- FALSE
        } else {
          mapns_cond <- sum(sapply(apns_cond, `[[`, "estimate"))
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
            sub_dm_list <- stats::setNames(
              lapply(pref_groups, function(g) dm[!is.na(dm$.pg) & dm$.pg == g, , drop = FALSE]),
              pref_groups
            )

            apns_cond <- list(); acmce_a <- list()
            for (q in seq_along(levs)) for (p in seq_along(levs)) {
              if (q >= p) next
              pair <- paste0(levs[q], " vs ", levs[p])
              # Direct pairwise estimate for this pair, within each
              # preference-group subclass (generalises Theorem 2 to K groups).
              grp_amces <- sapply(pref_groups, function(g) {
                sub_dm <- sub_dm_list[[g]]
                if (nrow(sub_dm) == 0) return(0)
                pw <- estimate_pairwise(sub_dm, outcome_var, a, levs[q], levs[p],
                                        id_var = id_var, task_var = task_var,
                                        informative = informative,
                                        profile_var = profile_var,
                                        weights_col = weights_col)
                if (is.na(pw$estimate)) 0 else pw$estimate
              })
              names(grp_amces) <- pref_groups
              apns_cond[[pair]] <- list(tq = levs[q], tp = levs[p],
                estimate = sum(pi_vals * abs(grp_amces)), assumption = "conditional")
              acmce_a[[pair]] <- list(groups = grp_amces, pi = pi_vals)
            }
            mapns_cond <- .mapns_from_pairwise_cond(apns_cond, Dl)
            pi_hat[[a]] <- pi_vals; acmce[[a]] <- acmce_a
          } else {
            pi_val <- mean(dm$.pg, na.rm = TRUE)
            dm_pro <- dm[dm$.pg == 1, , drop = FALSE]
            dm_con <- dm[dm$.pg == 0, , drop = FALSE]

            apns_cond <- list(); acmce_a <- list()
            for (q in seq_along(levs)) for (p in seq_along(levs)) {
              if (q >= p) next
              pair <- paste0(levs[q], " vs ", levs[p])
              # Direct pairwise estimate for this pair, within the pro/con
              # subclass (Proposition 2/8), not reconstructed from base-level
              # regression coefficients.
              v_pro <- estimate_pairwise(dm_pro, outcome_var, a, levs[q], levs[p],
                                         id_var = id_var, task_var = task_var,
                                         informative = informative,
                                         profile_var = profile_var,
                                         weights_col = weights_col)$estimate
              v_con <- estimate_pairwise(dm_con, outcome_var, a, levs[q], levs[p],
                                         id_var = id_var, task_var = task_var,
                                         informative = informative,
                                         profile_var = profile_var,
                                         weights_col = weights_col)$estimate
              if (is.na(v_pro)) v_pro <- 0
              if (is.na(v_con)) v_con <- 0
              apns_cond[[pair]] <- list(tq = levs[q], tp = levs[p],
                estimate = pi_val * abs(v_pro) + (1 - pi_val) * abs(v_con),
                assumption = "conditional")
              acmce_a[[pair]] <- list(pro = v_pro, con = v_con, pi = pi_val)
            }
            mapns_cond <- .mapns_from_pairwise_cond(apns_cond, Dl)
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
      for (ep in names(pt$acmce[[a]])) {
        v <- pt$acmce[[a]][[ep]]
        if (!is.null(v$estimate)) {
          # ranking type: single ACMCE for this extreme-pair group
          out <- c(out, stats::setNames(v$estimate,
                   paste0("acmce.", make.names(ep), ".", a)))
        } else if (!is.null(v$pro)) {
          out <- c(out, stats::setNames(v$pro, paste0("acmce.pro.", a, ".", ep)))
          out <- c(out, stats::setNames(v$con, paste0("acmce.con.", a, ".", ep)))
        } else {
          for (g in names(v$groups))
            out <- c(out, stats::setNames(v$groups[[g]], paste0("acmce.", g, ".", a, ".", ep)))
        }
      }
    }
  }
  out
}
