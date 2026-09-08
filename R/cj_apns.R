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
#'     Preference-ordering groups with zero informative-task observations
#'     are excluded and the remaining groups' weights renormalized (see
#'     `mapns_coverage`, `missing_groups`, `capns_table` in the return value).
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
#' @param split_sample Logical, default `FALSE`. Under separable monotonicity,
#'   the separability MAPNS is the *maximum* of several noisy pairwise APNS
#'   estimates (Proposition 3/7); because `max(.)` is convex, this plug-in
#'   max is upward-biased for the true max (Jensen's inequality), especially
#'   when several pairs are near-tied. When `TRUE`, `mapns[[attribute]]`'s
#'   separability entry is replaced with a split-sample estimate for every
#'   attribute: across `n_splits` random respondent-level halves, one half
#'   selects the argmax pair and the other estimates that (already-fixed)
#'   pair's APNS, and the median across splits is reported in place of the
#'   ordinary plug-in max. This only affects `estimand = "mapns"` under
#'   `assumption` `"separability"` or `"both"`; it is a no-op (with a
#'   warning) otherwise. It is a bias-reduction device, not an unbiased
#'   estimator of the true max — see \code{\link{mapns_split_test}}, which
#'   reports the full plug-in-vs-split-sample comparison for diagnostic use,
#'   and whose method this implements internally
#'   (\code{\link{.split_sample_mapns}}). The split-sample estimate has its
#'   own SE/CI, reported separately in `mapns_split` (a Rubin's-rules
#'   combination of within- and between-split variance — see
#'   \code{\link{.split_sample_mapns}} for the formula and its
#'   approximations). The main `se`/`ci` output (`se != "none"`) is still
#'   computed around the ordinary plug-in max, not the split-sample estimate
#'   now reported in `mapns` — combining the two triggers a warning.
#' @param n_splits Number of independent random respondent-level splits used
#'   when `split_sample = TRUE`. Default 200. Ignored otherwise.
#' @param prop Share of respondents assigned to the *selection* half of each
#'   split when `split_sample = TRUE`, the remainder going to the estimation
#'   half. Must be strictly between 0 and 1; default 0.3, i.e. a 30/70
#'   selection/estimation split. Shifting more data to the estimation half
#'   (lower `prop`) tightens the confidence interval on the selected pair,
#'   at the cost of a noisier — but still consistent — argmax selection
#'   step; raising it stabilises which pair wins at the cost of a wider
#'   interval. Matches the `prop` argument of
#'   \code{\link{mapns_split_test}}, so a diagnostic run and an estimation
#'   run can be held to the same split geometry. Ignored otherwise.
#' @param seed Integer seed making the whole call reproducible, default 123.
#'   This covers every stochastic step: the respondent-level splits drawn
#'   when `split_sample = TRUE` (without a seed, two identical calls return
#'   different MAPNS values) and the resampling draws behind
#'   `se = "parametric"` and `se = "bootstrap"`. The caller's RNG stream is
#'   saved and restored on exit, so `cj_apns()` never disturbs randomness
#'   elsewhere in the session. Pass `seed = NULL` to draw from the ambient
#'   stream instead — results then vary call to call.
#'
#' @return An object of class `"cj_apns"` with components:
#'   * `estimand`, `assumption`, `se_method`: as requested.
#'   * `amce`: AMCE estimates per attribute.
#'   * `apns`, `mapns`: causal attribution estimates (if requested).
#'   * `acmce`, `pi_hat`: conditional AMCEs and group shares (if conditional).
#'   * `mapns_coverage`: named vector (per attribute), only present when
#'     `assumption` is `"conditional"` or `"both"`. The share (in \[0, 1\])
#'     of the preference-ordering-group population (by \eqn{\pi_g} mass)
#'     that had informative-task coverage and so contributed to the
#'     reported conditional MAPNS. Groups lacking coverage are excluded and
#'     the remaining groups' weights renormalized (rather than treating an
#'     unidentified group's contribution as zero, which would be downward
#'     biased); `mapns_coverage` tells you how much of the population that
#'     renormalization is silent about. Equals 1 when every group had
#'     coverage.
#'   * `missing_groups`: data.frame (`item`, `contrast`, `pi_g`)
#'     listing every preference-ordering group with zero informative-task
#'     observations for its own extreme pair — excluded from the
#'     conditional MAPNS sum entirely. `NULL` if none, or if
#'     `assumption = "separability"`.
#'   * `capns_table`: data.frame (`item`, `contrast`, `estimate`, `pi_g`,
#'     `n_informative`) listing every conditional CAPNS estimate that *was*
#'     identified, across all attributes, sorted by `n_informative`
#'     descending. Use alongside `missing_groups` to see the full group
#'     structure feeding each attribute's conditional MAPNS. `NULL` under
#'     `assumption = "separability"`.
#'   * `mapns_split`: data.frame (`item`, `plugin_mapns`, `mapns`, `se`,
#'     `lower`, `upper`, `n_valid_splits`, `n_splits`), only present when
#'     `split_sample = TRUE`. `plugin_mapns` is the ordinary full-sample
#'     plug-in max; `mapns` is the split-sample estimate now stored in the
#'     main `mapns` output for that attribute; `se`/`lower`/`upper` are a
#'     Rubin's-rules combination of within- and between-split variance (see
#'     \code{\link{.split_sample_mapns}}), `NA` if fewer than two splits
#'     yielded a usable within-split SE; `n_valid_splits` counts how many of
#'     the `n_splits` replicates yielded both a valid point estimate and a
#'     valid, finite within-split SE (low values mean the check itself is
#'     underpowered for that attribute, not that there is no winner's-curse
#'     bias). `NULL` if `split_sample = FALSE`.
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
                     design = "uniform",
                     split_sample = FALSE, n_splits = 200, prop = 0.3,
                     seed = 123) {

  cl <- match.call()
  estimand <- match.arg(estimand)
  assumption <- match.arg(assumption)
  se <- match.arg(se)
  informative <- match.arg(informative)
  stopifnot(is.logical(split_sample), length(split_sample) == 1, !is.na(split_sample))
  stopifnot(is.numeric(n_splits), length(n_splits) == 1, n_splits > 0)
  n_splits <- as.integer(n_splits)
  stopifnot(is.numeric(prop), length(prop) == 1, !is.na(prop), prop > 0, prop < 1)

  # Reproducibility: seed every stochastic step of the call at once -- the
  # respondent-level splits drawn by .split_sample_mapns() below, and the
  # resampling draws in .se_parametric()/.se_bootstrap() further down. The
  # caller's RNG stream is saved here and restored on exit, so a non-NULL
  # default seed makes cj_apns() reproducible without silently fixing the
  # randomness of whatever the user runs after it.
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1, !is.na(seed))
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
    } else {
      # RNG was uninitialised before this call; set.seed() creates
      # .Random.seed, so drop it again to leave the session as we found it.
      on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
    }
    set.seed(seed)
  }

  if (split_sample && assumption == "conditional")
    warning("'split_sample = TRUE' has no effect when assumption = \"conditional\": ",
            "the winner's-curse concern it addresses only arises from the max-selection ",
            "step under separable monotonicity.", call. = FALSE)
  if (split_sample && estimand != "mapns")
    warning("'split_sample = TRUE' has no effect when estimand = \"", estimand,
            "\": it only adjusts the separability MAPNS (the argmax over pairwise APNS).",
            call. = FALSE)
  if (split_sample && se != "none")
    warning("'split_sample = TRUE' with se != \"none\": the main 'ci'/'se_detail' output is ",
            "still computed around the ordinary plug-in max, not the split-sample estimate ",
            "now reported in mapns[[attribute]]$separability (or mapns[[attribute]] under ",
            "assumption = \"separability\"). The split-sample estimate's own SE/CI (a ",
            "Rubin's-rules combination across splits) is reported separately in ",
            "'mapns_split' ('se', 'lower', 'upper' columns).", call. = FALSE)

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
                        profile_var, design, split_sample, n_splits, prop, alpha)
  .warn_mapns_not_identified(pt, assumption, preferences)

  # mapns_split: diagnostic detail when split_sample = TRUE -- the ordinary
  # plug-in max alongside the split-based estimate now reported in mapns,
  # its own SE/CI (Rubin's-rules combination across splits -- see
  # .split_sample_mapns()), and how many of the n_splits replicates yielded
  # both a valid point estimate and a valid, finite within-split SE.
  mapns_split <- NULL
  if (!is.null(pt$split_detail) && length(pt$split_detail) > 0) {
    mapns_split <- do.call(rbind, lapply(names(pt$split_detail), function(a) {
      sd  <- pt$split_detail[[a]]
      val <- pt$mapns[[a]]
      data.frame(item = a, plugin_mapns = sd$plugin_mapns,
                 mapns = if (is.list(val)) val$separability else val,
                 se = sd$se, lower = sd$lower, upper = sd$upper,
                 n_valid_splits = sd$n_valid_splits, n_splits = sd$n_splits,
                 stringsAsFactors = FALSE)
    }))
    row.names(mapns_split) <- NULL
  }

  # ---- conditional-group diagnostics -------------------------------------
  # missing_groups: ordering groups with zero informative-task coverage,
  # excluded from the renormalized conditional MAPNS (see .renormalized_mapns).
  # capns_table: every group with an identified CAPNS, its share, and its
  # informative-task sample size, across all attributes.
  missing_groups <- NULL; capns_table <- NULL
  if (!is.null(pt$cond_groups)) {
    all_records <- do.call(rbind, lapply(names(pt$cond_groups), function(a) {
      recs <- pt$cond_groups[[a]]
      if (length(recs) == 0) return(NULL)
      do.call(rbind, lapply(recs, function(r) data.frame(
        item = a, contrast = r$contrast, pi_g = r$pi_g,
        estimate = r$estimate, n_respondents = r$n_respondents,
        n_informative = r$n_informative, stringsAsFactors = FALSE)))
    }))
    if (!is.null(all_records)) {
      missing_groups <- all_records[is.na(all_records$estimate),
                                    c("item", "contrast", "pi_g")]
      row.names(missing_groups) <- NULL

      capns_table <- all_records[!is.na(all_records$estimate),
                                 c("item", "contrast", "estimate", "pi_g", "n_informative")]
      capns_table <- capns_table[order(-capns_table$n_informative), ]
      row.names(capns_table) <- NULL
    }
  }
  mapns_coverage <- if (!is.null(pt$mapns_coverage) && length(pt$mapns_coverage) > 0)
    unlist(pt$mapns_coverage) else NULL

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
         mapns_coverage = mapns_coverage,
         missing_groups = missing_groups,
         capns_table = capns_table,
         mapns_split = mapns_split,
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

#' Renormalized weighted mean of conditional group CAPNS estimates
#'
#' Implements a "drop and renormalize" treatment of ordering groups with no
#' informative-task coverage: groups with an NA estimate are excluded from
#' the weighted sum entirely, and the remaining groups' shares are
#' renormalized to sum to 1, rather than treating an unidentified group's
#' contribution as zero. Zero-filling is downward biased -- CAPNS is always
#' non-negative, so forcing an unidentified group's unknown true value to 0
#' can only pull the aggregate down, never up (see \code{mapns_sparsity_issues.tex}).
#'
#' Note this changes the estimand slightly relative to the full-population
#' MAPNS: it targets the sub-population whose ordering group has
#' informative-task coverage in this sample. \code{coverage} reports what
#' share of that population (by \eqn{\pi_g} mass) was actually used, so
#' callers can judge how much of the population the estimate is silent about.
#'
#' @param pi_vec Numeric vector of group shares (need not sum to 1; only
#'   groups with a non-NA estimate are used in the renormalized denominator).
#' @param est_vec Numeric vector of (already absolute-valued) group CAPNS
#'   estimates, same length/order as \code{pi_vec}, with NA for groups
#'   lacking informative-task coverage.
#' @return A list with \code{mapns} (renormalized weighted mean, or NA if
#'   every group is missing or \code{pi_vec} sums to 0) and \code{coverage}
#'   (share of \code{pi_vec}'s total mass with a non-NA estimate, in [0,1]).
#' @keywords internal
.renormalized_mapns <- function(pi_vec, est_vec) {
  obs <- !is.na(est_vec) & !is.na(pi_vec)
  total_pi <- sum(pi_vec, na.rm = TRUE)
  if (!any(obs) || total_pi <= 0)
    return(list(mapns = NA_real_, coverage = 0))
  observed_pi <- sum(pi_vec[obs])
  list(
    mapns = sum(pi_vec[obs] * est_vec[obs]) / observed_pi,
    coverage = observed_pi / total_pi
  )
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
                            profile_var = NULL, design = "uniform",
                            split_sample = FALSE, n_splits = 200,
                            prop = 0.3, alpha = 0.05) {

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
  cond_groups <- list(); mapns_coverage <- list()
  split_detail <- list()
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

      # Optional winner's-curse mitigation (Jensen's inequality: the plug-in
      # max is upward-biased for the true max, worst when several pairs are
      # near-tied). Replaces mapns_sep with the median, across n_splits
      # respondent-level splits, of estimating an argmax pair (selected on
      # one random half) on the other, independent half -- see
      # .split_sample_mapns() / mapns_split_test() for the method and its
      # own limitations (a bias-reduction device, not an unbiased estimator).
      # Individual pairwise apns_sep entries are left untouched: winner's
      # curse only affects the max operation, not the pairwise estimates
      # feeding it.
      if (split_sample && estimand == "mapns") {
        ss <- .split_sample_mapns(data, outcome_var, a, levs, id_var, task_var,
                                   informative, profile_var, weights_col, n_splits,
                                   prop = prop, alpha = alpha)
        split_detail[[a]] <- list(plugin_mapns = mapns_sep,
                                   se = ss$se, lower = ss$lower, upper = ss$upper,
                                   n_valid_splits = ss$n_valid_splits,
                                   n_splits = ss$n_splits)
        mapns_sep <- ss$estimate
      }
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
        apns_cond <- list(); acmce_a <- list(); group_records <- list()

        for (ep in unique_eps) {
          in_pref  <- !is.na(ep_label) & ep_label == ep
          pi_g     <- mean(in_pref, na.rm = TRUE)
          ep_top   <- top_vals[which(in_pref)[1]]
          ep_bot   <- bot_vals[which(in_pref)[1]]
          n_resp   <- sum(in_pref)

          data$.pg <- as.integer(!is.na(ep_row) & ep_row == ep)
          dm_g     <- data[data$.pg == 1L, , drop = FALSE]

          # Direct pairwise estimate at this group's own extreme pair
          # (ep_top, ep_bot), restricted to informative tasks (Proposition 8).
          v_g <- NA_real_; n_inf <- 0L
          if (nrow(dm_g) > 0) {
            pw_g <- estimate_pairwise(dm_g, outcome_var, a, ep_top, ep_bot,
                                      id_var = id_var, task_var = task_var,
                                      informative = informative,
                                      profile_var = profile_var,
                                      weights_col = weights_col)
            v_g   <- pw_g$estimate
            n_inf <- if (!is.na(pw_g$n_tq)) pw_g$n_tq + pw_g$n_tp else 0L
          }

          group_records[[ep]] <- list(
            contrast = ep, pi_g = pi_g,
            estimate = if (is.na(v_g)) NA_real_ else abs(v_g),
            n_respondents = n_resp, n_informative = n_inf)

          if (!is.na(v_g)) {
            apns_cond[[ep]] <- list(tq = ep_top, tp = ep_bot,
              estimate = pi_g * abs(v_g), assumption = "conditional")
            acmce_a[[ep]] <- list(estimate = v_g, pi = pi_g,
                                  tq = ep_top, tp = ep_bot)
          }
        }
        data$.pg <- NULL
        cond_groups[[a]] <- group_records

        if (length(apns_cond) == 0) {
          do_cond_a <- FALSE
        } else {
          # Proposition 8, renormalized (see .renormalized_mapns): groups
          # with no informative-task coverage are excluded rather than
          # zero-filled, and the remaining groups' shares are renormalized.
          agg <- .renormalized_mapns(
            sapply(group_records, `[[`, "pi_g"),
            sapply(group_records, `[[`, "estimate"))
          mapns_cond <- agg$mapns
          mapns_coverage[[a]] <- agg$coverage
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
            n_resp_grp <- stats::setNames(
              sapply(pref_groups, function(g) length(unique(sub_dm_list[[g]][[id_var]]))),
              pref_groups
            )

            apns_cond <- list(); acmce_a <- list(); group_records <- list()
            for (q in seq_along(levs)) for (p in seq_along(levs)) {
              if (q >= p) next
              pair <- paste0(levs[q], " vs ", levs[p])
              # Direct pairwise estimate for this pair, within each
              # preference-group subclass (generalises Theorem 2 to K groups).
              grp_amces <- stats::setNames(rep(NA_real_, length(pref_groups)), pref_groups)
              n_inf_grp  <- stats::setNames(rep(0L, length(pref_groups)), pref_groups)
              for (g in pref_groups) {
                sub_dm <- sub_dm_list[[g]]
                if (nrow(sub_dm) == 0) next
                pw <- estimate_pairwise(sub_dm, outcome_var, a, levs[q], levs[p],
                                        id_var = id_var, task_var = task_var,
                                        informative = informative,
                                        profile_var = profile_var,
                                        weights_col = weights_col)
                grp_amces[g] <- pw$estimate
                n_inf_grp[g] <- if (!is.na(pw$n_tq)) pw$n_tq + pw$n_tp else 0L
              }
              # Renormalized (see .renormalized_mapns): groups with no
              # informative-task coverage for this pair are excluded rather
              # than zero-filled.
              agg <- .renormalized_mapns(pi_vals, abs(grp_amces))
              apns_cond[[pair]] <- list(tq = levs[q], tp = levs[p],
                estimate = agg$mapns, assumption = "conditional")
              acmce_a[[pair]] <- list(groups = grp_amces, pi = pi_vals)
              if (Dl == 2) mapns_coverage[[a]] <- agg$coverage

              for (g in pref_groups) {
                group_records[[paste0(pair, " [", g, "]")]] <- list(
                  contrast = paste0(pair, " (group: ", g, ")"), pi_g = pi_vals[[g]],
                  estimate = if (is.na(grp_amces[[g]])) NA_real_ else abs(grp_amces[[g]]),
                  n_respondents = n_resp_grp[[g]], n_informative = n_inf_grp[[g]])
              }
            }
            mapns_cond <- .mapns_from_pairwise_cond(apns_cond, Dl)
            pi_hat[[a]] <- pi_vals; acmce[[a]] <- acmce_a
            cond_groups[[a]] <- group_records
          } else {
            pi_val <- mean(dm$.pg, na.rm = TRUE)
            dm_pro <- dm[dm$.pg == 1, , drop = FALSE]
            dm_con <- dm[dm$.pg == 0, , drop = FALSE]
            n_resp_pro <- length(unique(dm_pro[[id_var]]))
            n_resp_con <- length(unique(dm_con[[id_var]]))

            apns_cond <- list(); acmce_a <- list(); group_records <- list()
            for (q in seq_along(levs)) for (p in seq_along(levs)) {
              if (q >= p) next
              pair <- paste0(levs[q], " vs ", levs[p])
              # Direct pairwise estimate for this pair, within the pro/con
              # subclass (Proposition 2/8), not reconstructed from base-level
              # regression coefficients.
              pw_pro <- estimate_pairwise(dm_pro, outcome_var, a, levs[q], levs[p],
                                          id_var = id_var, task_var = task_var,
                                          informative = informative,
                                          profile_var = profile_var,
                                          weights_col = weights_col)
              pw_con <- estimate_pairwise(dm_con, outcome_var, a, levs[q], levs[p],
                                          id_var = id_var, task_var = task_var,
                                          informative = informative,
                                          profile_var = profile_var,
                                          weights_col = weights_col)
              v_pro <- pw_pro$estimate; v_con <- pw_con$estimate
              n_inf_pro <- if (!is.na(pw_pro$n_tq)) pw_pro$n_tq + pw_pro$n_tp else 0L
              n_inf_con <- if (!is.na(pw_con$n_tq)) pw_con$n_tq + pw_con$n_tp else 0L

              # Renormalized (see .renormalized_mapns): if either the pro or
              # con group has no informative-task coverage for this pair, it
              # is excluded rather than zero-filled.
              agg <- .renormalized_mapns(
                c(pi_val, 1 - pi_val),
                c(if (is.na(v_pro)) NA_real_ else abs(v_pro),
                  if (is.na(v_con)) NA_real_ else abs(v_con)))
              apns_cond[[pair]] <- list(tq = levs[q], tp = levs[p],
                estimate = agg$mapns, assumption = "conditional")
              acmce_a[[pair]] <- list(pro = v_pro, con = v_con, pi = pi_val)
              if (Dl == 2) mapns_coverage[[a]] <- agg$coverage

              group_records[[paste0(pair, " [pro]")]] <- list(
                contrast = paste0(pair, " (pro: ", levs[q], " > ", levs[p], ")"),
                pi_g = pi_val, estimate = if (is.na(v_pro)) NA_real_ else abs(v_pro),
                n_respondents = n_resp_pro, n_informative = n_inf_pro)
              group_records[[paste0(pair, " [con]")]] <- list(
                contrast = paste0(pair, " (con: ", levs[p], " > ", levs[q], ")"),
                pi_g = 1 - pi_val, estimate = if (is.na(v_con)) NA_real_ else abs(v_con),
                n_respondents = n_resp_con, n_informative = n_inf_con)
            }
            mapns_cond <- .mapns_from_pairwise_cond(apns_cond, Dl)
            pi_hat[[a]] <- pi_val; acmce[[a]] <- acmce_a
            cond_groups[[a]] <- group_records
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
       acmce = if (do_cond) acmce else NULL,
       cond_groups = if (do_cond) cond_groups else NULL,
       mapns_coverage = if (do_cond) mapns_coverage else NULL,
       split_detail = if (split_sample && length(split_detail) > 0) split_detail else NULL)
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
