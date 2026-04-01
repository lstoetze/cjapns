# ══════════════════════════════════════════════════════════════════════════════
# Ranking data preprocessing
# ══════════════════════════════════════════════════════════════════════════════

#' Parse a rank value to integer (handles strings like "1 - höchste Präferenz")
#' @keywords internal
.parse_rank_value <- function(x) {
  as.integer(sub("^(\\d+).*", "\\1", trimws(as.character(x))))
}


#' Preprocess raw ranking survey data into a structured ranking object
#'
#' Converts raw survey ranking columns into a \code{"cj_rank_data"} object
#' that can be passed directly to \code{make_preferences}. This is the
#' recommended preprocessing step before calling \code{make_preferences} with
#' \code{type = "ranking"} (Proposition 8).
#'
#' @param data A data.frame with \strong{one row per respondent} containing the
#'   ID variable and all raw ranking columns. If multiple rows per respondent
#'   are present (e.g. the full long-format conjoint data), deduplicate first
#'   with e.g. \code{data[!duplicated(data$respondent_id), ]}.
#' @param id A one-sided formula specifying the respondent ID variable
#'   (e.g. \code{~ respondent_id}).
#' @param rank_spec A named list with one entry per conjoint attribute. The
#'   \strong{name} of each entry must match the attribute column name used in
#'   the conjoint formula (e.g. \code{"geschlecht"}, not \code{"gender"}).
#'   Each entry is a list with two fields:
#'   \describe{
#'     \item{\code{rank_cols}}{Named character vector mapping level labels to
#'       raw column names in \code{data}:
#'       \code{c("Weiblich" = "rank_genderr1", "Maennlich" = "rank_genderr2")}.
#'       The names become the level labels used throughout the analysis.}
#'     \item{\code{rank_values}}{(Optional) Character or numeric vector
#'       specifying the possible raw values in the rank columns \emph{in order
#'       from most preferred (rank 1) to least preferred}. The position of a
#'       raw value in this vector becomes its integer rank. If omitted, the
#'       leading integer is parsed from the raw value automatically (handles
#'       strings such as \code{"1 - höchste Präferenz"}).}
#'   }
#'
#' @return An object of class \code{"cj_rank_data"} with components:
#'   \describe{
#'     \item{\code{rank_data}}{Named list of data frames, one per attribute.
#'       Each data frame has the ID column plus one column per level containing
#'       integer rank positions (1 = most preferred).}
#'     \item{\code{id_var}}{Name of the respondent ID variable.}
#'     \item{\code{attributes}}{Character vector of attribute names.}
#'   }
#'
#' @examples
#' \dontrun{
#' rank_spec <- list(
#'   geschlecht = list(
#'     rank_cols   = c("Weiblich" = "rank_genderr1", "Maennlich" = "rank_genderr2"),
#'     rank_values = 1:2
#'   ),
#'   asylzeugnis = list(
#'     rank_cols   = c("Keine Unstimmigkeiten"     = "rank_asyltestr1",
#'                     "Geringe Unstimmigkeiten"   = "rank_asyltestr2",
#'                     "Groessere Unstimmigkeiten" = "rank_asyltestr3"),
#'     rank_values = c("1 - höchste Präferenz", "2", "3 - niedrigste Präferenz")
#'   )
#' )
#'
#' pref_dat <- tidy_ranking_data(
#'   data      = conjoint_df[!duplicated(conjoint_df$respondent_id), ],
#'   id        = ~ respondent_id,
#'   rank_spec = rank_spec
#' )
#' pref_dat
#' }
#'
#' @export
tidy_ranking_data <- function(data, id, rank_spec) {
  stopifnot(is.data.frame(data), is.list(rank_spec), !is.null(names(rank_spec)))
  id_var <- all.vars(id)
  stopifnot(length(id_var) == 1, id_var %in% names(data))

  if (anyDuplicated(data[[id_var]]))
    stop("'data' contains duplicate values of '", id_var, "'. ",
         "Please pass one row per respondent, e.g. ",
         "data[!duplicated(data$", id_var, "), ].")

  rank_data <- list()
  for (a in names(rank_spec)) {
    spec      <- rank_spec[[a]]
    rank_cols <- spec$rank_cols

    if (is.null(rank_cols) || !is.character(rank_cols) || is.null(names(rank_cols)))
      stop("rank_spec[['", a, "']]$rank_cols must be a named character vector ",
           "of the form c(\"Level label\" = \"column_name\", ...).")

    missing_cols <- setdiff(rank_cols, names(data))
    if (length(missing_cols) > 0)
      stop("Column(s) not found in data for attribute '", a, "': ",
           paste(missing_cols, collapse = ", "))

    rd        <- data.frame(as.character(data[[id_var]]), stringsAsFactors = FALSE)
    names(rd) <- id_var

    for (lbl in names(rank_cols)) {
      col <- rank_cols[[lbl]]
      raw <- data[[col]]

      if (!is.null(spec$rank_values)) {
        rank_pos <- match(as.character(raw), as.character(spec$rank_values))
        n_unmatched <- sum(is.na(rank_pos) & !is.na(raw))
        if (n_unmatched > 0)
          warning(n_unmatched, " value(s) in column '", col, "' for attribute '",
                  a, "' were not found in rank_values and will be NA.")
      } else {
        rank_pos <- .parse_rank_value(raw)
      }

      rd[[lbl]] <- rank_pos
    }

    rank_data[[a]] <- rd
  }

  structure(
    list(rank_data = rank_data, id_var = id_var, attributes = names(rank_spec)),
    class = "cj_rank_data"
  )
}


#' @export
print.cj_rank_data <- function(x, ...) {
  n_resp <- nrow(x$rank_data[[x$attributes[1]]])
  cat("Conjoint Ranking Data\n")
  cat("  Respondents:", n_resp, "\n")
  cat("  Attributes: ", paste(x$attributes, collapse = ", "), "\n\n")
  for (a in x$attributes) {
    rd         <- x$rank_data[[a]]
    level_cols <- setdiff(names(rd), x$id_var)
    mean_ranks <- sapply(level_cols, function(lc) mean(rd[[lc]], na.rm = TRUE))
    cat("  ", a, " (", length(level_cols), " levels):\n", sep = "")
    for (lc in level_cols)
      cat(sprintf("    %-45s mean rank = %.2f\n", lc, mean_ranks[[lc]]))
  }
  invisible(x)
}


# ══════════════════════════════════════════════════════════════════════════════
# Preferences
# ══════════════════════════════════════════════════════════════════════════════

#' Create preference data for conditional separable monotonicity
#'
#' Converts respondent-level preference data into the format needed for APNS
#' estimation under conditional separable monotonicity.
#'
#' For \strong{multi-level attributes}, use \code{type = "ranking"} together
#' with a \code{\link{tidy_ranking_data}} object. This correctly implements
#' the pair-specific binary indicator \eqn{p^{qp}_{il}} from Proposition 8:
#' for each pair \eqn{(t_q, t_p)}, respondents who ranked \eqn{t_q} above
#' \eqn{t_p} form the pro group, and those who ranked \eqn{t_p} above
#' \eqn{t_q} form the con group.
#'
#' @param data For \code{type \%in\% c("binary", "scale", "multilevel")}: a
#'   data.frame with one row per respondent containing the ID column and one
#'   column per attribute. For \code{type = "ranking"}: a
#'   \code{"cj_rank_data"} object produced by \code{\link{tidy_ranking_data}}.
#' @param id A one-sided formula specifying the respondent ID variable
#'   (e.g. \code{~ ResponseId}).
#' @param type How preferences are encoded:
#'   \describe{
#'     \item{\code{"binary"}}{Values are \code{"pro"}/\code{"con"},
#'       \code{1}/\code{0}, or \code{TRUE}/\code{FALSE}.}
#'     \item{\code{"scale"}}{Numeric values (e.g. Likert), dichotomised at
#'       \code{cutpoint}.}
#'     \item{\code{"multilevel"}}{Raw character values stored as-is; one
#'       group per unique value. Does \emph{not} implement Proposition 8
#'       correctly for multi-level attributes.}
#'     \item{\code{"ranking"}}{Pair-specific binary indicator derived from a
#'       full level ranking. Pass a \code{"cj_rank_data"} object as
#'       \code{data}. This is the correct type for Proposition 8.}
#'   }
#' @param cutpoint For \code{type = "scale"}, the threshold at or above which
#'   a respondent is classified as "pro". Defaults to the midpoint of the
#'   observed range.
#'
#' @return An object of class \code{"cj_preferences"}.
#'
#' @examples
#' \dontrun{
#' data(cnj_cand)
#'
#' # --- Binary attributes ---
#' pref_dat <- cnj_cand[!duplicated(cnj_cand$ResponseId),
#'                       c("ResponseId", "Att_borders", "Att_eurobonds",
#'                         "Att_immucard", "Att_schools")]
#' names(pref_dat) <- gsub("^Att_", "", names(pref_dat))
#' prefs <- make_preferences(pref_dat, id = ~ ResponseId, type = "binary")
#'
#' # --- Multi-level attributes (ranking, Proposition 8) ---
#' rank_spec <- list(
#'   asylzeugnis = list(
#'     rank_cols   = c("Keine Unstimmigkeiten"     = "rank_asyltestr1",
#'                     "Geringe Unstimmigkeiten"   = "rank_asyltestr2",
#'                     "Groessere Unstimmigkeiten" = "rank_asyltestr3"),
#'     rank_values = c("1 - höchste Präferenz", "2", "3 - niedrigste Präferenz")
#'   ),
#'   geschlecht = list(
#'     rank_cols   = c("Weiblich" = "rank_genderr1", "Maennlich" = "rank_genderr2"),
#'     rank_values = 1:2
#'   )
#' )
#' pref_dat_rank <- tidy_ranking_data(
#'   conjoint_df[!duplicated(conjoint_df$respondent_id), ],
#'   id        = ~ respondent_id,
#'   rank_spec = rank_spec
#' )
#' prefs_rank <- make_preferences(pref_dat_rank, id = ~ respondent_id,
#'                                type = "ranking")
#' }
#'
#' @export
make_preferences <- function(data, id, type = c("binary", "scale", "multilevel", "ranking"),
                             cutpoint = NULL) {
  type   <- match.arg(type)
  id_var <- all.vars(id)
  stopifnot(length(id_var) == 1)

  # ── ranking type: data must be a cj_rank_data object ─────────────────────
  if (type == "ranking") {
    if (!inherits(data, "cj_rank_data"))
      stop("For type = \"ranking\", 'data' must be a \"cj_rank_data\" object ",
           "produced by tidy_ranking_data().")
    if (data$id_var != id_var)
      stop("The id variable in the cj_rank_data object ('", data$id_var,
           "') does not match the id formula ('", id_var, "').")

    out        <- data.frame(data$rank_data[[1]][[id_var]], stringsAsFactors = FALSE)
    names(out) <- id_var

    return(structure(
      list(data = out, id_var = id_var, attributes = data$attributes,
           type = "ranking", rank_data = data$rank_data),
      class = "cj_preferences"
    ))
  }

  # ── binary / scale / multilevel ──────────────────────────────────────────
  stopifnot(id_var %in% names(data))
  attr_cols <- setdiff(names(data), id_var)
  if (length(attr_cols) == 0) stop("No attribute columns found.")

  groups <- list()
  for (a in attr_cols) {
    vals <- data[[a]]
    if (type == "binary") {
      if (is.factor(vals)) vals <- as.character(vals)
      if (is.character(vals)) {
        vals_l    <- tolower(vals)
        indicator <- ifelse(
          vals_l %in% c("pro", "1", "true", "yes", "favour", "in favour"),
          1L, ifelse(is.na(vals), NA_integer_, 0L))
      } else {
        indicator <- as.integer(as.logical(vals))
      }
    } else if (type == "multilevel") {
      if (is.factor(vals)) vals <- as.character(vals)
      indicator <- as.character(vals)
    } else {
      stopifnot(is.numeric(vals))
      cp <- if (is.null(cutpoint)) {
        (min(vals, na.rm = TRUE) + max(vals, na.rm = TRUE)) / 2
      } else cutpoint
      indicator <- ifelse(vals >= cp, 1L, 0L)
    }
    groups[[a]] <- indicator
  }

  out <- data.frame(data[[id_var]], stringsAsFactors = FALSE)
  names(out) <- id_var
  for (a in attr_cols) out[[a]] <- groups[[a]]

  structure(list(data = out, id_var = id_var, attributes = attr_cols,
                 type = type),
            class = "cj_preferences")
}


#' @export
print.cj_preferences <- function(x, ...) {
  cat("Conjoint Preferences\n")
  cat("  Respondents:", nrow(x$data), "\n")
  cat("  Attributes: ", paste(x$attributes, collapse = ", "), "\n")
  cat("  Type:       ", x$type, "\n")
  if (x$type == "ranking") {
    for (a in x$attributes) {
      rd         <- x$rank_data[[a]]
      level_cols <- setdiff(names(rd), x$id_var)
      mean_ranks <- sapply(level_cols, function(lc) mean(rd[[lc]], na.rm = TRUE))
      cat("    ", a, ": ", length(level_cols), " levels",
          " (mean ranks: ",
          paste(level_cols, round(mean_ranks, 2), sep = " = ", collapse = ", "),
          ")\n", sep = "")
    }
  } else {
    for (a in x$attributes) {
      v <- x$data[[a]]
      if (is.character(v)) {
        tbl <- table(v, useNA = "no")
        cat("    ", a, ": ", paste(names(tbl), tbl, sep = " = ", collapse = ", "), "\n", sep = "")
      } else {
        cat("    ", a, ": ", sum(v == 1, na.rm = TRUE), " pro / ",
            sum(v == 0, na.rm = TRUE), " con\n", sep = "")
      }
    }
  }
  invisible(x)
}


# ══════════════════════════════════════════════════════════════════════════════
# Conjoint design
# ══════════════════════════════════════════════════════════════════════════════

#' Create a conjoint design object
#'
#' Specifies non-uniform randomization for conjoint experiments. When level
#' probabilities are not uniform or some profiles are restricted, pass this
#' object to \code{cj_apns} via the \code{design} argument.
#'
#' @param level_probs Named list of named numeric vectors. Each element is an
#'   attribute; inner names are levels, values are marginal probabilities
#'   (must sum to 1).
#' @param constraints Optional list of lists. Each inner list maps attribute
#'   names to forbidden level values for a restricted profile combination.
#'
#' @return An object of class \code{"cj_design"}.
#'
#' @examples
#' d <- make_design(
#'   level_probs = list(
#'     borders = c("0" = 0.5, "1" = 0.5),
#'     schools = c("0" = 0.5, "1" = 0.5)
#'   )
#' )
#' d
#'
#' @export
make_design <- function(level_probs, constraints = NULL) {
  stopifnot(is.list(level_probs), !is.null(names(level_probs)))
  for (a in names(level_probs)) {
    p <- level_probs[[a]]
    if (!is.numeric(p) || is.null(names(p)))
      stop("level_probs[['", a, "']] must be a named numeric vector.")
    if (abs(sum(p) - 1) > 1e-8)
      stop("level_probs[['", a, "']] must sum to 1.")
  }
  structure(list(level_probs = level_probs, constraints = constraints),
            class = "cj_design")
}

#' @export
print.cj_design <- function(x, ...) {
  cat("Conjoint Design\n")
  for (a in names(x$level_probs)) {
    p <- x$level_probs[[a]]
    cat("  ", a, ": ", paste(names(p), "=", round(p, 3), collapse = ", "), "\n")
  }
  if (!is.null(x$constraints))
    cat("  Constraints:", length(x$constraints), "restricted profiles\n")
  invisible(x)
}
