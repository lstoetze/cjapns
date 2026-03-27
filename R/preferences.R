#' Create preference data for conditional separable monotonicity
#'
#' Converts respondent-level preference data (e.g., from a pre-conjoint
#' survey question) into the format needed for APNS estimation under
#' conditional separable monotonicity. This maps each respondent's stated
#' preference direction to the binary group indicator \eqn{p^{qp}_{il}}.
#'
#' @param data A data.frame with one row per respondent containing the
#'   respondent ID and one column per attribute with preference values.
#' @param id A one-sided formula specifying the respondent ID variable
#'   (e.g., `~ ResponseId`).
#' @param type How preferences are encoded:
#'   * `"binary"`: Values are `"pro"`/`"con"`, `1`/`0`, or `TRUE`/`FALSE`.
#'   * `"scale"`: Numeric values (e.g., Likert), dichotomised at `cutpoint`.
#' @param cutpoint For `type = "scale"`, the threshold at or above which a
#'   respondent is classified as "pro". Defaults to the midpoint of the range.
#'
#' @return An object of class `"cj_preferences"`.
#'
#' @examples
#' \dontrun{
#' data(cnj_cand)
#'
#' # Extract preference data: one row per respondent
#' pref_dat <- cnj_cand[!duplicated(cnj_cand$ResponseId),
#'                       c("ResponseId", "Att_borders", "Att_eurobonds",
#'                         "Att_immucard", "Att_schools", "Att_tracingapp")]
#'
#' # Strip Att_ prefix so names match the conjoint attributes
#' names(pref_dat) <- gsub("^Att_", "", names(pref_dat))
#'
#' prefs <- make_preferences(pref_dat, id = ~ ResponseId, type = "binary")
#' prefs
#' }
#'
#' @export
make_preferences <- function(data, id, type = c("binary", "scale"),
                             cutpoint = NULL) {
  type <- match.arg(type)
  id_var <- all.vars(id)
  stopifnot(length(id_var) == 1, id_var %in% names(data))

  attr_cols <- setdiff(names(data), id_var)
  if (length(attr_cols) == 0) stop("No attribute columns found.")

  groups <- list()
  for (a in attr_cols) {
    vals <- data[[a]]
    if (type == "binary") {
      if (is.factor(vals)) vals <- as.character(vals)
      if (is.character(vals)) {
        vals_l <- tolower(vals)
        indicator <- ifelse(
          vals_l %in% c("pro", "1", "true", "yes", "favour", "in favour"),
          1L, ifelse(is.na(vals), NA_integer_, 0L))
      } else {
        indicator <- as.integer(as.logical(vals))
      }
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
  for (a in x$attributes) {
    v <- x$data[[a]]
    cat("    ", a, ": ", sum(v == 1, na.rm = TRUE), " pro / ",
        sum(v == 0, na.rm = TRUE), " con\n", sep = "")
  }
  invisible(x)
}


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
#' @return An object of class `"cj_design"`.
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
    cat("  ", a, ": ", paste(names(p), "=", round(p, 3), collapse = ", "),
        "\n")
  }
  if (!is.null(x$constraints))
    cat("  Constraints:", length(x$constraints), "restricted profiles\n")
  invisible(x)
}
