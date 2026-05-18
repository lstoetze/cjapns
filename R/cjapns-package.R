#' @keywords internal
"_PACKAGE"

# Suppress R CMD check NOTEs for ggplot2 non-standard evaluation
utils::globalVariables(c(".data", "estimate", "label", "assumption",
                          "lower", "upper", "attribute"))

#' @title cjapns: Causal Attribution in Conjoint Experiments
#'
#' @description
#' Implements the Average Probability of Necessary and Sufficient conditions
#' (APNS) for measuring attribute relevance in conjoint experiments. The main
#' function \code{cj_apns} estimates AMCEs, APNS, and MAPNS under separable or
#' conditional separable monotonicity assumptions.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{cj_apns}: Estimate AMCE, APNS, or MAPNS.
#'   \item \code{make_preferences}: Prepare preference data for conditional
#'     monotonicity.
#'   \item \code{make_design}: Specify non-uniform randomization designs.
#' }
#'
#' @references
#' Stoetzer, L.F. and Magazinnik, A. (2026). Measuring Attribute Relevance
#' in Conjoint Analysis.
#'
#' Hainmueller, J., Hopkins, D.J. and Yamamoto, T. (2014). Causal Inference
#' in Conjoint Analysis. *Political Analysis*, 22(1), 1--30.
#'
#' @name cjapns-package
#' @aliases cjapns
NULL
