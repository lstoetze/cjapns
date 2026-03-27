#' Print method for cj_apns objects
#'
#' @param x An object of class `"cj_apns"`.
#' @param ... Ignored.
#' @export
print.cj_apns <- function(x, ...) {
  cat("Causal Attribution in Conjoint Analysis\n")
  cat(strrep("-", 50), "\n")
  cat("Estimand:   ", x$estimand, "\n")
  cat("Assumption: ", x$assumption, "\n")
  cat("SE method:  ", x$se_method, "\n\n")

  if (x$estimand == "amce") {
    cat("Average Marginal Component Effects:\n\n")
    for (a in names(x$amce)) {
      .print_amce_attr(x, a)
    }
  }

  if (x$estimand %in% c("epns", "eapns") && !is.null(x$eapns)) {
    cat("Expected Avg. Probability of Necessary & Sufficient (EAPNS):\n\n")
    for (a in names(x$eapns)) {
      val <- x$eapns[[a]]
      if (x$assumption == "both" && is.list(val)) {
        cat("  ", a, ":\n", sep = "")
        for (asn in names(val)) if (!is.null(val[[asn]]))
          cat(sprintf("    %-22s %7.4f%s\n", asn, val[[asn]],
                      .ci_str(x, paste0("eapns.", asn, ".", a))))
      } else {
        cat(sprintf("  %-27s %7.4f%s\n", a, val,
                    .ci_str(x, paste0("eapns.", x$assumption, ".", a))))
      }
    }
    cat("\n")
  }

  if (!is.null(x$pi_hat) && length(x$pi_hat) > 0) {
    cat("Preference group shares (pi_hat):\n")
    for (a in names(x$pi_hat))
      cat(sprintf("  %-27s pro: %.3f  con: %.3f\n",
                  a, x$pi_hat[[a]], 1 - x$pi_hat[[a]]))
    cat("\n")
  }
  invisible(x)
}


#' Summary method for cj_apns objects
#'
#' Prints a detailed summary including pairwise EPNS and conditional AMCEs.
#'
#' @param object An object of class `"cj_apns"`.
#' @param ... Ignored.
#' @export
summary.cj_apns <- function(object, ...) {
  cat("Causal Attribution in Conjoint Analysis\n")
  cat(strrep("=", 50), "\n")
  cat("Call:\n"); print(object$call); cat("\n")
  print.cj_apns(object)

  if (object$estimand == "eapns" && !is.null(object$epns)) {
    cat("Pairwise EPNS:\n")
    for (a in names(object$epns)) {
      cat("  ", a, ":\n", sep = "")
      pairs <- object$epns[[a]]
      if (object$assumption == "both") {
        for (asn in names(pairs)) if (!is.null(pairs[[asn]])) {
          cat("    [", asn, "]\n")
          for (p in names(pairs[[asn]]))
            cat(sprintf("      %-20s  %.4f\n", p, pairs[[asn]][[p]]$estimate))
        }
      } else {
        for (p in names(pairs))
          cat(sprintf("    %-20s  %.4f\n", p, pairs[[p]]$estimate))
      }
    }
    cat("\n")
  }

  if (!is.null(object$acmce)) {
    cat("Conditional AMCEs by preference group:\n")
    for (a in names(object$acmce)) {
      cat("  ", a, ":\n", sep = "")
      for (p in names(object$acmce[[a]])) {
        v <- object$acmce[[a]][[p]]
        cat(sprintf("    %-20s  pro: %7.4f  con: %7.4f  pi: %.3f\n",
                    p, v$pro, v$con, v$pi))
      }
    }
    cat("\n")
  }
  invisible(object)
}


#' Convert cj_apns results to a data frame
#'
#' Returns a tidy data frame of all estimates with standard errors and
#' confidence intervals (if available), suitable for custom ggplot2 plots.
#'
#' @param x An object of class `"cj_apns"`.
#' @param ... Ignored.
#' @return A data.frame with columns: `attribute`, `estimand`, `assumption`,
#'   `comparison`, `estimate`, `se`, `lower`, `upper`.
#'
#' @examples
#' \dontrun{
#' data(cnj_cand)
#' res <- cj_apns(vote ~ borders + schools, data = cnj_cand,
#'                id = ~ ResponseId, se = "none")
#' as.data.frame(res)
#' }
#'
#' @export
as.data.frame.cj_apns <- function(x, ...) {
  rows <- list()

  # AMCEs
  for (a in names(x$amce)) {
    amce_a <- x$amce[[a]]
    for (i in seq_along(amce_a$estimate)) {
      lev <- names(amce_a$estimate)[i]
      nm <- paste0("amce.", a, ".", lev)
      rows <- c(rows, list(.make_row(a, "amce", NA, paste(lev, "vs", amce_a$base_level),
                                      amce_a$estimate[i], x, nm)))
    }
  }

  # EAPNS
  if (!is.null(x$eapns)) for (a in names(x$eapns)) {
    val <- x$eapns[[a]]
    if (x$assumption == "both" && is.list(val)) {
      for (asn in names(val)) if (!is.null(val[[asn]])) {
        nm <- paste0("eapns.", asn, ".", a)
        rows <- c(rows, list(.make_row(a, "eapns", asn, "all pairs",
                                        val[[asn]], x, nm)))
      }
    } else {
      nm <- paste0("eapns.", x$assumption, ".", a)
      rows <- c(rows, list(.make_row(a, "eapns", x$assumption, "all pairs",
                                      val, x, nm)))
    }
  }

  # EPNS
  if (!is.null(x$epns)) for (a in names(x$epns)) {
    pairs <- x$epns[[a]]
    if (x$assumption == "both") {
      for (asn in names(pairs)) if (!is.null(pairs[[asn]])) {
        for (p in names(pairs[[asn]])) {
          nm <- paste0("epns.", asn, ".", a, ".", p)
          rows <- c(rows, list(.make_row(a, "epns", asn, p,
                                          pairs[[asn]][[p]]$estimate, x, nm)))
        }
      }
    } else {
      for (p in names(pairs)) {
        nm <- paste0("epns.", x$assumption, ".", a, ".", p)
        rows <- c(rows, list(.make_row(a, "epns", x$assumption, p,
                                        pairs[[p]]$estimate, x, nm)))
      }
    }
  }

  # Conditional AMCEs
  if (!is.null(x$acmce)) for (a in names(x$acmce)) if (!is.null(x$acmce[[a]])) {
    for (p in names(x$acmce[[a]])) {
      v <- x$acmce[[a]][[p]]
      nm_pro <- paste0("acmce.pro.", a, ".", p)
      nm_con <- paste0("acmce.con.", a, ".", p)
      rows <- c(rows, list(.make_row(a, "acmce_pro", "conditional", p,
                                      v$pro, x, nm_pro)))
      rows <- c(rows, list(.make_row(a, "acmce_con", "conditional", p,
                                      v$con, x, nm_con)))
    }
  }

  do.call(rbind, rows)
}


#' Plot method for cj_apns objects
#'
#' Produces a dot-and-whisker plot of EAPNS, EPNS, or AMCE estimates,
#' replicating the style of Figure 2 in Stoetzer & Magazinnik (2026).
#' Requires \pkg{ggplot2}.
#'
#' @param x An object of class `"cj_apns"`.
#' @param what Which estimand to plot: `"eapns"`, `"epns"`, or `"amce"`.
#'   Defaults to the estimand stored in `x`.
#' @param ... Passed to ggplot2 functions.
#' @return A ggplot2 object (invisibly).
#'
#' @examples
#' \dontrun{
#' data(cnj_cand)
#' res <- cj_apns(vote ~ borders + eurobonds + immucard + schools,
#'                data = cnj_cand, id = ~ ResponseId,
#'                se = "parametric", B = 200)
#' plot(res)
#' }
#'
#' @export
plot.cj_apns <- function(x, what = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.")

  if (is.null(what)) what <- x$estimand
  df <- as.data.frame(x)
  df <- df[df$estimand == what, , drop = FALSE]
  if (nrow(df) == 0) stop("No estimates for estimand = '", what, "'.")

  gg <- ggplot2::ggplot

  label_col <- if (what == "eapns") "attribute" else "comparison"
  df$label <- df[[label_col]]

  has_ci <- !all(is.na(df$lower))

  if (x$assumption == "both" && !all(is.na(df$assumption))) {
    p <- gg(df, ggplot2::aes(x = .data$estimate, y = .data$label,
                              shape = .data$assumption, col = .data$assumption))
  } else {
    p <- gg(df, ggplot2::aes(x = .data$estimate, y = .data$label))
  }

  p <- p + ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(0.3))

  if (has_ci) {
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
      height = 0.15, position = ggplot2::position_dodge(0.3))
  }

  p <- p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
    ggplot2::scale_color_grey() +
    ggplot2::labs(x = "Estimate", y = "",
                  title = switch(what,
                    eapns = "Expected Avg. Probability of Necessary & Sufficient",
                    epns = "Expected Probability of Necessary & Sufficient",
                    amce = "Average Marginal Component Effects")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = "bottom")

  if (!is.null(df$attribute) && length(unique(df$attribute)) > 1 && what != "eapns")
    p <- p + ggplot2::facet_wrap(~ attribute, scales = "free_y")

  print(p)
  invisible(p)
}


# ── helpers ──────────────────────────────────────────────────────────────────

#' @keywords internal
.print_amce_attr <- function(x, a) {
  amce_a <- x$amce[[a]]
  cat("  ", a, " (base: ", amce_a$base_level, ")\n", sep = "")
  for (i in seq_along(amce_a$estimate))
    cat(sprintf("    %-20s  %7.4f  (SE: %.4f)%s\n",
                names(amce_a$estimate)[i], amce_a$estimate[i], amce_a$se[i],
                .ci_str(x, paste0("amce.", a, ".", names(amce_a$estimate)[i]))))
  cat("\n")
}

#' @keywords internal
.ci_str <- function(x, nm) {
  if (is.null(x$ci) || !nm %in% x$ci$parameter) return("")
  r <- x$ci[x$ci$parameter == nm, ]
  sprintf("  [%.4f, %.4f]", r$lower, r$upper)
}

#' @keywords internal
.make_row <- function(attribute, estimand, assumption, comparison, estimate, x, nm) {
  se_val <- NA_real_; lo <- NA_real_; hi <- NA_real_
  if (!is.null(x$ci) && nm %in% x$ci$parameter) {
    r <- x$ci[x$ci$parameter == nm, ]
    se_val <- r$se; lo <- r$lower; hi <- r$upper
  }
  data.frame(attribute = attribute, estimand = estimand, assumption = assumption,
             comparison = comparison, estimate = unname(estimate),
             se = se_val, lower = lo, upper = hi,
             stringsAsFactors = FALSE, row.names = NULL)
}
