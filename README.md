# cjapns: Causal Attribution in Conjoint Experiments

<!-- badges: start -->
<!-- badges: end -->

R package implementing the Average Probability of Necessary and Sufficient conditions (APNS) for measuring attribute relevance in conjoint experiments (Stoetzer and Magazinnik, 2026).

## Installation

```r
# Install from GitHub
remotes::install_github("lstoetze/cjapns")
```

## Overview

Standard conjoint analysis estimates Average Marginal Component Effects (AMCEs), which measure how changing an attribute level affects profile selection. However, the AMCE is not well suited for assessing **attribute relevance** — a zero AMCE can arise even when an attribute is highly important, if respondents disagree on the preferred direction.

The APNS estimand captures the probability that a change in an attribute level is both *necessary* and *sufficient* for a respondent's choice — a direct measure of attribute relevance. Under separable monotonicity, the EAPNS simplifies to the average absolute AMCE. Under the more realistic conditional separable monotonicity, it is a weighted sum of absolute conditional AMCEs by preference group.

## Quick start

```r
library(cjapns)
data(cnj_cand)

# ── 1. Prepare preference data (one row per respondent) ──────────────

pref_dat <- cnj_cand[!duplicated(cnj_cand$ResponseId),
                      c("ResponseId", "Att_borders", "Att_eurobonds",
                        "Att_immucard", "Att_schools", "Att_tracingapp")]
names(pref_dat) <- gsub("^Att_", "", names(pref_dat))

prefs <- make_preferences(pref_dat, id = ~ ResponseId, type = "binary")

# ── 2. Estimate EAPNS under both assumptions ─────────────────────────

res <- cj_apns(
  vote ~ borders + eurobonds + immucard + schools + tracingapp,
  data      = cnj_cand,
  id        = ~ ResponseId,
  assumption = "both",
  preferences = prefs,
  se         = "parametric",
  B          = 500
)

# View results
summary(res)

# ── 3. Plot ──────────────────────────────────────────────────────────

plot(res)

# Or build a custom ggplot from the tidy data frame
library(ggplot2)
df <- as.data.frame(res)
df_eapns <- df[df$estimand == "eapns", ]
df_eapns$assumption <- factor(df_eapns$assumption,
  levels = c("conditional", "separability"),
  labels = c("APNS cond. sep. mono.", "APNS sep. mono."))

ggplot(df_eapns, aes(x = reorder(attribute, estimate), y = estimate,
                      col = assumption,
                      ymin = lower, ymax = upper)) +
  geom_pointrange(position = position_dodge(0.3)) +
  coord_flip() +
  scale_color_grey() +
  geom_hline(yintercept = 0, col = "red", alpha = 0.3) +
  labs(x = "", y = "Estimate") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom")
```

## Main functions

| Function | Description |
|----------|-------------|
| `cj_apns()` | Estimate AMCE, EPNS, or EAPNS with multiple SE methods |
| `make_preferences()` | Prepare preference data for conditional monotonicity |
| `make_design()` | Specify non-uniform randomization designs |
| `print()`, `summary()` | Display results |
| `plot()` | Dot-and-whisker plot (requires ggplot2) |
| `as.data.frame()` | Tidy export for custom plots |

## Standard error methods

| Method | `se =` | Description |
|--------|--------|-------------|
| Parametric bootstrap | `"parametric"` (default) | Resamples AMCEs from N(est, se²), applies plugin. Fast. |
| Nonparametric bootstrap | `"bootstrap"` | Resamples respondents, re-estimates full pipeline. |
| Folded normal | `"folded_normal"` | Analytical SE from the folded normal distribution of \|AMCE\|. |
| Jackknife | `"jackknife"` | Leave-one-respondent-out. Exact but slow for large N. |
| None | `"none"` | Point estimates only. |

## Generating documentation

After editing roxygen2 comments in the `R/*.R` files:

```r
# Generate .Rd files and update NAMESPACE
devtools::document()

# Check the package
devtools::check()

# Build
devtools::build()
```

## References

Stoetzer, L.F. and Magazinnik, A. (2026). Measuring Attribute Relevance in Conjoint Analysis.

Hainmueller, J., Hopkins, D.J. and Yamamoto, T. (2014). Causal Inference in Conjoint Analysis. *Political Analysis*, 22(1), 1-30.

## License

MIT
