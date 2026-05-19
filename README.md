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

The APNS estimand captures the probability that a change in an attribute level is both *necessary* and *sufficient* for a respondent's choice — a direct measure of attribute relevance. Under separable monotonicity, the MAPNS simplifies to the average absolute AMCE. Under the more realistic conditional separable monotonicity, it is a weighted sum of absolute conditional AMCEs by preference group.

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

# ── 2. Estimate MAPNS under both assumptions ─────────────────────────
#    tasks = ~ time restricts AMCEs to informative tasks (Definition 3)

res <- cj_apns(
  vote ~ borders + eurobonds + immucard + schools + tracingapp,
  data        = cnj_cand,
  id          = ~ ResponseId,
  tasks       = ~ time,
  profile     = ~ profile,
  assumption  = "both",
  preferences = prefs,
  se          = "parametric",
  B           = 500
)

# View results
summary(res)

# ── 3. Plot ──────────────────────────────────────────────────────────

plot(res)

# Or build a custom ggplot from the tidy data frame
library(ggplot2)
df <- as.data.frame(res)
df_mapns <- df[df$estimand == "mapns", ]
df_mapns$assumption <- factor(df_mapns$assumption,
  levels = c("conditional", "separability"),
  labels = c("APNS cond. sep. mono.", "APNS sep. mono."))

ggplot(df_mapns, aes(x = reorder(attribute, estimate), y = estimate,
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
| `cj_apns()` | Estimate AMCE, APNS, or MAPNS with multiple SE methods |
| `make_preferences()` | Prepare preference data for conditional monotonicity |
| `make_design()` | Specify non-uniform randomization designs |
| `print()`, `summary()` | Display results |
| `plot()` | Dot-and-whisker plot (requires ggplot2) |
| `as.data.frame()` | Tidy export for custom plots |

## Informative task filtering

By default `cj_apns()` restricts AMCE estimation to *informative tasks* (Definition 3 in the paper): for each level pair (tq, tp), only tasks where exactly one profile shows tq and all remaining profiles show tp (or vice versa) are used. This is required for Proposition 1 to apply.

Supply the task-number and profile variables:

```r
res <- cj_apns(formula, data = dat, id = ~ id, tasks = ~ time, profile = ~ profile)
```

To use all tasks instead (e.g. for comparison):

```r
res_all <- cj_apns(formula, data = dat, id = ~ id,
                   tasks = ~ time, profile = ~ profile, informative = "all")
```

If `informative = "informative"` (the default) but `tasks` is not supplied, the function falls back to `"all"` with a warning.

## Standard error methods

| Method | `se =` | Description |
|--------|--------|-------------|
| Parametric bootstrap | `"parametric"` (default) | Resamples AMCEs from N(est, se²), applies plugin; CIs from empirical quantiles of bootstrap draws (Algorithm 1). |
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
