# Size and power of the randomization test against the paper's DGP. These
# replicate a few hundred simulated experiments, so they are skipped on CRAN.

# Rejection rate over `R` exact-test p-values, with a binomial tolerance band
# around the nominal level. The test is exact, so the rate should sit at or
# below alpha; a rate far below it would mean the statistic has collapsed.
expect_size <- function(p, alpha = 0.05, lower = 0.4, upper = 2.2) {
  rate <- mean(p <= alpha)
  expect_gt(rate, lower * alpha)
  expect_lt(rate, upper * alpha)
}

test_that("size is correct under the sharp null, separability APNS", {
  skip_on_cran()
  set.seed(101)
  p <- replicate(400, ri_p(
    sim_conjoint(N = 200, K = 6, L = 3, D = 2, w1 = 0, lambda = 2),
    estimand = "apns", assumption = "separability", B = 199))

  expect_size(p, 0.05)
  expect_size(p, 0.10)
  expect_gt(mean(p), 0.42)
  expect_lt(mean(p), 0.58)
})

test_that("size is correct under the sharp null, separability MAPNS over
           several pairs", {
  skip_on_cran()
  # D = 4 gives 6 candidate pairs, so the statistic is a genuine maximum. The
  # plug-in max is upward-biased, but the test re-takes the maximum on every
  # draw, so the reference distribution carries the same selection and size
  # stays at nominal.
  set.seed(102)
  p <- replicate(400, ri_p(
    sim_conjoint(N = 250, K = 8, L = 3, D = 4, w1 = 0, lambda = 2),
    estimand = "mapns", assumption = "separability", B = 199))

  expect_size(p, 0.05)
  expect_size(p, 0.10)
})

test_that("size is correct under the sharp null, conditional APNS", {
  skip_on_cran()
  set.seed(103)
  p <- replicate(300, {
    s <- sim_conjoint(N = 250, K = 8, L = 3, D = 2, w1 = 0, lambda = 2,
                      heterogeneous = TRUE)
    prefs <- make_preferences(s$pref_binary, id = ~ id, type = "binary")
    ri_p(s, estimand = "apns", assumption = "conditional",
         preferences = prefs, B = 199, min_informative = 0)
  })

  expect_size(p, 0.05)
  expect_size(p, 0.10)
})

test_that("size is correct under the sharp null, heterogeneous MAPNS", {
  skip_on_cran()
  set.seed(104)
  p <- replicate(300, {
    s <- sim_conjoint(N = 300, K = 10, L = 3, D = 3, w1 = 0, lambda = 2,
                      heterogeneous = TRUE)
    prefs <- make_preferences(s$pref_rank, id = ~ id, type = "ranking")
    ri_p(s, estimand = "mapns", assumption = "conditional",
         preferences = prefs, B = 199, min_informative = 0)
  })

  expect_size(p, 0.05)
  expect_size(p, 0.10)
})

test_that("power rises with the focal attribute's weight and with lambda", {
  skip_on_cran()
  set.seed(105)
  grid <- expand.grid(w1 = c(0.1, 0.5, 0.9), lambda = c(1, 5))
  rate <- mapply(function(w1, lambda) {
    mean(replicate(120, ri_p(
      sim_conjoint(N = 200, K = 6, L = 3, D = 2, w1 = w1, lambda = lambda),
      estimand = "apns", assumption = "separability", B = 199)) <= 0.05)
  }, grid$w1, grid$lambda)
  rate <- matrix(rate, nrow = 3, dimnames = list(c("0.1", "0.5", "0.9"),
                                                 c("lam1", "lam5")))

  # Monotone in w1 at each lambda ...
  expect_true(all(diff(rate[, "lam1"]) >= 0))
  expect_true(all(diff(rate[, "lam5"]) >= 0))
  # ... and in lambda at each w1 (sharper utility, more decisive swaps).
  expect_true(all(rate[, "lam5"] >= rate[, "lam1"]))
  # The strongest cell should be essentially always rejected.
  expect_gt(rate["0.9", "lam5"], 0.95)
})

test_that("heterogeneous orderings blunt the separability statistic while the
           conditional statistic recovers it", {
  skip_on_cran()
  # The paper's bottom row: opposing contributions cancel in the pooled
  # contrast, so the separability estimate stays near zero even when the focal
  # attribute drives the choice; conditioning on preference groups recovers it.
  # Cancellation is not exact in a finite sample -- the split between the two
  # orderings is only binomially balanced -- so the separability test still
  # rejects sometimes; what collapses is the magnitude, and with it the power
  # relative to the conditional test.
  set.seed(106)
  res <- replicate(150, {
    s <- sim_conjoint(N = 300, K = 10, L = 3, D = 2, w1 = 0.9, lambda = 5,
                      heterogeneous = TRUE)
    prefs <- make_preferences(s$pref_binary, id = ~ id, type = "binary")
    rs <- cj_ri_test(s$formula, data = s$data, id = ~ id, tasks = ~ task,
                     profile = ~ profile, attributes = "x1", estimand = "apns",
                     assumption = "separability", B = 199, checks = FALSE)
    rc <- cj_ri_test(s$formula, data = s$data, id = ~ id, tasks = ~ task,
                     profile = ~ profile, attributes = "x1", estimand = "apns",
                     assumption = "conditional", preferences = prefs,
                     B = 199, checks = FALSE, min_informative = 0)
    c(sep_stat = rs$results$statistic[1], sep_p = rs$results$p_value[1],
      cond_stat = rc$results$statistic[1], cond_p = rc$results$p_value[1])
  })

  # The conditional statistic recovers a large effect; the separability one
  # is an order of magnitude smaller.
  expect_gt(stats::median(res["cond_stat", ]), 0.5)
  expect_lt(stats::median(res["sep_stat", ]), 0.1)
  expect_gt(stats::median(res["cond_stat", ]) / stats::median(res["sep_stat", ]), 8)

  # Power follows: the conditional test always rejects, separability often
  # does not.
  expect_gt(mean(res["cond_p", ] <= 0.05), 0.95)
  expect_lt(mean(res["sep_p", ] <= 0.05), 0.75)
  expect_lt(mean(res["sep_p", ] <= 0.05), mean(res["cond_p", ] <= 0.05))
})
