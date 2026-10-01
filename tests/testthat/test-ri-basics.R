test_that("the observed statistic equals the estimate cj_apns() reports", {
  sim <- sim_conjoint(N = 200, K = 6, L = 3, D = 3, w1 = 0.6, seed = 11)

  est <- cj_apns(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
                 profile = ~ profile, estimand = "mapns",
                 assumption = "separability", se = "none")
  rt <- cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
                   profile = ~ profile, estimand = "mapns",
                   assumption = "separability", B = 49, checks = FALSE)

  for (a in names(est$mapns))
    expect_equal(rt$results$statistic[rt$results$item == a],
                 est$mapns[[a]], tolerance = 1e-12)
})

test_that("the observed APNS statistic equals the pairwise cj_apns() estimate", {
  sim <- sim_conjoint(N = 150, K = 6, L = 2, D = 3, w1 = 0.7, seed = 12)

  est <- cj_apns(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
                 profile = ~ profile, estimand = "apns",
                 assumption = "separability", se = "none")
  rt <- cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
                   profile = ~ profile, estimand = "apns",
                   assumption = "separability", B = 49, checks = FALSE)

  for (i in seq_len(nrow(rt$results))) {
    a  <- rt$results$item[i]
    pr <- rt$results$pair[i]
    expect_equal(rt$results$statistic[i], est$apns[[a]][[pr]]$estimate,
                 tolerance = 1e-12)
  }
})

test_that("the conditional observed statistic equals cj_apns()'s", {
  sim <- sim_conjoint(N = 250, K = 8, L = 2, D = 3, w1 = 0.6,
                      heterogeneous = TRUE, seed = 13)
  prefs <- make_preferences(sim$pref_rank, id = ~ id, type = "ranking")

  est <- cj_apns(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
                 profile = ~ profile, estimand = "mapns",
                 assumption = "conditional", preferences = prefs, se = "none")
  rt <- cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
                   profile = ~ profile, estimand = "mapns",
                   assumption = "conditional", preferences = prefs,
                   B = 49, checks = FALSE, min_informative = 0)

  for (a in names(est$mapns))
    expect_equal(rt$results$statistic[rt$results$item == a],
                 est$mapns[[a]], tolerance = 1e-12)
})

test_that("with a single preference group the conditional and separability
           statistics coincide", {
  # Homogeneous orderings put every respondent in one group, where the
  # conditional estimator collapses to the separability estimator
  # (Propositions 1 and 2).
  sim <- sim_conjoint(N = 250, K = 8, L = 2, D = 2, w1 = 0.7,
                      heterogeneous = FALSE, seed = 14)
  prefs <- make_preferences(sim$pref_rank, id = ~ id, type = "ranking")

  rt <- cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
                   profile = ~ profile, estimand = "mapns",
                   assumption = "both", preferences = prefs,
                   B = 99, checks = FALSE)

  sep  <- rt$results[rt$results$assumption == "separability", ]
  cond <- rt$results[rt$results$assumption == "conditional", ]
  expect_equal(sep$statistic, cond$statistic[match(sep$item, cond$item)],
               tolerance = 1e-12)
  expect_equal(sep$p_value, cond$p_value[match(sep$item, cond$item)])
})

test_that("seeded runs reproduce exactly and leave the caller's RNG alone", {
  sim <- sim_conjoint(N = 120, K = 5, L = 2, D = 2, w1 = 0.5, seed = 15)
  args <- list(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
               profile = ~ profile, B = 99, checks = FALSE, return_draws = TRUE)

  a <- do.call(cj_ri_test, c(args, list(seed = 7)))
  b <- do.call(cj_ri_test, c(args, list(seed = 7)))
  expect_identical(a$results, b$results)
  expect_identical(a$draws, b$draws)

  c3 <- do.call(cj_ri_test, c(args, list(seed = 8)))
  expect_false(isTRUE(all.equal(a$draws[[1]], c3$draws[[1]])))

  set.seed(99); before <- stats::runif(1)
  set.seed(99); invisible(do.call(cj_ri_test, c(args, list(seed = 7))))
  expect_equal(stats::runif(1), before)
})

test_that("p-values are in the valid range and respect the (B+1) grid", {
  sim <- sim_conjoint(N = 150, K = 6, L = 2, D = 2, w1 = 0.6, seed = 16)
  rt <- cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
                   profile = ~ profile, B = 199, checks = FALSE)
  p <- rt$results$p_value
  expect_true(all(p >= 1 / 200 & p <= 1))
  expect_equal(p * 200, round(p * 200), tolerance = 1e-10)
})

test_that("inputs are validated", {
  sim <- sim_conjoint(N = 60, K = 4, L = 2, D = 2, seed = 17)

  expect_error(
    cj_ri_test(sim$formula, data = sim$data, id = ~ id, profile = ~ profile),
    "'tasks' must be specified")

  expect_error(
    cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
               profile = ~ profile, assumption = "conditional"),
    "'preferences' must be provided")

  expect_error(
    cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
               profile = ~ profile, attributes = "nope"),
    "not in 'formula'")

  # J > 2 is out of scope
  three <- rbind(sim$data, transform(sim$data[sim$data$profile == "a", ],
                                     profile = "c"))
  expect_error(
    cj_ri_test(sim$formula, data = three, id = ~ id, tasks = ~ task,
               profile = ~ profile, B = 9, checks = FALSE),
    "two-profile designs only")
})

test_that("return_draws returns B draws whose tail reproduces the p-value", {
  sim <- sim_conjoint(N = 120, K = 5, L = 2, D = 2, w1 = 0.5, seed = 18)
  rt <- cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
                   profile = ~ profile, attributes = "x1", B = 199,
                   checks = FALSE, return_draws = TRUE)
  d <- rt$draws[[1]]
  expect_length(d, 199)
  expect_equal((1 + sum(d >= rt$results$statistic[1] - 1e-10)) / 200,
               rt$results$p_value[1])
})
