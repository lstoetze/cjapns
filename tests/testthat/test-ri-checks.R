test_that("orientation imbalance is detected", {
  sim <- sim_conjoint(N = 200, K = 6, L = 3, D = 2, w1 = 0.4, seed = 201)
  d <- sim$data

  # Force profile "a" to carry level "2" in every informative task, which is
  # exactly the position-dependent randomizer e_ik != 1/2 is meant to catch.
  key <- paste(d$id, d$task)
  ia <- which(d$profile == "a"); ib <- which(d$profile == "b")
  ia <- ia[order(key[ia])];      ib <- ib[order(key[ib])]
  inf <- d$x1[ia] != d$x1[ib]
  d$x1[ia[inf]] <- factor("2", levels = levels(d$x1))
  d$x1[ib[inf]] <- factor("1", levels = levels(d$x1))

  expect_warning(
    cj_ri_test(sim$formula, data = d, id = ~ id, tasks = ~ task,
               profile = ~ profile, attributes = "x1", B = 49, checks = TRUE),
    "Orientation imbalance")

  rt <- suppressWarnings(
    cj_ri_test(sim$formula, data = d, id = ~ id, tasks = ~ task,
               profile = ~ profile, attributes = "x1", B = 49, checks = TRUE))
  # Profile 1 always carries level "2", i.e. t_p, so D is 0 throughout.
  bal <- rt$checks$x1$balance
  expect_equal(bal$n_D1, 0L)
  expect_lt(bal$p_value, 1e-10)
})

test_that("a balanced design raises no balance warning", {
  sim <- sim_conjoint(N = 200, K = 6, L = 3, D = 2, w1 = 0.4, seed = 202)
  expect_no_warning(
    cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
               profile = ~ profile, attributes = "x1", B = 49, checks = TRUE))
})

test_that("dependence between the focal attribute and another is detected", {
  sim <- sim_conjoint(N = 200, K = 6, L = 2, D = 2, w1 = 0.5, seed = 203)
  d <- sim$data
  # x2 copies x1 90% of the time: the attribute is no longer assigned
  # independently of the others, so e_ik is not 1/2.
  copy <- stats::runif(nrow(d)) < 0.9
  d$x2[copy] <- d$x1[copy]

  expect_warning(
    cj_ri_test(sim$formula, data = d, id = ~ id, tasks = ~ task,
               profile = ~ profile, attributes = "x1", B = 49, checks = TRUE),
    "not independent across attributes")
})

test_that("a restricted design shows up as more missing swaps than sparsity
           explains", {
  sim <- sim_conjoint(N = 300, K = 8, L = 2, D = 2, w1 = 0.5, seed = 204)
  d <- sim$data
  # Design restriction: a profile showing x1 = "2" must show x2 = "2". Those
  # informative configurations then have no feasible swapped counterpart.
  d$x2[d$x1 == "2"] <- factor("2", levels = levels(d$x2))

  # The restriction also makes x1 and x2 dependent, so both checks fire.
  w <- character(0)
  rt <- withCallingHandlers(
    cj_ri_test(sim$formula, data = d, id = ~ id, tasks = ~ task,
               profile = ~ profile, attributes = "x1", B = 49, checks = TRUE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("than sparsity alone explains", w)))
  expect_true(any(grepl("not independent across attributes", w)))

  sw <- rt$checks$x1$swaps
  expect_gt(sw$share_infeasible, sw$share_expected)
  expect_lt(sw$p_value, 0.01)
})

test_that("sparse contexts do not trigger the swap warning", {
  # Five other attributes make almost every configuration unique, so nearly
  # every swap is unobserved -- but that is sparsity, not a restriction, and
  # the calibration against the orientation re-draw has to say so.
  sim <- sim_conjoint(N = 200, K = 6, L = 6, D = 3, w1 = 0.2, seed = 205)
  rt <- expect_no_warning(
    cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
               profile = ~ profile, attributes = "x1", B = 49, checks = TRUE))
  sw <- rt$checks$x1$swaps
  expect_gt(sw$share_infeasible, 0.5)       # raw share is high ...
  expect_equal(sw$share_infeasible, sw$share_expected, tolerance = 0.05)
})

test_that("the swap check abstains when sparsity leaves nothing to detect", {
  # Eight other attributes make every configuration unique, so the expected
  # share is pinned near 1 and its reference distribution collapses onto a
  # point: a difference of no practical size would otherwise land in the tail.
  sim <- sim_conjoint(N = 200, K = 8, L = 9, D = 4, w1 = 0.1, seed = 208)
  rt <- expect_no_warning(
    cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
               profile = ~ profile, attributes = "x1", B = 49, checks = TRUE))
  sw <- rt$checks$x1$swaps
  expect_gt(sw$share_expected, 0.9)
  expect_false(sw$informative)
})

test_that("thin preference groups are warned about", {
  sim <- sim_conjoint(N = 120, K = 3, L = 2, D = 4, w1 = 0.5,
                      heterogeneous = TRUE, seed = 206)
  prefs <- make_preferences(sim$pref_rank, id = ~ id, type = "ranking")
  expect_warning(
    cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
               profile = ~ profile, attributes = "x1", estimand = "mapns",
               assumption = "conditional", preferences = prefs,
               B = 49, checks = FALSE, min_informative = 10),
    "fewer than 10 informative tasks")
})

test_that("a conditional MAPNS that is not identified is skipped, not faked", {
  sim <- sim_conjoint(N = 200, K = 6, L = 2, D = 3, w1 = 0.5,
                      heterogeneous = TRUE, seed = 207)
  # Binary preferences do not carry each group's own extreme pair, so the
  # conditional MAPNS is unidentified for a 3-level attribute.
  prefs <- make_preferences(sim$pref_binary, id = ~ id, type = "binary")
  rt <- cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
                   profile = ~ profile, estimand = "mapns",
                   assumption = "conditional", preferences = prefs,
                   B = 49, checks = FALSE)

  expect_equal(nrow(rt$results), 0)
  expect_true(all(grepl("not identified", rt$skipped$reason)))
})
