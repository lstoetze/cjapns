# Proposition "Pooled Estimation under Full Randomization": under uniform,
# independent randomization the pooled difference in means targets the same
# estimand as the context-stratified estimator.

# Context-stratified pairwise APNS, computed directly from the long data:
# within each realized context t_[-l] -- here the remaining-attribute vectors
# of the two profiles, labelled by which one carries t_q -- take the
# difference in means, then average across contexts with a common weight.
strat_apns <- function(d, outcome, attribute, tq, tp, others, weights = NULL) {
  key <- paste(d$id, d$task)
  o <- order(key, d$profile)
  d <- d[o, ]; key <- key[o]
  i1 <- seq(1, nrow(d), by = 2); i2 <- i1 + 1

  l1 <- as.character(d[[attribute]])[i1]
  l2 <- as.character(d[[attribute]])[i2]
  keep <- (l1 == tq & l2 == tp) | (l1 == tp & l2 == tq)
  i1 <- i1[keep]; i2 <- i2[keep]; l1 <- l1[keep]

  ctx_of <- function(idx)
    do.call(paste, c(lapply(others, function(v) as.character(d[[v]])[idx]),
                     sep = "\r"))
  # the context records which remaining-attribute vector sits on the
  # t_q-carrying profile
  ctx_q <- ifelse(l1 == tq, ctx_of(i1), ctx_of(i2))
  ctx_p <- ifelse(l1 == tq, ctx_of(i2), ctx_of(i1))
  ctx   <- paste(ctx_q, ctx_p, sep = "\n")

  y <- as.numeric(d[[outcome]])
  y_q <- ifelse(l1 == tq, y[i1], y[i2])
  y_p <- ifelse(l1 == tq, y[i2], y[i1])

  cs <- split(seq_along(ctx), ctx)
  delta <- vapply(cs, function(j) mean(y_q[j]) - mean(y_p[j]), numeric(1))
  p_hat <- if (is.null(weights)) lengths(cs) / length(ctx) else weights[names(cs)]
  abs(sum(delta * p_hat))
}

test_that("the pooled estimator is the context-stratified estimator at the
           empirical context weights", {
  # With a common p_hat equal to the empirical context distribution among
  # informative tasks, the two coincide algebraically, not just in the limit:
  # each informative task contributes exactly one t_q profile and one t_p
  # profile, so both pooled means carry the same context weights.
  sim <- sim_conjoint(N = 250, K = 8, L = 3, D = 2, w1 = 0.5, seed = 301)
  pooled <- abs(estimate_pairwise(sim$data, "y", "x1", "1", "2",
                                  id_var = "id", task_var = "task",
                                  informative = "informative",
                                  profile_var = "profile")$estimate)
  expect_equal(pooled,
               strat_apns(sim$data, "y", "x1", "1", "2", c("x2", "x3")),
               tolerance = 1e-12)
})

test_that("stratifying at the design's context weights agrees with the pooled
           estimator under uniform independent randomization, and the gap
           shrinks with sample size", {
  skip_on_cran()
  # Here the common weight is the *design* distribution (uniform over the
  # realized contexts) rather than the empirical one, so the two estimators
  # differ in finite samples and agree only in the limit.
  gap <- function(N, seed) {
    sim <- sim_conjoint(N = N, K = 8, L = 3, D = 2, w1 = 0.5, seed = seed)
    pooled <- abs(estimate_pairwise(sim$data, "y", "x1", "1", "2",
                                    id_var = "id", task_var = "task",
                                    informative = "informative",
                                    profile_var = "profile")$estimate)
    # 2 other attributes x 2 levels, on each of two profiles: 16 contexts,
    # equally likely by design.
    ctxs <- apply(expand.grid(rep(list(c("1", "2")), 4)), 1,
                  function(v) paste(paste(v[1], v[2], sep = "\r"),
                                    paste(v[3], v[4], sep = "\r"), sep = "\n"))
    w <- stats::setNames(rep(1 / length(ctxs), length(ctxs)), ctxs)
    abs(pooled - strat_apns(sim$data, "y", "x1", "1", "2", c("x2", "x3"),
                            weights = w))
  }

  small <- mean(vapply(1:20, function(s) gap(100,  300 + s), numeric(1)))
  large <- mean(vapply(1:20, function(s) gap(1600, 400 + s), numeric(1)))
  expect_lt(large, small)
  expect_lt(large, 0.02)
})

test_that("declaring the true uniform marginals leaves the estimate
           essentially unchanged", {
  sim <- sim_conjoint(N = 600, K = 8, L = 3, D = 2, w1 = 0.5, seed = 310)
  des <- make_design(level_probs = list(x2 = c("1" = 0.5, "2" = 0.5),
                                        x3 = c("1" = 0.5, "2" = 0.5)))

  a <- cj_apns(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
               profile = ~ profile, estimand = "mapns",
               assumption = "separability", se = "none")
  b <- cj_apns(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
               profile = ~ profile, estimand = "mapns",
               assumption = "separability", se = "none", design = des)

  expect_equal(unlist(a$mapns), unlist(b$mapns), tolerance = 0.02)
})
