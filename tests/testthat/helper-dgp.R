# Simulated forced-choice conjoint matching the paper's DGP (Section
# "Recovering Attribute Relevance in Simulated Data"):
#
#   u_ijk = lambda_i * sum_l w_il * v_il(T_ijkl),   P(Y_ijk = 1) = softmax(u),
#
# with uniform, independent randomization of every attribute level in every
# profile, J = 2 profiles per task, and weights normalized to sum to 1.
#
# v_il maps the attribute's levels to equally spaced values in [0, 1] along
# respondent i's own ordering: homogeneous = everyone shares one ordering,
# heterogeneous = each respondent draws an independent random ordering.
#
# w1 = 0 makes attribute 1 irrelevant, so the sharp null of the randomization
# test holds exactly: utilities, and hence the choice distribution, do not
# depend on which profile carries which level of attribute 1.
sim_conjoint <- function(N = 300, K = 8, L = 3, D = 2, w1 = 1 / 3,
                         lambda = 2, heterogeneous = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  w <- c(w1, rep((1 - w1) / (L - 1), L - 1))
  levs <- as.character(seq_len(D))

  # Respondent-specific value functions: ord[[i]][[l]] is i's ordering of
  # attribute l's levels, worst to best, so v = (rank - 1) / (D - 1).
  ord <- lapply(seq_len(N), function(i)
    lapply(seq_len(L), function(l)
      if (heterogeneous) sample(levs) else levs))

  n_tasks <- N * K
  resp <- rep(seq_len(N), each = K * 2)
  task <- rep(rep(seq_len(K), each = 2), times = N)
  prof <- rep(c("a", "b"), times = n_tasks)

  X <- matrix(sample(levs, n_tasks * 2 * L, replace = TRUE), ncol = L)

  u <- numeric(nrow(X))
  for (l in seq_len(L)) {
    pos <- vapply(seq_along(u), function(r)
      match(X[r, l], ord[[resp[r]]][[l]]), integer(1))
    u <- u + w[l] * (pos - 1) / (D - 1)
  }
  u <- lambda * u

  ua <- u[seq(1, length(u), by = 2)]
  ub <- u[seq(2, length(u), by = 2)]
  pa <- exp(ua) / (exp(ua) + exp(ub))
  ya <- stats::rbinom(length(pa), 1L, pa)

  y <- numeric(length(u))
  y[seq(1, length(u), by = 2)] <- ya
  y[seq(2, length(u), by = 2)] <- 1L - ya

  dat <- data.frame(id = resp, task = task, profile = prof, y = y,
                    stringsAsFactors = FALSE)
  for (l in seq_len(L)) dat[[paste0("x", l)]] <- factor(X[, l], levels = levs)

  # Preference data: each respondent's most- and least-favored level per
  # attribute, measured without error.
  pref <- data.frame(id = seq_len(N), stringsAsFactors = FALSE)
  for (l in seq_len(L)) {
    pref[[paste0("x", l, "_top")]] <-
      vapply(seq_len(N), function(i) ord[[i]][[l]][D], character(1))
    pref[[paste0("x", l, "_bot")]] <-
      vapply(seq_len(N), function(i) ord[[i]][[l]][1], character(1))
  }
  # Binary form: 1 if the respondent prefers the attribute's second level.
  pref_bin <- data.frame(id = seq_len(N), stringsAsFactors = FALSE)
  for (l in seq_len(L))
    pref_bin[[paste0("x", l)]] <-
      as.integer(vapply(seq_len(N), function(i) ord[[i]][[l]][D],
                        character(1)) == levs[2])

  list(data = dat, pref_rank = pref, pref_binary = pref_bin,
       formula = stats::as.formula(paste("y ~", paste0("x", seq_len(L),
                                                       collapse = " + "))),
       w = w, lambda = lambda)
}

ri_p <- function(sim, ...) {
  rt <- cj_ri_test(sim$formula, data = sim$data, id = ~ id, tasks = ~ task,
                   profile = ~ profile, attributes = "x1", checks = FALSE, ...)
  rt$results$p_value[1]
}
