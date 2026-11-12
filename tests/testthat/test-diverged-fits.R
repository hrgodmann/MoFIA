# Regression tests for estimates at which the likelihood is not finite within a
# finite-difference step (diverged optimization): a clear estimation error instead of
# "'x' must be an array of at least two dimensions" from the standard-error step.

ns <- asNamespace("MoFIA")
edge_model <- function() {
  set.seed(2); N <- 120; eta <- rnorm(N); x <- rnorm(N)
  d <- data.frame(y1 = eta + rnorm(N, sd = 0.6), y2 = 0.8 * eta + rnorm(N, sd = 0.6),
                  y3 = 1.2 * eta + rnorm(N, sd = 0.6), x = x)
  mofia_model(d, items = c("y1", "y2", "y3"), moderators = "x", continuous = c("y1", "y2", "y3"))
}
# The normal density's constant -log(2 pi theta)/2 overflows for log(theta) above this value.
overflow_edge <- log(.Machine$double.xmax / (2 * pi))

test_that("a likelihood that overflows within a finite-difference step gives a clear estimation error", {
  m <- edge_model()
  th <- get(".mofia_start", ns)(m)
  th["omega0[y2]"] <- overflow_edge - 0.5 * 1e-5 * abs(overflow_edge)   # half a step below the edge
  expect_true(is.finite(mofia_loglik(m, th)))
  # iter.max = 0: the optimizer returns the start, so the standard-error step runs at it
  err <- expect_error(mofia(m, start = th, control = list(iter.max = 0)), class = "mofia_estimation_error")
  expect_match(conditionMessage(err), "^Estimation failed: the log-likelihood is not finite within a finite-difference step")
  expect_match(conditionMessage(err), "omega0[y2]", fixed = TRUE)
  expect_match(conditionMessage(err), "the optimizer stopped with", fixed = TRUE)
  # se = FALSE still returns the unconverged fit, as before
  expect_warning(f <- mofia(m, start = th, control = list(iter.max = 0), se = FALSE), "may not have converged")
  expect_false(f$converged)
  expect_null(f$status)
  # the information is not computed at estimates with a non-finite likelihood either (direct call)
  good <- mofia(m)
  expect_true(good$converged)
  expect_identical(good$status, "ok")
  bad <- good; bad$coefficients["omega0[y2]"] <- 800
  expect_error(get(".mofia_se", ns)(bad), "not finite at the estimates", class = "mofia_estimation_error")
  # a non-finite likelihood at the solution itself: same message as before, now with the class
  th["omega0[y2]"] <- 800
  expect_error(mofia(m, start = th, control = list(iter.max = 0)),
               "Estimation failed: the log-likelihood is not finite at the solution.", fixed = TRUE,
               class = "mofia_estimation_error")
})

test_that("a factor correlation at the boundary of the parameter space gives the same clear error", {
  set.seed(4); N <- 150; f <- matrix(rnorm(2 * N), N) %*% chol(matrix(c(1, 0.5, 0.5, 1), 2))
  d <- as.data.frame(cbind(f[, 1] + matrix(rnorm(3 * N, sd = 0.6), N), f[, 2] + matrix(rnorm(3 * N, sd = 0.6), N)))
  names(d) <- paste0("y", 1:6)
  m <- mofia_model(d, loadings = list(f1 = paste0("y", 1:3), f2 = paste0("y", 4:6)), continuous = names(d))
  th <- get(".mofia_start", ns)(m)
  z <- atanh(1 - 1e-10)                                   # .structural() refuses |rho| >= 1 - 1e-10
  th["zcor0[f1,f2]"] <- z - 0.5 * 1e-5 * z                # half a finite-difference step inside
  expect_true(is.finite(mofia_loglik(m, th)))
  err <- expect_error(mofia(m, start = th, control = list(iter.max = 0)), class = "mofia_estimation_error")
  expect_match(conditionMessage(err), "within a finite-difference step of the estimates of zcor0[f1,f2]", fixed = TRUE)
})

test_that("score tests and factor scores handle a likelihood that is not finite", {
  m <- edge_model(); f <- mofia(m)
  cand <- c("omega[y2,x]", "beta[y1,x]")
  expect_false(anyNA(mofia_score_test(f, candidates = cand)$statistic))
  # a residual variance at the edge: the omega[y2,x] step overflows it for persons with x > 0.1
  q <- f; q$coefficients["omega0[y2]"] <- overflow_edge - 1e-6
  st <- mofia_score_test(q, candidates = cand)
  expect_true(is.na(st$statistic[st$test == "omega[y2,x]"]))
  expect_true(is.na(st$p[st$test == "omega[y2,x]"]))
  expect_false(is.na(st$statistic[st$test == "beta[y1,x]"]))
  q$coefficients["omega0[y2]"] <- 800
  expect_error(mofia_score_test(q, candidates = cand), "not finite at the estimates")
  expect_error(mofia_scores(q), "not finite at these estimates")
})

test_that("the three-factor Custom model (JASP fixture F07), not identified, stops on its ridge instead of diverging", {
  skip_on_cran()                                     # about 15 s
  file <- test_path("fixtures", "F07.csv")           # provenance: fixtures/F07-provenance.txt
  expect_identical(unname(tools::md5sum(file)), "501dae50676fa53f902cba7073da78b5")
  raw <- utils::read.csv(file)
  items <- c("y1", "y3", "y4", "y5", "y6", "y7", "y8", "y2")
  load <- list(Factor1 = c("y1", paste0("y", 3:8)), Factor2 = c("y2", paste0("y", 3:8)), Factor3 = c("y1", "y2"))
  m <- mofia_model(raw[, c(items, "Age")], items = items, loadings = load, moderators = ~ Age, standardize = TRUE,
                   continuous = items, loading_mod = character(0), threshold_mod = character(0),
                   resvar_mod = character(0), mean_mod = list(), var_mod = list(), cor_mod = FALSE,
                   identification = "variance", cov_param = "fisherz", sparse = 5)
  # Not identified (lavaan: rank-deficient information, -2LL 5515.166943 on the ridge). Before E2 the optimizer
  # followed the error of its fixed adaptive rule and diverged until a residual variance was at the limit of double
  # precision, and the standard-error step failed (E1: an estimation error naming omega0[y8] after 20 outer
  # iterations). With the exact posterior and the monotone outer loop it stops where the re-adapted log-likelihood
  # stops increasing: not a maximum, but with the exact log-likelihood.
  f <- suppressWarnings(mofia(m, link = "logit"))
  expect_false(f$converged)
  expect_true(f$status %in% c("saddle", "singular"))
  expect_s3_class(f$optimizer$monotone, "data.frame")
  expect_gt(-2 * f$loglik, 5515.166943 - 1e-3)
  expect_lt(abs(f$loglik - get(".cont_loglik", ns)(f$coefficients, m)), 1e-6)
  expect_lt(max(abs(f$coefficients[m$par$block == "omega0"])), 10)
})
