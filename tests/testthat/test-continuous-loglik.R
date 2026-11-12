# Regression tests for continuous items with a residual variance near zero: the quadrature
# must agree with the closed-form normal log-likelihood, a disagreement at the solution is an
# estimation error, and the monotone outer loop does not let the optimizer follow the error
# of a fixed adaptive rule.

ns <- asNamespace("MoFIA")
# Independent closed form: y_i ~ N(nu_i + Lambda_i mu_i, Lambda_i C_i C_i' Lambda_i' + Theta_i), per person with chol().
mvn_loglik <- function(theta, m) {
  ei <- get(".engine_inputs", ns)(theta, m); P <- m$P; D <- m$D; lt <- get(".lt_index", ns)
  sum(vapply(seq_len(m$N), function(i) {
    C <- matrix(0, D, D); for (d in seq_len(D)) for (e in seq_len(d)) C[d, e] <- ei$CL[i, lt(d, e)]
    L <- matrix(ei$LAM[i, ], P, D); o <- !is.na(m$Y[i, ])
    S <- (tcrossprod(L %*% C) + diag(ei$RV[i, ], P))[o, o, drop = FALSE]
    R <- chol(S); z <- backsolve(R, (m$Y[i, ] - ei$NU[i, ] - drop(L %*% ei$st$mu[i, ]))[o], transpose = TRUE)
    -0.5 * (sum(o) * log(2 * pi) + sum(z^2)) - sum(log(diag(R)))
  }, 1))
}
# One factor (or two), y1 without residual variance in the population; seed 8 gives a negative unconstrained ML
# estimate of its residual variance (lavaan: -0.017), so the ML residual variance in MoFIA's space is 0.
heywood_data <- function(seed = 8, N = 300, two = FALSE) {
  set.seed(seed); age <- rnorm(N)
  if (!two) {
    eta <- 0.3 * age + exp(0.15 * age) * rnorm(N); lam <- c(1, 0.8, 0.7, 0.9, 0.6, 0.75); th <- c(0, 0.5, 0.5, 0.5, 0.6, 0.5)
    Y <- sapply(1:6, function(j) 0.2 * j + lam[j] * eta + rnorm(N, sd = sqrt(th[j])))
  } else {
    f <- matrix(rnorm(2 * N), N) %*% chol(matrix(c(1, 0.4, 0.4, 1), 2)); f[, 1] <- f[, 1] + 0.3 * age
    lam <- c(1, 0.8, 0.7, 0.9, 1, 0.8, 0.7, 0.9); th <- c(0, rep(0.5, 7))
    Y <- sapply(1:8, function(j) 0.2 * j + lam[j] * f[, if (j <= 4) 1 else 2] + rnorm(N, sd = sqrt(th[j])))
  }
  d <- as.data.frame(Y); names(d) <- paste0("y", seq_len(ncol(Y))); d$Age <- 40 + 12 * age; d
}
heywood_model <- function(d, loadings, resvar_mod = character(0), continuous = setdiff(names(d), "Age"))
  mofia_model(d, items = setdiff(names(d), "Age"), loadings = loadings, moderators = ~ Age, standardize = TRUE,
              continuous = continuous, identification = "variance", cov_param = "fisherz", sparse = 5,
              loading_mod = character(0), threshold_mod = character(0), resvar_mod = resvar_mod,
              mean_mod = stats::setNames(as.list(rep(TRUE, length(loadings))), names(loadings)),
              var_mod = stats::setNames(as.list(rep(TRUE, length(loadings))), names(loadings)), cor_mod = length(loadings) == 2)

test_that("the adaptive quadrature stays exact as a continuous residual variance approaches zero", {
  i6 <- paste0("y", 1:6); i8 <- paste0("y", 1:8)
  for (m in list(heywood_model(heywood_data(N = 120), list(F = i6)),
                 heywood_model(heywood_data(3, N = 120, two = TRUE), list(F1 = i8[1:4], F2 = i8[5:8])))) {
    th <- get(".mofia_start", ns)(m)
    for (v in c(1e-2, 1e-6, 1e-9, 1e-10, 1e-14, 1e-16)) {
      th["omega0[y1]"] <- log(v)
      expect_lt(abs(mofia_loglik(m, th) - mvn_loglik(th, m)), 1e-6)
    }
  }
})

test_that("the package's closed form agrees with an independent evaluation, with missing values", {
  d <- heywood_data(N = 80); d$y3[c(2, 5, 9)] <- NA; d$y1[c(5, 30)] <- NA
  m <- heywood_model(d, list(F = paste0("y", 1:6)), resvar_mod = "y2")
  th <- get(".mofia_start", ns)(m); th["omega[y2,Age]"] <- 0.3; th["gamma[F,Age]"] <- 0.2
  expect_equal(get(".cont_loglik", ns)(th, m), mvn_loglik(th, m), tolerance = 1e-12)
  expect_equal(sum(get(".cont_loglik", ns)(th, m, per_person = TRUE)), mvn_loglik(th, m), tolerance = 1e-12)
  expect_equal(mofia_loglik(m, th), mvn_loglik(th, m), tolerance = 1e-12)
  expect_error(get(".cont_loglik", ns)(th, mofia_model(parabelief, items = c("SpiritsGhosts", "Aliens", "Sorcery"))),
               "only continuous items")
})

test_that("a quadrature log-likelihood that disagrees with the closed form at the solution is an estimation error", {
  m <- heywood_model(heywood_data(N = 120), list(F = paste0("y", 1:6)))
  f <- suppressWarnings(mofia(m, se = FALSE))
  expect_equal(f$loglik, mvn_loglik(f$coefficients, m), tolerance = 1e-6 / abs(f$loglik))
  local_mocked_bindings(.cont_loglik = function(theta, spec, per_person = FALSE) mvn_loglik(theta, spec) + 2e-6)
  err <- expect_error(mofia(m), class = "mofia_estimation_error")
  expect_match(conditionMessage(err), "^Estimation failed: at the solution the log-likelihood from quadrature")
  expect_match(conditionMessage(err), "The smallest residual variance is [^ ]+ \\(item y[0-9]\\)")
  local_mocked_bindings(.cont_loglik = function(theta, spec, per_person = FALSE) NaN)
  expect_error(mofia(m, se = FALSE), "differs from its exact value \\(not computable\\)", class = "mofia_estimation_error")
})

# The monotone outer loop on a synthetic likelihood: true log-likelihood -sum((theta - 1)^2); a fixed
# "adaptation" centred at the point where it was made adds a bonus that the optimizer follows. With binary items
# (one threshold each, so the raw and natural parameters coincide) the loop only halves, undoes and retries; with
# continuous items it also switches to re-adapting at every evaluation once it has acted.
synthetic_loop <- function(bonus, grad_bonus, start, continuous = FALSE) {
  d <- heywood_data(N = 30)[, c("y1", "y2", "y3")]
  db <- as.data.frame(lapply(d, function(v) as.integer(v > stats::median(v))))
  m <- if (continuous) mofia_model(d, items = names(d), continuous = names(d)) else mofia_model(db, items = names(db))
  local_mocked_bindings(
    .mofia_adapt = function(theta, spec, link = "logit", nodes = NULL, adapt = NULL, ...) list(centre = unname(theta)),
    .mofia_ll = function(theta, spec, link = "logit", nodes = NULL, scores = FALSE, adapt = NULL, ...) {
      th <- unname(theta); out <- -sum((th - 1)^2) + bonus(th, adapt$centre)
      if (scores) attr(out, "scores") <- matrix(-2 * (th - 1) + grad_bonus(th, adapt$centre), 1)
      out
    })
  st <- stats::setNames(rep(start, nrow(m$par)), m$par$name)
  get(".mofia_optimize", ns)(m, "logit", list(nq = 5), st, list(eval.max = 3000, iter.max = 2000, rel.tol = 1e-10), FALSE, 1e-6)
}

test_that("the monotone outer loop halves an outer iteration that lowers the re-adapted log-likelihood", {
  # the fixed rule rewards moving away from its centre: nlminb overshoots the maximum
  r <- synthetic_loop(function(th, c) 0.6 * sum((th - c)^2), function(th, c) 1.2 * (th - c), start = 0)
  expect_s3_class(r$monotone, "data.frame")
  expect_true(all(r$monotone$loglik_end < r$monotone$loglik_start - 1e-6))
  expect_true(all(r$monotone$halvings >= 1))
  expect_true(all(r$monotone$loglik >= r$monotone$loglik_start - 1e-6))
  expect_equal(unname(r$theta), rep(1, length(r$theta)), tolerance = 1e-3)
  expect_equal(r$loglik, -sum((unname(r$theta) - 1)^2))                       # re-adapted, not the stale value
  expect_identical(r$monotone$nq[1], 5)
})

test_that("the monotone outer loop undoes and retries an iteration that no halving improves, and leaves exact fits alone", {
  # at the maximum, a tilted rule moves the optimizer away; every shorter step is worse too, and so are the retries
  # with a tenth of the iterations
  r <- synthetic_loop(function(th, c) 4 * sum(th - c), function(th, c) rep(4, length(th)), start = 1)
  expect_identical(nrow(r$monotone), 4L)
  expect_true(all(is.na(r$monotone$halvings)))
  expect_equal(r$monotone$iter_max, c(2000, 200, 20, 2))
  expect_identical(r$outer, 4L)
  expect_equal(unname(r$theta), rep(1, length(r$theta)))
  expect_equal(r$loglik, 0)
  # an exact rule: no action, nothing recorded
  r <- synthetic_loop(function(th, c) 0, function(th, c) 0 * th, start = 0)
  expect_null(r$monotone)
  expect_identical(r$exact_outer, 0L)
  expect_equal(unname(r$theta), rep(1, length(r$theta)), tolerance = 1e-6)
})

test_that("with continuous items the loop maximizes the exact log-likelihood once it has acted", {
  # the same tilted rule: the first iteration is undone, the next ones re-adapt at every evaluation and stay at the maximum
  r <- synthetic_loop(function(th, c) 4 * sum(th - c), function(th, c) rep(4, length(th)), start = 1, continuous = TRUE)
  expect_identical(nrow(r$monotone), 1L)
  expect_gt(r$exact_outer, 0L)
  expect_equal(unname(r$theta), rep(1, length(r$theta)), tolerance = 1e-6)
  expect_equal(r$loglik, 0, tolerance = 1e-10)
  # an overshooting rule: one halving, then the exact likelihood is maximized
  r <- synthetic_loop(function(th, c) 0.6 * sum((th - c)^2), function(th, c) 1.2 * (th - c), start = 0, continuous = TRUE)
  expect_identical(nrow(r$monotone), 1L)
  expect_gt(r$exact_outer, 0L)
  expect_equal(unname(r$theta), rep(1, length(r$theta)), tolerance = 1e-6)
})

test_that("an identified Heywood case reaches Mplus's -2LL with the exact log-likelihood (D = 1 and D = 2)", {
  skip_on_cran()                                       # about 5 s
  i6 <- paste0("y", 1:6); i8 <- paste0("y", 1:8)
  # residual variances not moderated: the ML is at the boundary RV[y1] = 0; Mplus (free, and with y1@0): -2LL 4165.460
  # (D = 1), 5897.228 (D = 2); before E2 MoFIA reported 2707.51 and 2653.05
  for (cs in list(list(m = heywood_model(heywood_data(), list(F = i6)), sup = 4165.460),
                  list(m = heywood_model(heywood_data(3, two = TRUE), list(F1 = i8[1:4], F2 = i8[5:8])), sup = 5897.228))) {
    f <- suppressWarnings(mofia(cs$m))
    expect_lt(abs(f$loglik - mvn_loglik(f$coefficients, cs$m)), 1e-6)
    expect_lt(abs(-2 * f$loglik - cs$sup), 2e-3)
    expect_lt(exp(f$coefficients[["omega0[y1]"]]), 1e-6)
    expect_gt(f$optimizer$exact_outer, 0L)
  }
  # residual variances moderated: an interior maximum, which Mplus also reaches (-2LL 4151.440; before E2: 3692.12)
  m <- heywood_model(heywood_data(), list(F = i6), resvar_mod = i6)
  f <- mofia(m)
  expect_identical(f$status, "ok")
  expect_true(f$converged)
  expect_lt(abs(-2 * f$loglik - 4151.440), 2e-3)
  expect_lt(abs(f$loglik - mvn_loglik(f$coefficients, m)), 1e-6)
})

test_that("JASP fixture F07c: Scalar and Strict report the exact log-likelihood, nested, or fail clearly", {
  skip_on_cran()                                       # about 40 s
  file <- test_path("fixtures", "F07.csv")             # provenance: fixtures/F07-provenance.txt
  raw <- utils::read.csv(file)
  items <- c("y1", "y3", "y4", "y5", "y6", "y7", "y8", "y2")
  load <- list(Factor1 = c("y1", paste0("y", 3:8)), Factor2 = c("y2", paste0("y", 3:8)), Factor3 = c("y1", "y2"))
  all <- list(Factor1 = TRUE, Factor2 = TRUE, Factor3 = TRUE); free <- setdiff(items, c("y1", "y2"))
  pre <- list(configural = list(loading_mod = free, threshold_mod = free, resvar_mod = items),
              scalar = list(loading_mod = character(0), threshold_mod = character(0), resvar_mod = items),
              strict = list(loading_mod = character(0), threshold_mod = character(0), resvar_mod = character(0)))
  m2LL <- c()
  for (p in names(pre)) {
    m <- do.call(mofia_model, c(list(data = raw[, c(items, "Age")], items = items, loadings = load, moderators = ~ Age,
                                     standardize = TRUE, continuous = items, mean_mod = all, var_mod = all, cor_mod = FALSE,
                                     identification = "variance", cov_param = "fisherz", sparse = 5), pre[[p]]))
    f <- tryCatch(suppressWarnings(mofia(m)), mofia_estimation_error = function(e) e)
    if (inherits(f, "error")) {
      # Strict: its maximum lies on the boundary of positive definiteness (as in Mplus, r13 = 1), so the
      # information cannot be computed (E1's error); before E2 a saddle with -2LL 5201.92 (exact: 5558.94)
      expect_identical(p, "strict")
      expect_match(conditionMessage(f), "within a finite-difference step of the estimates of zcor0", fixed = TRUE)
      next
    }
    expect_lt(abs(f$loglik - mvn_loglik(f$coefficients, m)), 1e-6)       # before E2: Scalar 4139.78 vs 5517.21
    if (p != "configural") { expect_false(f$converged); expect_true(f$status %in% c("saddle", "singular")) }
    m2LL[p] <- -2 * f$loglik
  }
  expect_lte(m2LL[["configural"]], m2LL[["scalar"]])
  if (!is.na(m2LL["strict"])) expect_lte(m2LL[["scalar"]], m2LL[["strict"]])
})

test_that("a mixed model pushed to a near-zero continuous residual variance is not silent", {
  skip_on_cran()                                       # about 5 s
  d <- heywood_data(N = 200)
  for (v in c("y5", "y6")) d[[v]] <- as.integer(cut(d[[v]], stats::quantile(d[[v]], c(0, 1/3, 2/3, 1)), include.lowest = TRUE))
  m <- heywood_model(d, list(F = paste0("y", 1:6)), continuous = paste0("y", 1:4))
  th <- get(".mofia_start", ns)(m); th["omega0[y1]"] <- log(1e-14)
  w <- character(0)
  f <- withCallingHandlers(tryCatch(mofia(m, start = th, control = list(iter.max = 0)), error = function(e) e),
                           warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_true(inherits(f, "error") || any(grepl("quadrature", w, ignore.case = TRUE)))
})
