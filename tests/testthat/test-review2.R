# Regression tests for the findings of the second external review (7 November 2026).

ns <- asNamespace("MoFIA")

cont_fixture <- function() {
  set.seed(1047); N <- 400; x <- rep(c(-1, 1), N / 2); eta <- rnorm(N, 0.2 * x, exp(0.1 * x))
  d <- data.frame(y1 = eta + rnorm(N), y2 = 1.2 * eta + rnorm(N), y3 = 0.8 * eta + rnorm(N),
                  y4 = 1.1 * eta + rnorm(N), x = x)
  list(d = d, m = mofia_model(d, items = paste0("y", 1:4), moderators = "x", continuous = paste0("y", 1:4)))
}

test_that("F1: score screening works with a single binary moderator; bundling equals separate tests", {
  skip_on_cran()
  fx <- cont_fixture(); f0 <- mofia(fx$m)
  sc <- mofia_score_test(f0)
  expect_true(nrow(sc) >= 12)
  expect_true(all(is.finite(sc$statistic)))
  b <- mofia_score_test(f0, paste0("beta[y", 1:4, ",x]"))
  one <- mofia_score_test(f0, "beta[y2,x]")
  expect_equal(b$statistic[b$test == "beta[y2,x]"], one$statistic, tolerance = 1e-10)
  expect_equal(one$statistic, 2.561826, tolerance = 1e-5)
})

test_that("F2: Wald and score tests refuse unconverged fits", {
  skip_on_cran()
  fx <- cont_fixture()
  fu <- suppressWarnings(mofia(fx$m, control = list(iter.max = 0, eval.max = 1)))
  expect_false(fu$converged)
  expect_error(mofia_wald(fu, "gamma[eta,x]"), "not converged")
  expect_error(mofia_score_test(fu, "beta[y2,x]"), "not converged")
  expect_output(print(summary(fu)), "has not converged")
  expect_warning(confint(fu), "not converged")
})

test_that("F3: models with only continuous items integrate exactly, also for concentrated posteriors", {
  set.seed(7310)
  for (D in 1:3) for (sd in c(0.1, 0.003)) {
    n <- 20; eta <- matrix(rnorm(n * D), n, D)
    d <- as.data.frame(sapply(1:(3 * D), function(j) eta[, ceiling(j / 3)] + rnorm(n, sd = sd)))
    its <- paste0("y", seq_len(3 * D)); names(d) <- its
    ld <- stats::setNames(split(its, rep(1:D, each = 3)), paste0("f", 1:D))
    m <- mofia_model(d, loadings = ld, continuous = its)
    th <- get(".mofia_start", ns)(m)
    th[m$par$block == "nu0"] <- 0; th[m$par$block == "lambda0"] <- 1
    th[m$par$block == "omega0"] <- 2 * log(sd); th[m$par$block == "r0"] <- 0
    L <- kronecker(diag(D), matrix(1, 3, 1)); C <- chol(L %*% t(L) + diag(sd^2, 3 * D))
    ex <- sum(apply(d, 1, function(y) -0.5 * (3 * D * log(2 * pi) + 2 * sum(log(diag(C))) +
                                                sum(backsolve(C, y, transpose = TRUE)^2))))
    expect_equal(mofia_loglik(m, th), ex, tolerance = 1e-9)
  }
})

test_that("F3: the quadrature check refines the rule for very discriminating items", {
  skip_on_cran()
  set.seed(771); n <- 400; eta <- rnorm(n)
  d <- as.data.frame(sapply(1:4, function(j) as.integer(4 * eta + rnorm(n) > 0))); names(d) <- paste0("y", 1:4)
  m <- mofia_model(d, items = names(d), sparse = 0)
  f <- mofia(m, link = "probit")
  expect_gt(f$nq, 31)
  expect_true(nrow(f$quadrature) >= 2)
  exact <- function(t) {
    lp <- t[m$par$block == "lambda0"]; tau <- t[m$par$block == "tau"]
    sum(apply(m$Y, 1, function(y) log(integrate(function(z) {
      p <- rep(1, length(z)); for (j in 1:4) p <- p * pnorm(tau[j] - lp[j] * z, lower.tail = y[j] == 0); p * dnorm(z)
    }, -Inf, Inf, rel.tol = 1e-10)$value)))
  }
  expect_lt(abs(f$loglik - exact(coef(f))), 1e-3)
  expect_output(print(f), "check: logLik changes")
})

test_that("F4: marginal predictions are exact (probit) or accurate (logit)", {
  data(parabelief)
  m <- mofia_model(parabelief, items = c("SpiritsGhosts", "Aliens", "Sorcery"), moderators = "Male", sparse = 0)
  th <- get(".mofia_start", ns)(m, "probit")
  th[m$par$block == "lambda0"] <- c(2.5, 3, 4); th["delta[eta,Male]"] <- 2; th["gamma[eta,Male]"] <- 0.4
  fake <- structure(list(coefficients = th, model = m, link = "probit", nq = 31), class = "mofia")
  p <- predict(fake, newdata = data.frame(Male = c(0, 1)), type = "prob")$Aliens
  tau <- th["tau[Aliens,1]"]; lam <- 3; mu <- c(0, 0.4); v <- c(1, exp(2))
  expect_equal(unname(1 - p[, 1]), unname(pnorm((lam * mu - tau) / sqrt(1 + lam^2 * v))), tolerance = 1e-12)
  fake$link <- "logit"
  p <- predict(fake, newdata = data.frame(Male = c(0, 1)), type = "prob")$Aliens
  ref <- vapply(1:2, function(i) integrate(function(z) plogis(lam * mu[i] - tau + lam * sqrt(v[i]) * z) * dnorm(z),
                                           -Inf, Inf, rel.tol = 1e-12)$value, 1)
  expect_equal(unname(1 - p[, 1]), ref, tolerance = 1e-6)
  expect_equal(unname(rowSums(p)), c(1, 1), tolerance = 1e-12)
})

test_that("F5: the SB2010 test refuses models not nested by parameter names", {
  skip_on_cran()
  fx <- cont_fixture(); its <- paste0("y", 1:4)
  f0 <- mofia(fx$m)
  fm <- mofia(mofia_model(fx$d, items = its, moderators = "x", continuous = its, identification = "marker",
                          threshold_mod = "y2"))
  expect_true(is.finite(mofia_lrt(f0, fm, check_nested = FALSE)$p_value))
  expect_error(mofia_lrt(f0, fm, check_nested = FALSE, sb2010 = TRUE), "same identification")
})

test_that("F6: item plots work for items on any factor and for cross-loadings", {
  set.seed(98107); N <- 100
  d <- data.frame(x = seq(-1, 1, length.out = N)); eta <- matrix(rnorm(N * 2), N, 2)
  for (j in 1:6) d[[paste0("y", j)]] <- eta[, if (j <= 3) 1 else 2] + rnorm(N, sd = 0.6)
  L <- cbind(f1 = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE), f2 = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
  rownames(L) <- paste0("y", 1:6)
  f <- mofia(mofia_model(d, loadings = L, moderators = "x", continuous = paste0("y", 1:6)), se = FALSE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(o4 <- plot(f, "items", items = "y4"))
  expect_equal(nrow(o4), 121 * 3)
  expect_silent(plot(f, "items", items = "y6"))
})

test_that("F7: factor scores agree between the C++ and the R reference engine", {
  data(parabelief)
  f <- mofia(mofia_model(parabelief, items = c("SpiritsGhosts", "GnomesElves", "Aliens"), moderators = "Male", sparse = 0))
  s1 <- mofia_scores(f)
  old <- options(mofia.engine = "R"); on.exit(options(old))
  s2 <- mofia_scores(f)
  expect_equal(s2, s1, tolerance = 1e-8)
})
