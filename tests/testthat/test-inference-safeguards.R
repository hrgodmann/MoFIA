# Regression tests for inference safeguards, model validation and numerical accuracy.

ns <- asNamespace("MoFIA")
st0 <- function(m, link = "logit") get(".mofia_start", ns)(m, link)

test_that("the robust sandwich uses the exact information by default", {
  skip_on_cran()
  # misspecified residual-variance function: true variance 1 + 3 x^2, fitted exp(w0 + w x)
  set.seed(122); nb <- 200; eta <- rnorm(nb)
  base <- data.frame(y1 = as.integer(eta + rlogis(nb) > 0), y2 = as.integer(1.2 * eta + rlogis(nb) > 0),
                     y3 = as.integer(0.9 * eta + rlogis(nb) > 0))
  grid <- expand.grid(i = seq_len(nb), x = c(-1, 0, 1), sgn = c(-1, 1))
  dat <- base[grid$i, ]; dat$x <- grid$x; dat$y4 <- grid$sgn * sqrt(1 + 3 * grid$x^2)
  m <- mofia_model(dat, items = paste0("y", 1:4), moderators = "x", continuous = "y4", resvar_mod = "y4",
                   mean_mod = FALSE, var_mod = FALSE, sparse = 0)
  th <- st0(m); th["lambda0[y4,eta]"] <- 0
  f <- mofia(m, start = th)
  p <- "omega[y4,x]"
  expect_identical(f$information, "exact")
  expect_equal(vcov(f)[p, p], vcov(f, "MLR", "exact")[p, p])
  # the reduced bread is 5/9 instead of 4/9 per observation: SE ratio 0.8 (analytic)
  expect_equal(sqrt(vcov(f, "MLR", "reduced")[p, p] / vcov(f, "MLR", "exact")[p, p]), 0.8, tolerance = 1e-4)
})

test_that("exactly unidentified moderation patterns are refused, weak ones warned", {
  dat <- data.frame(y1 = rep(0:2, 20), y2 = rep(c(0, 2, 1), 20), y3 = rep(c(1, 1, 0, 2), 15), x = rep(0:1, 30),
                    z = round(sin(1:60), 3))
  # factor mean and all thresholds moderated by x
  expect_error(mofia_model(dat, items = c("y1", "y2", "y3"), moderators = "x", threshold_mod = TRUE, sparse = 0),
               "not identified.*mean of factor")
  # an anchor item makes it identified
  expect_silent(mofia_model(dat, items = c("y1", "y2", "y3"), moderators = "x", threshold_mod = c("y2", "y3"), sparse = 0))
  # no mean moderation: fine
  expect_silent(mofia_model(dat, items = c("y1", "y2", "y3"), moderators = "x", threshold_mod = TRUE,
                            mean_mod = FALSE, sparse = 0))
  # binary moderator: factor variance and all loadings moderated
  expect_error(mofia_model(dat, items = c("y1", "y2", "y3"), moderators = "x", loading_mod = TRUE,
                           mean_mod = FALSE, sparse = 0), "variance of factor")
  # continuous moderator: identified only through the functional form -> warning
  expect_warning(mofia_model(dat, items = c("y1", "y2", "y3"), moderators = "z", loading_mod = TRUE,
                             mean_mod = FALSE, sparse = 0), "functional form")
  # marker identification: the marker loading is unmoderated, so no scale problem
  expect_silent(mofia_model(dat, items = c("y1", "y2", "y3"), moderators = "x", loading_mod = c("y2", "y3"),
                            mean_mod = FALSE, identification = "marker", sparse = 0))
  # why: the binary-scale case is exactly likelihood-equivalent
  m <- suppressWarnings(mofia_model(dat, items = c("y1", "y2", "y3"), moderators = "z", loading_mod = TRUE,
                                    mean_mod = FALSE, sparse = 0))
  m$X[, 1] <- dat$x
  t1 <- st0(m); t1[m$par$block == "lambda0"] <- c(0.8, 1.2, 1.1); t1[m$par$block == "lambda1"] <- c(0.2, -0.1, 0.1)
  t1["delta[eta,z]"] <- 0.3; t2 <- t1; cc <- 1.5
  t2[m$par$block == "lambda1"] <- (t1[m$par$block == "lambda0"] + t1[m$par$block == "lambda1"]) / cc - t1[m$par$block == "lambda0"]
  t2["delta[eta,z]"] <- t1["delta[eta,z]"] + 2 * log(cc)
  expect_equal(mofia_loglik(m, t2), mofia_loglik(m, t1), tolerance = 1e-10)
})

test_that("degenerate moderators are refused", {
  dat <- data.frame(y1 = rep(0:1, 10), y2 = rep(c(0, 1, 1, 0), 5), x = 1, z = 1:20, w = 2 * (1:20) + 1)
  expect_error(mofia_model(dat, items = c("y1", "y2"), moderators = "x", sparse = 0), "Constant")
  expect_error(mofia_model(dat, items = c("y1", "y2"), moderators = c("z", "w"), sparse = 0), "linearly dependent")
  dat$v <- c(Inf, 1:19)
  expect_error(mofia_model(dat, items = c("y1", "y2"), moderators = "v", sparse = 0), "finite")
})

test_that("a saddle point is not reported as converged and gets no standard errors", {
  set.seed(5); N <- 300; eta <- rnorm(N)
  dat <- as.data.frame(sapply(1:4, function(j) eta + rnorm(N, sd = 0.6))); names(dat) <- paste0("y", 1:4)
  m <- mofia_model(dat, items = names(dat), continuous = names(dat))
  th <- st0(m); th[m$par$block == "lambda0"] <- 0
  th[m$par$block == "nu0"] <- colMeans(dat)
  th[m$par$block == "omega0"] <- log(colMeans(sweep(as.matrix(dat), 2, colMeans(dat))^2))
  expect_warning(f <- mofia(m, start = th), "saddle point")
  expect_false(f$converged)
  expect_identical(f$status, "saddle")
  expect_lt(f$information_eigen[["min"]], 0)
  expect_error(vcov(f), "saddle point")
  expect_true(all(is.na(summary(f)$coefficients$se)))
  good <- mofia(m)
  expect_true(good$converged); expect_identical(good$status, "ok")
  expect_gt(good$loglik - f$loglik, 100)
  expect_output(print(f), "saddle point")
})

test_that("non-positive variances give NA, never zero standard errors", {
  se_fun <- get(".se_from_vcov", ns)
  expect_identical(se_fun(diag(c(4, -1, 0, NA))), c(2, NA, NA, NA))
})

test_that("tail probabilities and scores are evaluated stably", {
  ll <- get(".mofia_ll", ns); nodes <- get(".mofia_nodes", ns)
  # probit: thresholds at 10, zero loadings -> product of known normal tail probabilities
  dat <- data.frame(y1 = rep(0:1, 5), y2 = rep(0:1, 5), y3 = rep(0:1, 5))
  m <- mofia_model(dat, items = names(dat), sparse = 0)
  th <- st0(m, "probit"); th[m$par$block == "lambda0"] <- 0; th[m$par$block == "tau"] <- 10
  ref <- sum(vapply(as.numeric(as.matrix(dat)), function(y) pnorm(10, lower.tail = y == 0, log.p = TRUE), 1))
  expect_equal(mofia_loglik(m, th, link = "probit"), ref, tolerance = 1e-10)
  ev <- ll(th, m, "probit", nodes(31, 1), scores = TRUE)
  expect_true(all(is.finite(attr(ev, "scores"))))
  expect_lt(max(abs(colSums(attr(ev, "scores")))), 1e3)
  expect_equal(as.numeric(ll(th, m, "probit", nodes(31, 1), engine = "R")), ref, tolerance = 1e-10)
  # logit and probit: a common shift of thresholds and predictor leaves the likelihood unchanged
  dat <- data.frame(y1 = rep(0:2, 4), y2 = rep(0:2, 4), y3 = rep(0:2, 4), x = rep(0:1, 6))
  m <- mofia_model(dat, items = c("y1", "y2", "y3"), moderators = "x", threshold_mod = TRUE, mean_mod = FALSE,
                   var_mod = FALSE, sparse = 0)
  m$X[] <- 1
  th <- st0(m); th[m$par$block == "tau"] <- rep(c(-0.5, 0.5), 3); th[m$par$block == "lambda0"] <- 0.1
  sh <- th; sh[m$par$block == "tau"] <- sh[m$par$block == "tau"] + 400; sh[m$par$block == "beta"] <- 400
  for (lk in c("logit", "probit")) {
    a <- ll(th, m, lk, nodes(41, 1), scores = TRUE); b <- ll(sh, m, lk, nodes(41, 1), scores = TRUE)
    expect_equal(as.numeric(b), as.numeric(a), tolerance = 1e-10)
    expect_true(all(is.finite(attr(b, "scores"))))
    bR <- ll(sh, m, lk, nodes(41, 1), scores = TRUE, engine = "R")
    expect_equal(as.numeric(bR), as.numeric(b), tolerance = 1e-10)
    expect_equal(attr(bR, "scores"), attr(b, "scores"), tolerance = 1e-8)
  }
  # predictions stay valid probabilities in the tails
  cp <- get(".cat_probs", ns)(matrix(c(-50, 0, 50), 1), c(-1, 0, 1), "probit")
  expect_equal(Reduce(`+`, cp), matrix(1, 1, 3), tolerance = 1e-12)
  expect_true(all(vapply(cp, function(p) all(p >= 0), TRUE)))
})

test_that("likelihood-ratio tests refuse incompatible comparisons", {
  skip_on_cran()
  data(parabelief); its <- c("SpiritsGhosts", "GnomesElves", "Aliens", "Sorcery")
  f0 <- mofia(mofia_model(parabelief, items = its, moderators = "Male"))
  dat <- parabelief; dat$Aliens <- rev(dat$Aliens)
  f1d <- mofia(mofia_model(dat, items = its, moderators = "Male", threshold_mod = "Aliens"))
  expect_error(mofia_lrt(f0, f1d), "item responses differ")
  expect_error(anova(f0, f1d), "item responses differ")
  f1p <- mofia(mofia_model(parabelief, items = its, moderators = "Male", threshold_mod = "Aliens"), link = "probit")
  expect_error(mofia_lrt(f0, f1p), "different links")
  f1 <- mofia(mofia_model(parabelief, items = its, moderators = "Male", threshold_mod = "Aliens"))
  expect_true(is.finite(mofia_lrt(f0, f1)$p_value))
  bad <- f1; bad$loglik <- f0$loglik - 5
  expect_error(mofia_lrt(f0, bad), "lower log-likelihood")
  expect_warning(a <- anova(f0, bad), "lower log-likelihood")
  expect_true(is.na(a$Chisq[2]))
  fz <- mofia(mofia_model(parabelief, items = its, moderators = "Age", standardize = TRUE))
  fA <- mofia(mofia_model(parabelief, items = its, moderators = c("Age", "Male"), threshold_mod = "Aliens"))
  expect_error(mofia_lrt(fz, fA), "moderator values differ")
  unconv <- f1; unconv$converged <- FALSE
  expect_error(mofia_lrt(f0, unconv), "not converged")
})

test_that("cases without any observed response are removed and do not change N or BIC", {
  skip_on_cran()
  data(parabelief); its <- c("SpiritsGhosts", "GnomesElves", "Aliens", "Sorcery")
  f <- mofia(mofia_model(parabelief, items = its, moderators = "Male"))
  d <- parabelief; d[, its] <- NA; d <- rbind(parabelief, d)
  expect_warning(m2 <- mofia_model(d, items = its, moderators = "Male"), "600 case")
  f2 <- mofia(m2)
  expect_identical(nobs(f2), nobs(f))
  expect_equal(BIC(f2), BIC(f), tolerance = 1e-8)
})

test_that("named specification matrices are aligned by name; wrong shapes are refused", {
  dat <- data.frame(y1 = rep(0:2, 10), y2 = rep(c(0, 2, 1), 10), y3 = rep(c(1, 0, 2), 10),
                    x = rep(0:1, 15), z = round(cos(1:30), 3))
  tm <- matrix(c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), 3, 2, dimnames = list(c("y1", "y2", "y3"), c("x", "z")))
  tm_rev <- tm[3:1, 2:1]
  m1 <- mofia_model(dat, items = c("y1", "y2", "y3"), moderators = c("x", "z"), threshold_mod = tm, sparse = 0)
  m2 <- mofia_model(dat, items = c("y1", "y2", "y3"), moderators = c("x", "z"), threshold_mod = tm_rev, sparse = 0)
  expect_identical(m1$par, m2$par)
  expect_true("beta[y1,z]" %in% m1$par$name)
  bad <- tm; colnames(bad) <- c("x", "w")
  expect_error(mofia_model(dat, items = c("y1", "y2", "y3"), moderators = c("x", "z"), threshold_mod = bad, sparse = 0),
               "do not match")
  expect_error(mofia_model(dat, items = c("y1", "y2", "y3"), loadings = matrix(TRUE, 1, 1), sparse = 0), "dimensions")
  expect_error(mofia_model(dat, items = c("y1", "y2", "y3"), moderators = c("x", "z"), threshold_mod = list("x"),
                           sparse = 0), "must be named")
  expect_error(mofia_model(dat, items = c("y1", "y2", "y3"), moderators = c("x", "z"), mean_mod = c(TRUE, FALSE, TRUE),
                           sparse = 0), "one value per moderator")
  lm3 <- array(FALSE, c(3, 1, 2), dimnames = list(c("y1", "y2", "y3"), "eta", c("x", "z"))); lm3["y2", "eta", "z"] <- TRUE
  m3 <- mofia_model(dat, items = c("y1", "y2", "y3"), moderators = c("x", "z"), loading_mod = lm3[, , 2:1, drop = FALSE],
                    sparse = 0)
  expect_true("lambda1[y2,eta,z]" %in% m3$par$name)
})

test_that("data-dependent formula transformations keep their training values", {
  dd <- data.frame(y1 = rep(0:2, 10), y2 = rep(c(0, 1, 2, 2, 1, 0), 5), y3 = rep(c(1, 0, 2), 10), x = seq_len(30))
  dx <- get(".design_X", ns)
  m <- mofia_model(dd, items = c("y1", "y2", "y3"), moderators = ~ scale(x), sparse = 0)
  expect_equal(dx(m, dd[1:5, ]), m$X[1:5, , drop = FALSE], ignore_attr = TRUE)
  expect_equal(dx(m, dd[1, , drop = FALSE]), m$X[1, , drop = FALSE], ignore_attr = TRUE)
  mp <- mofia_model(dd, items = c("y1", "y2", "y3"), moderators = ~ poly(x, 2), sparse = 0)
  expect_false(any(grepl(",", mp$moderators)))
  expect_equal(dx(mp, dd[3:7, ]), mp$X[3:7, ], ignore_attr = TRUE)
})

test_that("summaries print for every standard-error type", {
  skip_on_cran()
  data(parabelief)
  f <- mofia(mofia_model(parabelief, items = c("SpiritsGhosts", "GnomesElves", "Aliens"), moderators = "Male"))
  for (s in c("MLR", "ML", "MLF")) expect_output(print(summary(f, se = s)), "Standard errors")
  expect_output(print(summary(f, information = "reduced")), "reduced information")
})

test_that("simulated data keep the original moderator variables", {
  skip_on_cran()
  set.seed(3); n <- 300
  d <- data.frame(age = sample(10:39, n, TRUE), g = factor(sample(c("a", "b"), n, TRUE)))
  d$y1 <- sample(0:2, n, TRUE); d$y2 <- sample(0:2, n, TRUE); d$y3 <- sample(0:2, n, TRUE)
  m <- mofia_model(d, items = c("y1", "y2", "y3"), moderators = ~ age * g, standardize = TRUE)
  th <- st0(m)
  sim <- mofia_simulate(m, th, seed = 1)
  expect_true(all(c("age", "g") %in% names(sim)))
  expect_equal(range(sim$age), range(d$age))
  m2 <- mofia_model(sim, items = c("y1", "y2", "y3"), moderators = ~ age * g, standardize = TRUE, sparse = 0)
  expect_equal(m2$X, m$X, ignore_attr = TRUE)
  sim2 <- mofia_simulate(m, th, newdata = data.frame(age = c(12, 30), g = c("a", "b")), seed = 2)
  expect_equal(sim2$age, c(12, 30))
  sim3 <- mofia_simulate(m, th, X = m$X[1:4, ], seed = 2)
  expect_true(all(m$moderators %in% names(sim3)))
})
