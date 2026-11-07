test_that("parameters are recovered and methods work", {
  s <- sim_two_factor(N = 1000, seed = 11)
  m <- do.call(mofia_model, c(list(data = s$data), s$args))
  fit <- mofia(m)
  expect_s3_class(fit, "mofia")
  expect_true(fit$converged)
  z <- (coef(fit) - s$theta[names(coef(fit))]) / sqrt(diag(vcov(fit)))
  expect_lt(max(abs(z)), 4)
  expect_equal(dim(vcov(fit, "ML", "exact")), c(fit$npar, fit$npar))
  expect_equal(attr(logLik(fit), "df"), fit$npar)
  expect_equal(stats::AIC(fit), -2 * fit$loglik + 2 * fit$npar)
  expect_equal(nobs(fit), 1000)
  expect_output(print(summary(fit)), "Loading moderation")

  # nested comparison: drop the DIF effect
  fit0 <- mofia(mofia_drop(m, "beta[y3,x2]"), start = coef(fit))
  a <- anova(fit0, fit)
  t1 <- mofia_lrt(fit0, fit)
  expect_equal(a$Chisq[2], t1$statistic)
  expect_equal(a$Df[2], 1)
  expect_equal(mofia_lrt(fit0, fit, scaled = FALSE)$LR, 2 * (fit$loglik - fit0$loglik))
  w <- mofia_wald(fit, "beta[y3,x2]")
  expect_equal(w$df, 1)
  expect_lt(w$p_value, 0.05)
  expect_error(mofia_lrt(fit, fit0), "fewer parameters")
})

test_that("simulated data keep the original coding", {
  data(parabelief)
  m <- mofia_model(parabelief, items = c("Aliens", "Sorcery", "HauntedPlaces"), moderators = "Male")
  th <- stats::setNames(numeric(nrow(m$par)), m$par$name)
  th[grepl("^tau", names(th))] <- rep(c(-1, 0, 1, 2), 3); th[grepl("^lambda0", names(th))] <- 1
  sim <- mofia_simulate(m, th, seed = 3)
  expect_equal(dim(sim), c(600, 4))
  expect_true(all(unlist(sim[1:3]) %in% 1:5))
  expect_equal(dim(attr(sim, "eta")), c(600, 1))
})

test_that("the paranormal-belief metric model reproduces its reference log-likelihood", {
  skip_on_cran()
  data(parabelief)
  eha <- c("HealingHands", "LifeCoursePossessions", "Graphology", "Dowsing", "MindReading", "Willpower",
           "SeeingFuture", "SeeingEvents")
  ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession", "Sorcery", "HauntedPlaces")
  d <- parabelief; d$Age <- as.numeric(scale(d$Age))
  m <- mofia_model(d, loadings = list(EHA = eha, UB = ub), moderators = c("Age", "Male"),
                   threshold_mod = TRUE, mean_mod = FALSE, var_mod = TRUE, cor_mod = TRUE)
  fit <- mofia(m, threads = 2)
  expect_lt(abs(fit$loglik + 8394.7028), 1e-3)
  expect_equal(fit$npar, 105)
})
