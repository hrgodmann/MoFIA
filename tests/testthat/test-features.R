ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession", "Sorcery", "HauntedPlaces")

test_that("measurement syntax, moderator formulas and standardization", {
  data(parabelief)
  d <- parabelief; d$grp <- factor(rep(c("a", "b", "c"), 200)); d$AgeYears <- 50 + 15 * d$Age
  m1 <- mofia_model(d, loadings = "UB =~ SpiritsGhosts + GnomesElves + Aliens", moderators = c("Age", "Male"))
  m2 <- mofia_model(d, loadings = list(UB = c("SpiritsGhosts", "GnomesElves", "Aliens")), moderators = c("Age", "Male"))
  expect_identical(m1$par, m2$par)
  m <- mofia_model(d, items = ub[1:3], moderators = ~ AgeYears * Male + grp, standardize = TRUE)
  expect_equal(m$moderators, c("AgeYears", "Male", "grpb", "grpc", "AgeYears:Male"))
  expect_equal(mean(m$X[, "AgeYears"]), 0, tolerance = 1e-12)
  expect_equal(stats::sd(m$X[, "AgeYears"]), 1, tolerance = 1e-12)
  expect_equal(m$X[, "AgeYears:Male"], m$X[, "AgeYears"] * m$X[, "Male"])
  # new data on the original scale reproduce the model's moderator matrix
  expect_equal(MoFIA:::.design_X(m, d), m$X, ignore_attr = TRUE)
  expect_error(mofia_model(d, items = ub[1:3], moderators = "grp"), "formula")
})

test_that("sparse and unused categories are reported", {
  data(parabelief)
  d <- parabelief; d$GnomesElves[d$GnomesElves == 3] <- 2
  expect_warning(expect_warning(mofia_model(d, items = c("GnomesElves", "Aliens"), sparse = 10),
                                "never chosen"), "fewer than 10")
})

test_that("mofia_collapse merges sparse categories", {
  data(parabelief)
  d <- mofia_collapse(parabelief, c("GnomesElves", "Aliens", "Sorcery"), min_count = 20)
  for (it in c("GnomesElves", "Aliens")) {
    tab <- table(d[[it]])
    expect_true(all(tab >= 20) || length(tab) == 2)
    expect_equal(as.numeric(names(tab)), seq_along(tab))
  }
  map <- attr(d, "collapsed")
  expect_true(all(diff(map$to[map$item == "GnomesElves"]) >= 0))
})

test_that("EAP factor scores match direct numerical integration", {
  data(parabelief)
  fit <- mofia(mofia_model(parabelief, items = ub, moderators = c("Age", "Male"), threshold_mod = "Aliens"))
  sc <- mofia_scores(fit)
  expect_equal(dim(sc), c(600, 2))
  m <- fit$model; pr <- MoFIA:::.unpack(coef(fit), m); st <- MoFIA:::.structural(pr, m)
  g <- seq(-12, 12, length.out = 3001)
  for (i in c(1, 77, 450)) {
    eta <- st$mu[i, 1] + st$s[i, 1] * g; lik <- rep(1, length(g))
    for (j in seq_len(m$P)) {
      A <- MoFIA:::.item_lp(pr, j, m$X[i, , drop = FALSE], list(matrix(eta, 1)), m$loadings)
      lik <- lik * MoFIA:::.cat_probs(A, pr$tau[[j]], "logit")[[m$Y[i, j] + 1]][1, ]
    }
    w <- lik * stats::dnorm(g); w <- w / sum(w); mu <- sum(w * eta)
    expect_equal(sc$eta[i], mu, tolerance = 1e-4)
    expect_equal(sc$eta_sd[i], sqrt(sum(w * (eta - mu)^2)), tolerance = 1e-4)
  }

  # predictions, intervals, tidiers
  p <- predict(fit, newdata = data.frame(Age = c(-1, 1), Male = c(0, 1)), type = "prob")
  expect_equal(unname(rowSums(p$Aliens)), c(1, 1), tolerance = 1e-10)
  e <- predict(fit, newdata = data.frame(Age = c(-1, 1), Male = c(0, 1)))
  expect_equal(unname(e[, "Aliens"]), unname(drop(p$Aliens %*% 1:5)))
  expect_equal(dim(fitted(fit)), c(600, 6))
  ci <- confint(fit)
  expect_true(all(ci[, 1] < coef(fit) & coef(fit) < ci[, 2]))
  td <- tidy.mofia(fit, conf.int = TRUE)
  expect_equal(td$estimate, unname(coef(fit)))
  expect_equal(glance.mofia(fit)$logLik, fit$loglik)
})

test_that("score tests agree with likelihood-ratio tests", {
  data(parabelief)
  fit0 <- mofia(mofia_model(parabelief, items = ub, moderators = c("Age", "Male")))
  st <- mofia_score_test(fit0)
  expect_true(all(c("beta[Aliens,Age]", "lambda1[Sorcery,eta,Age]") %in% st$test))
  fit1 <- mofia(mofia_model(parabelief, items = ub, moderators = c("Age", "Male"), threshold_mod = list(Aliens = "Age")))
  lr <- 2 * (fit1$loglik - fit0$loglik)
  r <- st[st$test == "beta[Aliens,Age]", ]
  expect_equal(r$statistic_ML, lr, tolerance = 0.1)
  expect_equal(r$epc, unname(coef(fit1)["beta[Aliens,Age]"]), tolerance = 0.1)
  gi <- mofia_score_test(fit0, by = "item")
  expect_true(all(gi$df == 2))
  expect_error(mofia_score_test(fit0, "gamma[eta,Age]"), "Already in the model")
})

test_that("random starts and plots run", {
  data(parabelief)
  set.seed(2)
  fit <- mofia(mofia_model(parabelief, items = ub[1:4], moderators = ~ Age + Male, loading_mod = "Aliens"), starts = 2)
  expect_equal(nrow(fit$starts), 2)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_s3_class(plot(fit, "factors", moderator = "Age"), "data.frame")
  expect_s3_class(plot(fit, "loadings", moderator = "Male"), "data.frame")
  expect_s3_class(plot(fit, "items", moderator = "Age", items = "Aliens"), "data.frame")
})
