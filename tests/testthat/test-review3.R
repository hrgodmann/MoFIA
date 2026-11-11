# Regression tests for the engine review of 11 November 2026 (audit of 0.1.0 at
# commit 08b9acb). R1: unit-invariant information checks; R2: exact mean/DIF aliases
# with moderated loadings; R3: row identity of scores and predictions; R4: exclusions
# before moderator classification and scaling; R5: thread capability and citation.

ns <- asNamespace("MoFIA")
info_inverse <- function(...) get(".info_inverse", ns)(...)
with_check_off <- function(expr) {
  st <- get(".mofia_state", ns); old <- st$check_identification
  st$check_identification <- FALSE; on.exit(st$check_identification <- old)
  expr
}

# ---- R1 ------------------------------------------------------------------------------
test_that("R1: information classification and inverse do not depend on parameter units", {
  S <- diag(c(1, 1e5, 1e-4, 1e3))                       # change of parameter units
  A <- matrix(c(4, 1, 0.5, 0.2, 1, 3, 0.3, 0.1, 0.5, 0.3, 2, 0.4, 0.2, 0.1, 0.4, 1), 4)
  As <- S %*% A %*% S
  # the raw eigenvalue ratio is far below the tolerance, yet As is positive definite
  ev <- eigen(As, only.values = TRUE)$values
  expect_lt(min(ev) / max(ev), 1e-8)
  r <- info_inverse(As)
  expect_identical(r$status, "ok")
  expect_equal(diag(S) * r$inverse * rep(diag(S), each = 4), solve(A), tolerance = 1e-10)
  expect_equal(r$eigen, info_inverse(A)$eigen, tolerance = 1e-12)   # same scaled spectrum
  # genuine failures stay failures, whatever the units
  set.seed(1); L <- matrix(stats::rnorm(12), 4, 3)
  sing <- S %*% tcrossprod(L) %*% S                         # rank 3
  expect_identical(info_inverse(sing)$status, "singular")
  expect_null(info_inverse(sing)$inverse)
  near <- S %*% tcrossprod(cbind(L, 1e-7 * stats::rnorm(4))) %*% S   # numerically rank 3
  expect_identical(info_inverse(near)$status, "singular")
  ind <- A; ind[2, 2] <- -1
  expect_identical(info_inverse(S %*% ind %*% S)$status, "saddle")
  zero <- A; zero[3, ] <- zero[, 3] <- 0
  expect_identical(info_inverse(zero)$status, "singular")
  zero_off <- zero; zero_off[3, 1] <- zero_off[1, 3] <- 0.5      # zero curvature but coupled: indefinite
  expect_identical(info_inverse(zero_off)$status, "saddle")
  nonfin <- A; nonfin[1, 2] <- nonfin[2, 1] <- NaN
  expect_identical(info_inverse(nonfin)$status, "singular")
  expect_identical(info_inverse(matrix(Inf, 2, 2))$status, "singular")
})

test_that("R1: a moderator in other units gives the same fit, transformed covariances and tests", {
  skip_on_cran()
  data(parabelief)
  its <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession", "Sorcery")
  d <- parabelief; d$AgeK <- d$Age * 1e4                     # same moderator, other units
  f1 <- mofia(mofia_model(d, items = its, moderators = c("Age", "Male"), threshold_mod = "Aliens",
                          loading_mod = list(Sorcery = "Age")))
  f2 <- mofia(mofia_model(d, items = its, moderators = c("AgeK", "Male"), threshold_mod = "Aliens",
                          loading_mod = list(Sorcery = "AgeK")))
  expect_identical(f2$status, "ok"); expect_true(f2$converged)
  # the raw spectrum alone would have called it singular (kept as a diagnostic)
  expect_lt(f2$information_eigen[["min"]] / f2$information_eigen[["max"]], 1e-8)
  expect_equal(f2$loglik, f1$loglik, tolerance = 1e-8)
  k <- ifelse(grepl("Age", names(coef(f1))), 1e4, 1)          # theta_AgeK = theta_Age / 1e4
  # optimizer precision: both fits stop with Newton decrement < 1e-6 (log-likelihood
  # units), i.e. a few 1e-4 standard errors from the maximum
  expect_equal(unname(coef(f2) * k), unname(coef(f1)), tolerance = 1e-4)
  for (se in c("ML", "MLR", "MLF")) {
    V1 <- vcov(f1, se); V2 <- vcov(f2, se) * outer(k, k)
    expect_equal(unname(sqrt(diag(V2))), unname(sqrt(diag(V1))), tolerance = 1e-4)
  }
  expect_equal(f2$scaling, f1$scaling, tolerance = 1e-4)
  expect_lt(f2$newton_decrement, 1e-6)
  w1 <- mofia_wald(f1, c("beta[Aliens,Age]", "beta[Aliens,Male]"))
  w2 <- mofia_wald(f2, c("beta[Aliens,AgeK]", "beta[Aliens,Male]"))
  expect_equal(w2$statistic, w1$statistic, tolerance = 1e-4)
  s1 <- mofia_score_test(f1, c("beta[GnomesElves,Age]", "lambda1[Aliens,eta,Age]"))
  s2 <- mofia_score_test(f2, c("beta[GnomesElves,AgeK]", "lambda1[Aliens,eta,AgeK]"))
  expect_equal(s2$statistic, s1$statistic, tolerance = 1e-4)
  expect_equal(s2$statistic_ML, s1$statistic_ML, tolerance = 1e-4)
  expect_equal(sort(s2$epc * 1e4), sort(s1$epc), tolerance = 1e-4)
  lr1 <- mofia_lrt(mofia(mofia_drop(f1, c("beta[Aliens,Age]", "beta[Aliens,Male]"))), f1)
  lr2 <- mofia_lrt(mofia(mofia_drop(f2, c("beta[Aliens,AgeK]", "beta[Aliens,Male]"))), f2)
  expect_equal(lr2$statistic, lr1$statistic, tolerance = 1e-4)
})

test_that("R1: age in years and its square match the standardized model after transformation", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  data(parabelief)
  ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession", "Sorcery", "HauntedPlaces")
  d <- parabelief; d$AgeY <- round(48 + 17 * d$Age)
  fr <- mofia(mofia_model(d, items = ub, moderators = ~ AgeY + I(AgeY^2) + Male, threshold_mod = "DevilPossession"))
  fs <- mofia(mofia_model(d, items = ub, moderators = ~ AgeY + I(AgeY^2) + Male, standardize = TRUE,
                          threshold_mod = "DevilPossession"))
  expect_identical(fr$status, "ok"); expect_true(fr$converged)      # was "singular" in 0.1.0
  expect_equal(fr$loglik, fs$loglik, tolerance = 1e-6)
  # standardized parameters as a function of the raw ones (origin and scale of the
  # factor move with the reference age; see the derivation in the review notes)
  m <- mean(d$AgeY); s <- stats::sd(d$AgeY)
  g <- function(t) {
    o <- t; G <- function(n) t[[n]]
    cmu <- G("gamma[eta,AgeY]") * m + G("gamma[eta,I(AgeY^2)]") * m^2
    k <- exp((G("delta[eta,AgeY]") * m + G("delta[eta,I(AgeY^2)]") * m^2) / 2)
    o["gamma[eta,AgeY]"] <- (G("gamma[eta,AgeY]") + 2 * G("gamma[eta,I(AgeY^2)]") * m) * s / k
    o["gamma[eta,I(AgeY^2)]"] <- G("gamma[eta,I(AgeY^2)]") * s^2 / k
    o["gamma[eta,Male]"] <- G("gamma[eta,Male]") / k
    o["delta[eta,AgeY]"] <- (G("delta[eta,AgeY]") + 2 * G("delta[eta,I(AgeY^2)]") * m) * s
    o["delta[eta,I(AgeY^2)]"] <- G("delta[eta,I(AgeY^2)]") * s^2
    for (it in ub) {
      la <- G(sprintf("lambda0[%s,eta]", it)); o[sprintf("lambda0[%s,eta]", it)] <- la * k; cb <- 0
      if (it == "DevilPossession") {
        b1 <- G("beta[DevilPossession,AgeY]"); b2 <- G("beta[DevilPossession,I(AgeY^2)]"); cb <- b1 * m + b2 * m^2
        o["beta[DevilPossession,AgeY]"] <- (b1 + 2 * b2 * m) * s; o["beta[DevilPossession,I(AgeY^2)]"] <- b2 * s^2
      }
      tn <- grep(sprintf("^tau\\[%s,", it), names(t), value = TRUE); o[tn] <- t[tn] - cb - la * cmu
    }
    o
  }
  th <- g(coef(fr))
  expect_equal(th, coef(fs)[names(th)], tolerance = 1e-3)
  J <- numDeriv::jacobian(g, coef(fr))
  for (se in c("ML", "MLR")) {
    Vt <- J %*% vcov(fr, se) %*% t(J)
    expect_equal(sqrt(diag(Vt)), unname(sqrt(diag(vcov(fs, se)[names(th), names(th)]))), tolerance = 1e-3)
  }
  expect_equal(fr$scaling, fs$scaling, tolerance = 1e-4)
  expect_equal(mofia_wald(fr, "^beta\\[DevilPossession,", regex = TRUE)$statistic,
               mofia_wald(fs, "^beta\\[DevilPossession,", regex = TRUE)$statistic, tolerance = 1e-3)
})

test_that("R1: a truly unidentified model is still reported as singular", {
  skip_on_cran()
  data(parabelief)
  ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession", "Sorcery", "HauntedPlaces")
  m <- with_check_off(mofia_model(parabelief, items = ub, moderators = "Male", var_mod = FALSE,
                                  threshold_mod = TRUE, loading_mod = TRUE))
  expect_warning(f <- mofia(m), "singular")
  expect_identical(f$status, "singular")
  expect_false(f$converged)
  expect_error(vcov(f), "singular")
})

test_that("R1: score tests flag candidates without information given the model, in any units", {
  skip_on_cran()
  data(parabelief)
  its <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession")
  d <- parabelief; d$MaleK <- d$Male * 1e4
  for (mod in c("Male", "MaleK")) {
    f <- mofia(mofia_model(d, items = its, moderators = mod, threshold_mod = its[-1]))
    # adding the last threshold shift makes the mean shift unidentified
    st <- mofia_score_test(f, c(sprintf("beta[%s,%s]", its[1], mod), sprintf("lambda1[%s,eta,%s]", its[1], mod)))
    expect_true(is.na(st$statistic[st$test == sprintf("beta[%s,%s]", its[1], mod)]))
    expect_false(is.na(st$statistic[st$test == sprintf("lambda1[%s,eta,%s]", its[1], mod)]))
  }
})

# ---- R2 ------------------------------------------------------------------------------
id_data <- function(N = 240, seed = 3) {
  set.seed(seed)
  data.frame(y1 = sample(0:2, N, TRUE), y2 = sample(0:2, N, TRUE), y3 = sample(0:2, N, TRUE),
             y4 = sample(0:2, N, TRUE), y5 = sample(0:2, N, TRUE),
             z = rep(0:1, N / 2), a = rep(c(0, 0, 1, 1), N / 4), w = round(stats::rnorm(N), 3))
}
id_cond <- function(expr) tryCatch(withCallingHandlers(expr, warning = function(w) {
  if (inherits(w, "mofia_identification_warning")) stop(w) }), mofia_identification_error = function(e) e,
  mofia_identification_warning = function(w) w)
its3 <- c("y1", "y2", "y3")

test_that("R2: mean + all thresholds + loadings by the same binary moderator is refused (exact alias)", {
  dat <- id_data()
  e <- id_cond(mofia_model(dat, items = its3, moderators = "z", threshold_mod = TRUE, loading_mod = TRUE,
                           var_mod = FALSE, sparse = 0))
  expect_s3_class(e, "mofia_identification_error")
  expect_identical(e$issues$status, "not_identified")
  expect_identical(e$issues$parameter, "mean")
  expect_identical(e$issues$moderators, "z")
  expect_match(conditionMessage(e), "offset exactly")
  # some loadings moderated is no different: z * z = z
  e <- id_cond(mofia_model(dat, items = its3, moderators = "z", threshold_mod = TRUE, loading_mod = "y2",
                           var_mod = FALSE, sparse = 0))
  expect_s3_class(e, "mofia_identification_error")
  # two binary moderators with their interaction (closed under products)
  e <- id_cond(mofia_model(dat, items = its3, moderators = ~ z * a, threshold_mod = TRUE, loading_mod = TRUE,
                           var_mod = FALSE, sparse = 0))
  expect_s3_class(e, "mofia_identification_error")
  expect_match(e$issues$moderators, "z, a, z:a")
  # x and x^2 in the threshold moderation absorb the product x * x
  e <- id_cond(mofia_model(dat, items = its3, moderators = ~ w + I(w^2), mean_mod = "w", threshold_mod = TRUE,
                           loading_mod = list(y1 = "w"), var_mod = FALSE, sparse = 0))
  expect_s3_class(e, "mofia_identification_error")
  expect_identical(e$issues$moderators, "w")
})

test_that("R2: the refused patterns are exact likelihood invariances", {
  dat <- id_data()
  for (spec in list(list(mods = "z", lm = TRUE), list(mods = ~ z * a, lm = TRUE))) {
    m <- with_check_off(mofia_model(dat, items = its3, moderators = spec$mods, threshold_mod = TRUE,
                                    loading_mod = spec$lm, var_mod = FALSE, sparse = 0))
    set.seed(9); th <- get(".mofia_start", ns)(m)
    b <- m$par$block; th[b == "lambda1"] <- stats::rnorm(sum(b == "lambda1"), 0, 0.3)
    th[b %in% c("gamma", "beta")] <- stats::rnorm(sum(b %in% c("gamma", "beta")), 0, 0.3)
    X <- m$X; cc <- stats::rnorm(ncol(X)); th2 <- th
    th2[b == "gamma"] <- th[b == "gamma"] + cc
    pr <- get(".unpack", ns)(th, m)
    for (j in seq_len(m$P)) {
      shift <- (pr$L0[j, 1] + drop(X %*% pr$L1[j, 1, ])) * drop(X %*% cc)   # change of the linear predictor
      co <- stats::lm.fit(cbind(1, X), shift)$coefficients
      expect_lt(max(abs(cbind(1, X) %*% co - shift)), 1e-10)                # absorbed exactly
      th2[m$par$block == "tau" & m$par$item %in% j] <- th[m$par$block == "tau" & m$par$item %in% j] + co[1]
      th2[m$par$block == "beta" & m$par$item %in% j] <- th[m$par$block == "beta" & m$par$item %in% j] - co[-1]
    }
    expect_equal(mofia_loglik(m, th2), mofia_loglik(m, th), tolerance = 1e-10)
  }
})

test_that("R2: weakly identified and identified patterns are not refused", {
  dat <- id_data()
  # continuous moderator: the product w * w is not absorbed -> identified through the loadings only
  w <- id_cond(mofia_model(dat, items = its3, moderators = "w", threshold_mod = TRUE, loading_mod = TRUE,
                           var_mod = FALSE, sparse = 0))
  expect_s3_class(w, "mofia_identification_warning")
  expect_identical(w$issues$status, "identified_by_loading_moderation")
  expect_match(conditionMessage(w), "only through the moderated loadings")
  m <- suppressWarnings(mofia_model(dat, items = its3, moderators = "w", threshold_mod = TRUE, loading_mod = TRUE,
                                    var_mod = FALSE, sparse = 0))
  expect_identical(m$identification_issues$status, "identified_by_loading_moderation")
  # binary mean moderator, loadings moderated by a continuous one: z * w is not absorbed
  w <- id_cond(mofia_model(dat, items = its3, moderators = c("z", "w"), mean_mod = "z",
                           threshold_mod = list(y1 = "z", y2 = "z", y3 = "z"), loading_mod = list(y1 = "w"),
                           var_mod = FALSE, sparse = 0))
  expect_s3_class(w, "mofia_identification_warning")
  expect_identical(w$issues$moderators, "z")
  # anchor item, no mean moderation, or mean on another moderator: no issue
  expect_silent(m <- mofia_model(dat, items = its3, moderators = "z", threshold_mod = c("y2", "y3"),
                                 loading_mod = TRUE, var_mod = FALSE, sparse = 0))
  expect_identical(nrow(m$identification_issues), 0L)
  expect_silent(mofia_model(dat, items = its3, moderators = "z", threshold_mod = TRUE, loading_mod = TRUE,
                            mean_mod = FALSE, var_mod = FALSE, sparse = 0))
  expect_silent(mofia_model(dat, items = its3, moderators = c("z", "a"), mean_mod = "a",
                            threshold_mod = list(y1 = "z", y2 = "z", y3 = "z"), loading_mod = TRUE,
                            var_mod = FALSE, sparse = 0))
  # the existing functional-form warning (scale) keeps its status
  w <- id_cond(mofia_model(dat, items = its3, moderators = "w", loading_mod = TRUE, mean_mod = FALSE, sparse = 0))
  expect_identical(w$issues$status, "identified_by_functional_form")
  # and the binary scale alias stays an error
  e <- id_cond(mofia_model(dat, items = its3, moderators = "z", loading_mod = TRUE, mean_mod = FALSE, sparse = 0))
  expect_s3_class(e, "mofia_identification_error")
  expect_identical(e$issues$parameter, "variance")
})

test_that("R2: cross-loadings are checked per factor", {
  dat <- id_data()
  ld <- list(f1 = c("y1", "y2", "y3"), f2 = c("y3", "y4", "y5"))
  e <- id_cond(mofia_model(dat, loadings = ld, moderators = "z", threshold_mod = TRUE, loading_mod = TRUE,
                           var_mod = FALSE, sparse = 0))
  expect_s3_class(e, "mofia_identification_error")
  expect_setequal(e$issues$factor, c("f1", "f2"))
  # an anchor for f1 only: f2 is still aliased
  e <- id_cond(mofia_model(dat, loadings = ld, moderators = "z", threshold_mod = c("y2", "y3", "y4", "y5"),
                           loading_mod = TRUE, var_mod = FALSE, sparse = 0))
  expect_identical(e$issues$factor, "f2")
})

# ---- R3 / R4 -------------------------------------------------------------------------
row_data <- function() {
  set.seed(11); N <- 60; eta <- stats::rnorm(N)
  d <- data.frame(y1 = as.integer(eta + stats::rlogis(N) > 0), y2 = as.integer(1.2 * eta + stats::rlogis(N) > 0),
                  y3 = as.integer(0.8 * eta + stats::rlogis(N) > 0), y4 = as.integer(eta + stats::rlogis(N) > 0.3),
                  age = round(stats::runif(N, 20, 70)), g = rep(0:1, N / 2))
  rownames(d) <- sprintf("id%03d", seq_len(N) * 7)
  d$age[c(2, 7)] <- NA                                  # missing moderator
  d[c(4, 9, 31), c("y1", "y2", "y3", "y4")] <- NA       # no response
  d$age[31] <- 1000                                     # ... including an outlier in age
  d$y2[5] <- NA; d$y3[12] <- NA                         # partial responses stay
  d
}

test_that("R3: original row positions and names are kept through both exclusions and rebuilds", {
  d <- row_data(); its <- c("y1", "y2", "y3", "y4")
  expect_warning(expect_warning(m <- mofia_model(d, items = its, moderators = c("age", "g"), standardize = TRUE),
                                "2 case\\(s\\) with missing moderator"), "3 case\\(s\\) without any observed")
  keep <- setdiff(seq_len(60), c(2, 7, 4, 9, 31))
  expect_identical(m$cases$index, keep)
  expect_identical(m$cases$names, rownames(d)[keep])
  expect_identical(m$cases$n_input, 60L)
  expect_identical(m$cases$excluded_moderator, c(2L, 7L))
  expect_identical(m$cases$excluded_response, c(4L, 9L, 31L))
  expect_identical(rownames(m$Y), rownames(d)[keep])
  expect_identical(m$N, length(keep))
  m0 <- mofia_drop(m, "gamma[eta,g]")
  expect_identical(m0$cases, m$cases)
})

test_that("R3: scores and predictions carry the identities of their rows", {
  skip_on_cran()
  d <- row_data(); its <- c("y1", "y2", "y3", "y4")
  m <- suppressWarnings(mofia_model(d, items = its, moderators = c("age", "g"), standardize = TRUE, var_mod = FALSE))
  f <- mofia(m)
  sc <- mofia_scores(f)
  expect_s3_class(sc, "data.frame"); expect_identical(nrow(sc), m$N)          # unchanged type and count
  expect_identical(rownames(sc), m$cases$names)
  expect_identical(attr(sc, "row_index"), m$cases$index)
  # the same cases fitted from manually prefiltered data give the same scores, row by row
  dp <- d[m$cases$index, ]
  fp <- mofia(mofia_model(dp, items = its, moderators = c("age", "g"), standardize = TRUE, var_mod = FALSE))
  expect_equal(coef(f), coef(fp), tolerance = 1e-6)
  scp <- mofia_scores(fp)
  expect_equal(sc$eta[match(rownames(dp), rownames(sc))], scp$eta, tolerance = 1e-6)
  # filling a full-length column by index
  full <- rep(NA_real_, nrow(d)); full[attr(sc, "row_index")] <- sc$eta
  expect_equal(full[m$cases$index], sc$eta)
  expect_true(all(is.na(full[c(2, 7, 4, 9, 31)])))
  e <- predict(f); p <- predict(f, type = "prob")
  expect_identical(dim(e), c(m$N, 4L)); expect_identical(rownames(e), m$cases$names)
  expect_identical(attr(e, "row_index"), m$cases$index)
  expect_identical(rownames(p$y1), m$cases$names); expect_identical(attr(p$y1, "row_index"), m$cases$index)
  expect_equal(unclass(fitted(f))[, ], unclass(e)[, ])
  nd <- data.frame(age = c(30, 50, 40), g = c(0, 1, 1), row.names = c("a", "b", "c"))
  en <- predict(f, newdata = nd)
  expect_identical(rownames(en), c("a", "b", "c")); expect_identical(attr(en, "row_index"), 1:3)
  expect_equal(unname(en["b", ]), unname(predict(f, newdata = nd[2, ])[1, ]))
  expect_error(predict(f, newdata = data.frame(age = c(30, NA), g = 0:1)), "Missing moderator")
})

test_that("R3: posterior covariance rows are named like the score rows (two factors)", {
  skip_on_cran()
  fx <- sim_two_factor(N = 300, seed = 4)
  d <- fx$data; rownames(d) <- paste0("p", seq_len(nrow(d))); d[c(5, 50), paste0("y", 1:6)] <- NA
  m <- suppressWarnings(do.call(mofia_model, c(list(data = d, sparse = 0), fx$args)))
  f <- mofia(m, start = fx$theta, se = FALSE)
  sc <- mofia_scores(f)
  expect_identical(dimnames(attr(sc, "cov"))[[1]], rownames(sc))
  expect_identical(rownames(sc)[1:5], c("p1", "p2", "p3", "p4", "p6"))
})

test_that("R4: moderator classification and standardization use only the fitted cases", {
  d <- row_data(); its <- c("y1", "y2", "y3", "y4")
  m <- suppressWarnings(mofia_model(d, items = its, moderators = c("age", "g"), standardize = TRUE))
  kept <- d[m$cases$index, ]
  # the outlier age = 1000 in a case without responses no longer enters the scaling
  expect_equal(m$design$scaling$age[["mean"]], mean(kept$age))
  expect_equal(m$design$scaling$age[["sd"]], stats::sd(kept$age))
  expect_equal(mean(m$X[, "age"]), 0, tolerance = 1e-12); expect_equal(stats::sd(m$X[, "age"]), 1, tolerance = 1e-12)
  mp <- mofia_model(kept, items = its, moderators = c("age", "g"), standardize = TRUE)
  expect_identical(m$X, mp$X); expect_identical(m$design, mp$design)
  # a moderator whose third value occurs only in excluded cases is binary for the fitted cases
  d$h <- d$g; d$h[c(4, 7)] <- 5
  m <- suppressWarnings(mofia_model(d, items = its, moderators = c("h", "age"), standardize = TRUE))
  expect_identical(m$design$info$h$type, "binary")
  expect_null(m$design$scaling$h)
  # data-dependent formula terms (scale(), poly()) are learned from the fitted cases
  m <- suppressWarnings(mofia_model(d, items = its, moderators = ~ poly(age, 2) + scale(g)))
  mp <- mofia_model(d[m$cases$index, ], items = its, moderators = ~ poly(age, 2) + scale(g))
  expect_equal(m$X, mp$X)
  expect_equal(colMeans(m$X), c(0, 0, 0), tolerance = 1e-12, ignore_attr = TRUE)
})

# ---- R5 ------------------------------------------------------------------------------
test_that("R5: requested threads are reported, and give the same result where OpenMP is available", {
  skip_on_cran()
  data(parabelief)
  m <- mofia_model(parabelief, items = c("SpiritsGhosts", "GnomesElves", "Aliens"), moderators = "Male")
  omp <- get("mofia_openmp_threads", ns)()
  f1 <- mofia(m, threads = 1)
  expect_identical(f1$threads[["used"]], 1L)
  if (omp < 1) {
    # serial build: the request is reported, not silently ignored (fallback coverage only)
    expect_message(f2 <- mofia(m, threads = 2), "no OpenMP support")
    expect_identical(f2$threads[["used"]], 1L)
    expect_equal(f2$loglik, f1$loglik)
    skip("This build has no OpenMP support; multithreaded results were not compared.")
  }
  expect_silent(f2 <- mofia(m, threads = 2))
  expect_identical(f2$threads[["used"]], 2L)
  expect_equal(f2$loglik, f1$loglik, tolerance = 1e-10)
  expect_equal(coef(f2), coef(f1), tolerance = 1e-8)
})

test_that("R5: the citation year comes from release metadata, not the current date", {
  file <- system.file("CITATION", package = "MoFIA")
  skip_if(file == "")
  a <- utils::readCitationFile(file, meta = list(Version = "0.1.0", Date = "2025-03-04"))
  expect_identical(a$year, "2025")
  b <- utils::readCitationFile(file, meta = list(Version = "0.1.0"))
  expect_null(b$year)
  expect_match(format(b, style = "textVersion"), "n\\.d\\.")
})
