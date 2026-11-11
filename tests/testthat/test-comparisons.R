# Independent comparisons of the likelihood, estimates and standard errors (from the
# audit of 11 November 2026; full scripts and logs in validation/audit-2026-11-11).
# Reference software is optional: these tests are skipped when it is not installed.

test_that("continuous items, binary moderator: same fit and ML/MLR SEs as lavaan multigroup", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(2026); N <- 600; g <- rep(0:1, each = N / 2)
  eta <- stats::rnorm(N, 0.4 * g, exp(0.3 * g / 2))
  lam <- c(1, 0.8, 1.2, 0.9, 1.1); nu <- c(0, 0.3, -0.2, 0.1, 0.5)
  Y <- sapply(1:5, function(j) nu[j] + lam[j] * eta + 0.3 * (j == 4) * g + stats::rnorm(N, sd = 0.7))
  d <- data.frame(Y, g = g); its <- paste0("y", 1:5); names(d)[1:5] <- its
  fm <- mofia(mofia_model(d, items = its, moderators = "g", continuous = its, threshold_mod = "y4"))
  mod <- "f =~ NA*y1 + y2 + y3 + y4 + y5
          f ~~ c(1, NA)*f
          f ~ c(0, NA)*1"
  fl <- lavaan::cfa(mod, d, group = "g", group.label = c("0", "1"), group.equal = c("loadings", "intercepts", "residuals"),
                    group.partial = "y4~1", estimator = "ML", information = "observed",
                    observed.information = "hessian", meanstructure = TRUE)
  expect_equal(fm$loglik, unname(lavaan::fitMeasures(fl, "logl")), tolerance = 1e-8)
  pe <- lavaan::parameterEstimates(fl)
  sel <- function(p, lhs, op, rhs = NULL, grp) p[p$lhs == lhs & p$op == op & (is.null(rhs) | p$rhs %in% rhs) & p$group == grp, ]
  expect_equal(coef(fm)[["lambda0[y2,eta]"]], sel(pe, "f", "=~", "y2", 1)$est, tolerance = 1e-4)
  expect_equal(coef(fm)[["gamma[eta,g]"]], sel(pe, "f", "~1", NULL, 2)$est, tolerance = 1e-4)
  expect_equal(coef(fm)[["delta[eta,g]"]], log(sel(pe, "f", "~~", "f", 2)$est), tolerance = 1e-4)
  expect_equal(coef(fm)[["beta[y4,g]"]], sel(pe, "y4", "~1", NULL, 2)$est - sel(pe, "y4", "~1", NULL, 1)$est,
               tolerance = 1e-4)
  se <- sqrt(diag(vcov(fm, "ML")))
  expect_equal(se[["lambda0[y2,eta]"]], sel(pe, "f", "=~", "y2", 1)$se, tolerance = 1e-4)
  expect_equal(se[["gamma[eta,g]"]], sel(pe, "f", "~1", NULL, 2)$se, tolerance = 1e-4)
  per <- lavaan::parameterEstimates(stats::update(fl, estimator = "MLR"))
  se <- sqrt(diag(vcov(fm, "MLR")))
  expect_equal(se[["lambda0[y2,eta]"]], sel(per, "f", "=~", "y2", 1)$se, tolerance = 1e-4)
  expect_equal(se[["gamma[eta,g]"]], sel(per, "f", "~1", NULL, 2)$se, tolerance = 1e-4)
})

test_that("ordinal items, binary moderator: same fit as mirt multipleGroup (graded)", {
  skip_on_cran()
  skip_if_not_installed("mirt")
  set.seed(77); N <- 1000; g <- rep(0:1, each = N / 2)
  eta <- stats::rnorm(N, -0.3 * g, exp(-0.4 * g / 2))
  a <- c(1.5, 1.2, 1.0, 1.8, 0.9, 1.3)
  tau <- rbind(c(-1, 0.2, 1.3), c(-0.5, 0.5, 1.5), c(-1.5, -0.3, 0.8), c(-0.8, 0.4, 1.1), c(-1.2, 0, 1.4), c(-0.6, 0.6, 1.6))
  Y <- sapply(1:6, function(j) findInterval(a[j] * eta + stats::rlogis(N), tau[j, ], left.open = TRUE))
  d <- data.frame(Y, g = g); its <- paste0("y", 1:6); names(d)[1:6] <- its
  fm <- mofia(mofia_model(d, items = its, moderators = "g", sparse = 0))
  fmir <- mirt::multipleGroup(d[, its], 1, group = factor(d$g), itemtype = "graded",
                              invariance = c("slopes", "intercepts", "free_means", "free_var"),
                              SE = TRUE, SE.type = "Oakes", verbose = FALSE, technical = list(NCYCLES = 5000),
                              TOL = 1e-8, quadpts = 121)
  expect_equal(fm$loglik, as.numeric(mirt::extract.mirt(fmir, "logLik")), tolerance = 1e-8)
  cm <- mirt::coef(fmir, simplify = TRUE)
  expect_equal(coef(fm)[["lambda0[y1,eta]"]], cm$`0`$items["y1", "a1"], tolerance = 1e-4)
  expect_equal(unname(coef(fm)[sprintf("tau[y1,%d]", 1:3)]), unname(-cm$`0`$items["y1", c("d1", "d2", "d3")]),
               tolerance = 1e-4)
  expect_equal(coef(fm)[["gamma[eta,g]"]], unname(cm$`1`$means[1]), tolerance = 1e-4)
  expect_equal(exp(coef(fm)[["delta[eta,g]"]]), unname(cm$`1`$cov[1, 1]), tolerance = 1e-4)
  sem <- sqrt(diag(mirt::vcov(fmir))); sa <- sem[grep("^a1\\.", names(sem))][1]
  expect_equal(sqrt(diag(vcov(fm, "ML")))[["lambda0[y1,eta]"]], unname(sa), tolerance = 1e-3)
})

test_that("continuous items, continuous moderators, two factors: likelihood equals the closed form", {
  skip_on_cran()
  set.seed(4242); N <- 300
  x1 <- stats::rnorm(N); x2 <- stats::rbinom(N, 1, 0.5)
  r <- tanh(0.4 + 0.2 * x1); Z <- matrix(stats::rnorm(2 * N), N)
  e1 <- 0.3 * x1 - 0.2 * x2 + sqrt(exp(0.2 * x1)) * Z[, 1]
  e2 <- -0.2 * x1 + 0.4 * x2 + sqrt(exp(-0.3 * x2)) * (r * Z[, 1] + sqrt(1 - r^2) * Z[, 2])
  L <- c(1, 0.8, 1.2, 0.9, 1.1, 0.7); f <- c(1, 1, 1, 2, 2, 2)
  Y <- sapply(1:6, function(j) { et <- if (f[j] == 1) e1 else e2
    0.1 * j + 0.3 * (j == 5) * x1 + (L[j] + 0.25 * (j == 2) * x1) * et + stats::rnorm(N, sd = sqrt(0.5 * exp(0.3 * (j == 6) * x2))) })
  its <- paste0("y", 1:6); d <- data.frame(Y, x1, x2); names(d)[1:6] <- its
  m <- mofia_model(d, loadings = list(f1 = its[1:3], f2 = its[4:6]), moderators = c("x1", "x2"), continuous = its,
                   loading_mod = list(y2 = "x1"), threshold_mod = list(y5 = "x1"), resvar_mod = list(y6 = "x2"),
                   cor_mod = "x1")
  fit <- mofia(m)
  ll <- function(th) {
    g <- function(n) th[[n]]; X <- cbind(x1, x2); out <- 0
    mu <- cbind(X %*% c(g("gamma[f1,x1]"), g("gamma[f1,x2]")), X %*% c(g("gamma[f2,x1]"), g("gamma[f2,x2]")))
    s2 <- cbind(exp(X %*% c(g("delta[f1,x1]"), g("delta[f1,x2]"))), exp(X %*% c(g("delta[f2,x1]"), g("delta[f2,x2]"))))
    rho <- tanh(g("zcor0[f1,f2]") + x1 * g("zcor1[f1,f2,x1]"))
    for (i in seq_len(N)) {
      Lam <- matrix(0, 6, 2)
      for (j in 1:6) Lam[j, f[j]] <- g(sprintf("lambda0[%s,f%d]", its[j], f[j]))
      Lam[2, 1] <- Lam[2, 1] + x1[i] * g("lambda1[y2,f1,x1]")
      nu <- vapply(its, function(j) g(sprintf("nu0[%s]", j)), 1); nu[5] <- nu[5] + x1[i] * g("beta[y5,x1]")
      th_ <- vapply(its, function(j) exp(g(sprintf("omega0[%s]", j))), 1); th_[6] <- th_[6] * exp(x2[i] * g("omega[y6,x2]"))
      S <- diag(sqrt(s2[i, ])) %*% matrix(c(1, rho[i], rho[i], 1), 2) %*% diag(sqrt(s2[i, ]))
      V <- Lam %*% S %*% t(Lam) + diag(th_); C <- chol(V)
      z <- backsolve(C, Y[i, ] - nu - Lam %*% mu[i, ], transpose = TRUE)
      out <- out - 0.5 * (6 * log(2 * pi) + 2 * sum(log(diag(C))) + sum(z^2))
    }
    out
  }
  th <- coef(fit)
  expect_equal(fit$loglik, ll(th), tolerance = 1e-10)
  set.seed(1); th2 <- th + stats::rnorm(length(th), 0, 0.02)
  expect_equal(mofia_loglik(m, th2), ll(th2), tolerance = 1e-10)
})

test_that("ordinal items, continuous moderators: likelihood equals adaptive integrate()", {
  skip_on_cran()
  set.seed(99); N <- 150; x1 <- stats::rnorm(N); x2 <- stats::rbinom(N, 1, 0.5)
  eta <- stats::rnorm(N, 0.3 * x1 - 0.3 * x2, exp((0.25 * x1) / 2))
  lam <- c(1.4, 1.1, 0.9, 1.6, 1.2)
  Y <- sapply(1:5, function(j) findInterval((lam[j] + 0.3 * (j == 2) * x1) * eta + 0.4 * (j == 4) * x2 + stats::rlogis(N),
                                            c(-1, 0, 1) + 0.2 * j, left.open = TRUE))
  its <- paste0("y", 1:5); d <- data.frame(Y, x1, x2); names(d)[1:5] <- its
  for (link in c("logit", "probit")) {
    m <- mofia_model(d, items = its, moderators = c("x1", "x2"), loading_mod = list(y2 = "x1"),
                     threshold_mod = list(y4 = c("x1", "x2")), sparse = 0)
    fit <- mofia(m, link = link)
    Fl <- if (link == "logit") stats::plogis else stats::pnorm
    th <- coef(fit); g <- function(n) th[[n]]
    ll <- 0
    for (i in seq_len(N)) {
      mu <- g("gamma[eta,x1]") * x1[i] + g("gamma[eta,x2]") * x2[i]
      sd <- exp((g("delta[eta,x1]") * x1[i] + g("delta[eta,x2]") * x2[i]) / 2)
      fi <- function(e) {
        p <- rep(1, length(e))
        for (j in 1:5) {
          l <- g(sprintf("lambda0[%s,eta]", its[j])) + if (j == 2) g("lambda1[y2,eta,x1]") * x1[i] else 0
          a <- l * e + if (j == 4) g("beta[y4,x1]") * x1[i] + g("beta[y4,x2]") * x2[i] else 0
          te <- c(-Inf, th[sprintf("tau[%s,%d]", its[j], 1:3)], Inf); k <- m$Y[i, j]
          p <- p * (Fl(te[k + 2] - a) - Fl(te[k + 1] - a))
        }
        p * stats::dnorm(e, mu, sd)
      }
      ll <- ll + log(stats::integrate(fi, mu - 10 * sd, mu + 10 * sd, rel.tol = 1e-12, subdivisions = 500)$value)
    }
    # absolute difference: MoFIA's quadrature (31 adaptive nodes) against integrate()
    expect_lt(abs(fit$loglik - ll), 1e-5)
  }
})
