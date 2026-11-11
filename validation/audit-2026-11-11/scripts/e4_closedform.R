# Independent implementation: continuous items, continuous + binary moderators, two factors.
# y_i ~ N(nu(x) + Lambda(x) mu(x), Lambda(x) Sigma(x) Lambda(x)' + Theta(x)). Fit by optim, compare with MoFIA.
suppressMessages(library(MoFIA))
set.seed(4242); N <- 500
x1 <- rnorm(N); x2 <- rbinom(N, 1, 0.5)
mu1 <- 0.3 * x1 - 0.2 * x2; mu2 <- -0.2 * x1 + 0.4 * x2
v1 <- exp(0.2 * x1); v2 <- exp(-0.3 * x2); r <- tanh(0.4 + 0.2 * x1)
Z <- matrix(rnorm(2 * N), N); e1 <- mu1 + sqrt(v1) * Z[, 1]; e2 <- mu2 + sqrt(v2) * (r * Z[, 1] + sqrt(1 - r^2) * Z[, 2])
L <- c(1, 0.8, 1.2, 0.9, 1.1, 0.7); f <- c(1, 1, 1, 2, 2, 2)
lamx <- c(0, 0.25, 0, 0, 0, 0)       # loading of y2 moderated by x1
Y <- sapply(1:6, function(j) { et <- if (f[j] == 1) e1 else e2
  0.1 * j + 0.3 * (j == 5) * x1 + (L[j] + lamx[j] * x1) * et + rnorm(N, sd = sqrt(0.5 * exp(0.3 * (j == 6) * x2))) })
its <- paste0("y", 1:6); d <- data.frame(Y, x1, x2); names(d)[1:6] <- its
m <- mofia_model(d, loadings = list(f1 = its[1:3], f2 = its[4:6]), moderators = c("x1", "x2"), continuous = its,
                 loading_mod = list(y2 = "x1"), threshold_mod = list(y5 = "x1"), resvar_mod = list(y6 = "x2"),
                 cor_mod = "x1")
fit <- mofia(m)
pt <- m$par
nll <- function(th) {
  g <- function(n) th[[n]]
  X <- cbind(x1, x2); ll <- 0
  mu <- cbind(X %*% c(g("gamma[f1,x1]"), g("gamma[f1,x2]")), X %*% c(g("gamma[f2,x1]"), g("gamma[f2,x2]")))
  s2 <- cbind(exp(X %*% c(g("delta[f1,x1]"), g("delta[f1,x2]"))), exp(X %*% c(g("delta[f2,x1]"), g("delta[f2,x2]"))))
  rho <- tanh(g("zcor0[f1,f2]") + x1 * g("zcor1[f1,f2,x1]"))
  out <- 0
  for (i in seq_len(N)) {
    Lam <- matrix(0, 6, 2)
    for (j in 1:6) Lam[j, f[j]] <- g(sprintf("lambda0[%s,f%d]", its[j], f[j]))
    Lam[2, 1] <- Lam[2, 1] + x1[i] * g("lambda1[y2,f1,x1]")
    nu <- sapply(its, function(j) g(sprintf("nu0[%s]", j))); nu[5] <- nu[5] + x1[i] * g("beta[y5,x1]")
    th_ <- sapply(its, function(j) exp(g(sprintf("omega0[%s]", j)))); th_[6] <- th_[6] * exp(x2[i] * g("omega[y6,x2]"))
    S <- diag(sqrt(s2[i, ])) %*% matrix(c(1, rho[i], rho[i], 1), 2) %*% diag(sqrt(s2[i, ]))
    V <- Lam %*% S %*% t(Lam) + diag(th_); mm <- nu + Lam %*% mu[i, ]
    C <- chol(V); z <- backsolve(C, Y[i, ] - mm, transpose = TRUE)
    out <- out - 0.5 * (6 * log(2 * pi) + 2 * sum(log(diag(C))) + sum(z^2))
  }
  -out
}
th0 <- coef(fit)
cat(sprintf("closed-form logLik at MoFIA estimates: %.8f | MoFIA: %.8f | diff %.2e\n", -nll(th0), fit$loglik, -nll(th0) - fit$loglik))
nll2 <- function(t) tryCatch(nll(t), error = function(e) 1e10); op <- optim(th0 + rnorm(length(th0), 0, 0.02), nll2, method = "BFGS", control = list(maxit = 2000, reltol = 1e-14))
cat(sprintf("independent optim from perturbed start: logLik %.6f, max |theta diff| vs MoFIA = %.2e\n",
            -op$value, max(abs(op$par - th0))))
H <- numDeriv::hessian(nll, th0)
se_cf <- sqrt(diag(solve(H))); se_mf <- sqrt(diag(vcov(fit, "ML")))
cat(sprintf("ML SEs: max relative difference closed-form vs MoFIA = %.2e\n", max(abs(se_cf / se_mf - 1))))
