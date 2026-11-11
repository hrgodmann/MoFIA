# Ordinal items, continuous + binary moderators acting on loadings, thresholds, mean, variance:
# per-person likelihood by adaptive integrate() (independent of MoFIA's quadrature), at MoFIA's estimates,
# plus the gradient of that independent likelihood (should be ~0 at the MLE).
suppressMessages(library(MoFIA))
set.seed(99); N <- 400; x1 <- rnorm(N); x2 <- rbinom(N, 1, 0.5)
eta <- rnorm(N, 0.3 * x1 - 0.3 * x2, exp((0.25 * x1) / 2))
tau <- c(-1, 0, 1); lam <- c(1.4, 1.1, 0.9, 1.6, 1.2)
Y <- sapply(1:5, function(j) findInterval((lam[j] + 0.3 * (j == 2) * x1) * eta + 0.4 * (j == 4) * x2 + rlogis(N),
                                          tau + 0.2 * j, left.open = TRUE))
its <- paste0("y", 1:5); d <- data.frame(Y, x1, x2); names(d)[1:5] <- its
for (link in c("logit", "probit")) {
  m <- mofia_model(d, items = its, moderators = c("x1", "x2"), loading_mod = list(y2 = "x1"),
                   threshold_mod = list(y4 = c("x1", "x2")), sparse = 0)
  fit <- mofia(m, link = link)
  F <- if (link == "logit") plogis else pnorm
  ll_ind <- function(th) {
    g <- function(n) th[[n]]
    out <- 0
    for (i in seq_len(N)) {
      mu <- g("gamma[eta,x1]") * x1[i] + g("gamma[eta,x2]") * x2[i]
      sd <- exp((g("delta[eta,x1]") * x1[i] + g("delta[eta,x2]") * x2[i]) / 2)
      f <- function(e) {
        p <- rep(1, length(e))
        for (j in 1:5) {
          it <- its[j]; l <- g(sprintf("lambda0[%s,eta]", it)) + if (j == 2) g("lambda1[y2,eta,x1]") * x1[i] else 0
          a <- l * e + if (j == 4) g("beta[y4,x1]") * x1[i] + g("beta[y4,x2]") * x2[i] else 0
          te <- c(-Inf, th[sprintf("tau[%s,%d]", it, 1:3)], Inf); k <- Y[i, j]
          p <- p * (F(te[k + 2] - a) - F(te[k + 1] - a))
        }
        p * dnorm(e, mu, sd)
      }
      out <- out + log(integrate(f, mu - 10 * sd, mu + 10 * sd, rel.tol = 1e-12, subdivisions = 500)$value)
    }
    out
  }
  th <- coef(fit)
  li <- ll_ind(th)
  gi <- numDeriv::grad(ll_ind, th, method.args = list(r = 4))
  cat(sprintf("%s: nq = %d  logLik MoFIA %.7f  integrate() %.7f  diff %.1e  | max |grad| of independent logLik at MoFIA MLE %.1e\n",
              link, fit$nq, fit$loglik, li, fit$loglik - li, max(abs(gi))))
}
