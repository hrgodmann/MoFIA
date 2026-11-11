# Calibration: CI coverage (ML, MLR) and type-I error of Wald / scaled LRT / score tests for a null DIF effect.
# True model: 1 factor, 6 ordinal (4-cat) items, x1 continuous, x2 binary; mean+var moderated, loading of y2 by x1.
# Fitted: M0 = true structure; M1 = M0 + uniform DIF of y4 on x1 and x2 (true value 0).
suppressMessages({ library(MoFIA); library(parallel) })
R <- as.integer(commandArgs(TRUE)[1]); N <- as.integer(commandArgs(TRUE)[2])
its <- paste0("y", 1:6)
one <- function(r) {
  set.seed(1000 + r)
  x1 <- rnorm(N); x2 <- rbinom(N, 1, 0.5)
  eta <- rnorm(N, 0.3 * x1 - 0.4 * x2, exp((0.2 * x1 + 0.3 * x2) / 2))
  lam <- c(1.5, 1.2, 1.0, 1.8, 0.9, 1.3)
  Y <- sapply(1:6, function(j) findInterval((lam[j] + 0.3 * (j == 2) * x1) * eta + rlogis(N),
                                            c(-1, 0.1, 1.2) + 0.15 * (j - 3), left.open = TRUE))
  d <- data.frame(Y, x1, x2); names(d)[1:6] <- its
  truth <- c("lambda0[y1,eta]" = 1.5, "lambda0[y2,eta]" = 1.2, "lambda1[y2,eta,x1]" = 0.3,
             "gamma[eta,x1]" = 0.3, "gamma[eta,x2]" = -0.4, "delta[eta,x1]" = 0.2, "delta[eta,x2]" = 0.3,
             "beta[y4,x1]" = 0, "beta[y4,x2]" = 0)
  res <- tryCatch({
    m0 <- mofia_model(d, items = its, moderators = c("x1", "x2"), loading_mod = list(y2 = "x1"), sparse = 0)
    m1 <- mofia_model(d, items = its, moderators = c("x1", "x2"), loading_mod = list(y2 = "x1"),
                      threshold_mod = list(y4 = c("x1", "x2")), sparse = 0)
    f0 <- suppressWarnings(mofia(m0)); f1 <- suppressWarnings(mofia(m1, start = coef(f0)))
    cov <- sapply(c("ML", "MLR"), function(s) {
      ci <- confint(f1, names(truth), se = s); as.numeric(ci[, 1] <= truth & truth <= ci[, 2]) })
    est <- coef(f1)[names(truth)]
    w <- mofia_wald(f1, c("beta[y4,x1]", "beta[y4,x2]"))$p_value
    lr <- mofia_lrt(f0, f1)$p_value; lru <- mofia_lrt(f0, f1, scaled = FALSE)$p_value
    sc <- mofia_score_test(f0, c("beta[y4,x1]", "beta[y4,x2]"), by = "item")
    list(ok = TRUE, cov_ML = cov[, "ML"], cov_MLR = cov[, "MLR"], est = est,
         p = c(wald = w, lrt_scaled = lr, lrt = lru, score_MLR = sc$p, score_ML = pchisq(sc$statistic_ML, 2, lower.tail = FALSE)),
         c0 = f0$scaling[["exact"]], conv = c(f0$converged, f1$converged))
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  res
}
t0 <- proc.time()[[3]]
out <- mclapply(seq_len(R), one, mc.cores = 8)
ok <- vapply(out, `[[`, TRUE, "ok"); cat(sprintf("N = %d, %d/%d replications ok, %.0f s\n", N, sum(ok), R, proc.time()[[3]] - t0))
if (any(!ok)) print(table(vapply(out[!ok], `[[`, "", "msg")))
o <- out[ok]
truth <- c(1.5, 1.2, 0.3, 0.3, -0.4, 0.2, 0.3, 0, 0)
E <- t(sapply(o, `[[`, "est"))
tab <- data.frame(truth = truth, mean_est = colMeans(E), rel_bias = ifelse(truth != 0, colMeans(E) / truth - 1, NA),
                  emp_sd = apply(E, 2, sd),
                  cover_ML = rowMeans(sapply(o, `[[`, "cov_ML")), cover_MLR = rowMeans(sapply(o, `[[`, "cov_MLR")))
print(round(tab, 3))
P <- t(sapply(o, `[[`, "p"))
cat("Rejection rate at alpha = .05 (null true; MC SE ~", round(sqrt(.05 * .95 / nrow(P)), 3), "):\n")
print(round(colMeans(P < 0.05), 3))
cat("mean scaling factor c0:", round(mean(sapply(o, `[[`, "c0")), 3), " all converged:", all(sapply(o, function(x) all(x$conv))), "\n")
saveRDS(out, sprintf("calib_N%d_R%d.rds", N, R))
