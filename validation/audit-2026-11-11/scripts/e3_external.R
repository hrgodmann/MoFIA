suppressMessages({ library(MoFIA); library(lavaan); library(mirt) })
set.seed(2026)
cat("=== E3a: continuous items, binary moderator vs lavaan multigroup (ML, observed Hessian information) ===\n")
N <- 600; g <- rep(0:1, each = N / 2)
eta <- rnorm(N, 0.4 * g, exp(0.3 * g / 2))
lam <- c(1, 0.8, 1.2, 0.9, 1.1); nu <- c(0, 0.3, -0.2, 0.1, 0.5)
Y <- sapply(1:5, function(j) nu[j] + lam[j] * eta + 0.3 * (j == 4) * g + rnorm(N, sd = 0.7))
d <- data.frame(Y, g = g); its <- paste0("y", 1:5); names(d)[1:5] <- its
# MoFIA: invariant loadings + residual variances; intercept DIF for y4; group-specific factor mean/var
fm <- mofia(mofia_model(d, items = its, moderators = "g", continuous = its, threshold_mod = "y4"))
# lavaan: equivalent multigroup model (group 0 = reference: mean 0, var 1)
mod <- "f =~ NA*y1 + y2 + y3 + y4 + y5
        f ~~ c(1, NA)*f
        f ~ c(0, NA)*1"
fl <- cfa(mod, d, group = "g", group.label = c("0", "1"), group.equal = c("loadings", "intercepts", "residuals"),
          group.partial = "y4~1", estimator = "ML", information = "observed", observed.information = "hessian",
          meanstructure = TRUE)
cat(sprintf("logLik  MoFIA %.6f   lavaan %.6f   diff %.2e\n", fm$loglik, logLik(fl), fm$loglik - as.numeric(logLik(fl))))
pe <- parameterEstimates(fl)
cmp <- data.frame(
  what = c("lambda y2", "beta y4 (intercept DIF)", "gamma (mean g=1)", "log var g=1"),
  mofia = c(coef(fm)["lambda0[y2,eta]"], coef(fm)["beta[y4,g]"], coef(fm)["gamma[eta,g]"], coef(fm)["delta[eta,g]"]),
  lavaan = c(pe$est[pe$lhs == "f" & pe$op == "=~" & pe$rhs == "y2" & pe$group == 1],
             pe$est[pe$lhs == "y4" & pe$op == "~1" & pe$group == 2] - pe$est[pe$lhs == "y4" & pe$op == "~1" & pe$group == 1],
             pe$est[pe$lhs == "f" & pe$op == "~1" & pe$group == 2],
             log(pe$est[pe$lhs == "f" & pe$op == "~~" & pe$group == 2])))
print(cmp, digits = 7)
# SEs for directly comparable parameters (loading, factor mean)
se_m <- sqrt(diag(vcov(fm, "ML")))
cat(sprintf("ML SE lambda y2: MoFIA %.6f lavaan %.6f | factor mean: MoFIA %.6f lavaan %.6f\n",
            se_m["lambda0[y2,eta]"], pe$se[pe$lhs == "f" & pe$op == "=~" & pe$rhs == "y2" & pe$group == 1],
            se_m["gamma[eta,g]"], pe$se[pe$lhs == "f" & pe$op == "~1" & pe$group == 2]))
flr <- update(fl, estimator = "MLR")
per <- parameterEstimates(flr); se_r <- sqrt(diag(vcov(fm, "MLR")))
cat(sprintf("MLR SE lambda y2: MoFIA %.6f lavaan %.6f | factor mean: MoFIA %.6f lavaan %.6f\n",
            se_r["lambda0[y2,eta]"], per$se[per$lhs == "f" & per$op == "=~" & per$rhs == "y2" & per$group == 1],
            se_r["gamma[eta,g]"], per$se[per$lhs == "f" & per$op == "~1" & per$group == 2]))
cat(sprintf("MLR scaling factor (MoFIA exact info): %.4f\n", fm$scaling[["exact"]]))

cat("\n=== E3b: ordinal items (logit graded), binary moderator vs mirt multipleGroup ===\n")
set.seed(77); N <- 1000; g <- rep(0:1, each = N / 2)
eta <- rnorm(N, -0.3 * g, exp(-0.4 * g / 2))
a <- c(1.5, 1.2, 1.0, 1.8, 0.9, 1.3); tau <- rbind(c(-1, 0.2, 1.3), c(-0.5, 0.5, 1.5), c(-1.5, -0.3, 0.8),
                                                    c(-0.8, 0.4, 1.1), c(-1.2, 0, 1.4), c(-0.6, 0.6, 1.6))
Y <- sapply(1:6, function(j) findInterval(a[j] * eta + 0.5 * (j == 3) * g + rlogis(N), tau[j, ], left.open = TRUE))
d <- data.frame(Y, g = g); its <- paste0("y", 1:6); names(d)[1:6] <- its
fm <- mofia(mofia_model(d, items = its, moderators = "g", sparse = 0))
fmir <- multipleGroup(d[, its], 1, group = factor(d$g), itemtype = "graded",
                      invariance = c("slopes", "intercepts", "free_means", "free_var"),
                      SE = TRUE, SE.type = "Oakes", verbose = FALSE, technical = list(NCYCLES = 5000),
                      TOL = 1e-8, quadpts = 121)
cat(sprintf("logLik  MoFIA %.5f   mirt %.5f   diff %.2e\n", fm$loglik, logLik(fmir), fm$loglik - logLik(fmir)))
cm <- coef(fmir, simplify = TRUE)
cat(sprintf("slope y1: MoFIA %.5f mirt %.5f | mean g1: MoFIA %.5f mirt %.5f | var g1: MoFIA %.5f mirt %.5f\n",
            coef(fm)["lambda0[y1,eta]"], cm$`0`$items["y1", "a1"], coef(fm)["gamma[eta,g]"], cm$`1`$means,
            exp(coef(fm)["delta[eta,g]"]), cm$`1`$cov))
cat(sprintf("tau[y1,] MoFIA: %s | -d mirt: %s\n", paste(round(coef(fm)[sprintf("tau[y1,%d]", 1:3)], 5), collapse = " "),
            paste(round(-cm$`0`$items["y1", c("d1", "d2", "d3")], 5), collapse = " ")))
se_mirt <- sqrt(diag(vcov(fmir)))
cat("mirt SE names (head):", paste(head(names(se_mirt), 3), collapse = ", "), "\n")
cat(sprintf("ML SE slope y1 MoFIA %.5f | mirt (Oakes) a1 SEs: %s\n", sqrt(diag(vcov(fm, "ML")))["lambda0[y1,eta]"],
            paste(round(se_mirt[grep("^a1\\.1$|a1", names(se_mirt))][1:2], 5), collapse = " ")))
# DIF: item 3 uniform threshold shift by g (MoFIA beta) vs mirt with free intercepts for y3 (category-specific)
fm3 <- mofia(mofia_model(d, items = its, moderators = "g", threshold_mod = "y3", sparse = 0))
cat(sprintf("MoFIA uniform DIF y3: beta = %.3f (true 0.5), LR vs invariant = %.2f on 1 df\n",
            coef(fm3)["beta[y3,g]"], 2 * (fm3$loglik - fm$loglik)))
