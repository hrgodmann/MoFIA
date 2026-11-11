# Finding 1 regression: raw Age (years) + Age^2 vs the standardized equivalent model.
.libPaths(c(commandArgs(TRUE)[1], .libPaths()))
suppressMessages(library(MoFIA))
data(parabelief)
ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession", "Sorcery", "HauntedPlaces")
d <- parabelief; d$AgeY <- round(48 + 17 * d$Age)
f_raw <- mofia(mofia_model(d, items = ub, moderators = ~ AgeY + I(AgeY^2) + Male, threshold_mod = "DevilPossession"))
f_std <- mofia(mofia_model(d, items = ub, moderators = ~ AgeY + I(AgeY^2) + Male, standardize = TRUE,
                           threshold_mod = "DevilPossession"))
cat("raw: status", f_raw$status, "converged", f_raw$converged, "raw eig", f_raw$information_eigen,
    "scaled eig", f_raw$information_eigen_scaled, "\n")
cat("std: status", f_std$status, "converged", f_std$converged, "\n")
cat(sprintf("logLik raw %.6f std %.6f diff %.2e\n", f_raw$loglik, f_std$loglik, f_raw$loglik - f_std$loglik))
m <- mean(d$AgeY); s <- sd(d$AgeY)
A <- "AgeY"; A2 <- "I(AgeY^2)"
g <- function(t) {
  o <- t; G <- function(n) t[[n]]
  cmu <- G("gamma[eta,AgeY]") * m + G("gamma[eta,I(AgeY^2)]") * m^2
  cv <- G("delta[eta,AgeY]") * m + G("delta[eta,I(AgeY^2)]") * m^2; k <- exp(cv / 2)
  o["gamma[eta,AgeY]"] <- (G("gamma[eta,AgeY]") + 2 * G("gamma[eta,I(AgeY^2)]") * m) * s / k
  o["gamma[eta,I(AgeY^2)]"] <- G("gamma[eta,I(AgeY^2)]") * s^2 / k
  o["gamma[eta,Male]"] <- G("gamma[eta,Male]") / k
  o["delta[eta,AgeY]"] <- (G("delta[eta,AgeY]") + 2 * G("delta[eta,I(AgeY^2)]") * m) * s
  o["delta[eta,I(AgeY^2)]"] <- G("delta[eta,I(AgeY^2)]") * s^2
  for (it in ub) {
    la <- G(sprintf("lambda0[%s,eta]", it)); o[sprintf("lambda0[%s,eta]", it)] <- la * k
    cb <- 0
    if (it == "DevilPossession") {
      b1 <- G("beta[DevilPossession,AgeY]"); b2 <- G("beta[DevilPossession,I(AgeY^2)]")
      cb <- b1 * m + b2 * m^2
      o["beta[DevilPossession,AgeY]"] <- (b1 + 2 * b2 * m) * s; o["beta[DevilPossession,I(AgeY^2)]"] <- b2 * s^2
    }
    tn <- grep(sprintf("^tau\\[%s,", it), names(t), value = TRUE)
    o[tn] <- t[tn] - cb - la * cmu
  }
  o
}
th_map <- g(coef(f_raw))
cat(sprintf("max |g(theta_raw) - theta_std| = %.2e\n", max(abs(th_map - coef(f_std)[names(th_map)]))))
J <- numDeriv::jacobian(g, coef(f_raw))
for (se in c("ML", "MLR")) {
  Vm <- J %*% vcov(f_raw, se) %*% t(J); Vs <- vcov(f_std, se)[names(th_map), names(th_map)]
  cat(sprintf("%s: max relative SE difference (J V_raw J' vs V_std) = %.2e\n", se, max(abs(sqrt(diag(Vm)) / sqrt(diag(Vs)) - 1))))
}
cat(sprintf("scaling factors raw %.6f std %.6f; Newton decrement raw %.2e std %.2e\n", f_raw$scaling[["exact"]], f_std$scaling[["exact"]],
            f_raw$newton_decrement, f_std$newton_decrement))
w1 <- mofia_wald(f_raw, "beta[DevilPossession,Male]"); w2 <- mofia_wald(f_std, "beta[DevilPossession,Male]")
cat(sprintf("Wald Male DIF raw %.6f std %.6f\n", w1$statistic, w2$statistic))
w1 <- mofia_wald(f_raw, "^beta\\[DevilPossession,", regex = TRUE); w2 <- mofia_wald(f_std, "^beta\\[DevilPossession,", regex = TRUE)
cat(sprintf("Wald all DIF raw %.6f std %.6f\n", w1$statistic, w2$statistic))
