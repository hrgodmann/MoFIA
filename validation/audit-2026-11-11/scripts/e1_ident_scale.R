suppressMessages(library(MoFIA))
data(parabelief)
ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession", "Sorcery", "HauntedPlaces")
cat("=== E1: binary moderator; mean + all thresholds + all loadings moderated by Male, var fixed ===\n")
m <- withCallingHandlers(
  mofia_model(parabelief, items = ub, moderators = "Male", mean_mod = TRUE, var_mod = FALSE,
              threshold_mod = TRUE, loading_mod = TRUE),
  warning = function(w) { cat("model warning:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
f <- withCallingHandlers(mofia(m),
  warning = function(w) { cat("fit warning:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
cat("status:", f$status, " converged:", f$converged, " eig:", f$information_eigen, "\n")
# demonstrate the exact non-identification: shift gamma by c, compensate beta_j by -(lambda0_j + lambda1_j) c
th <- coef(f); c0 <- 0.7; th2 <- th
th2["gamma[eta,Male]"] <- th["gamma[eta,Male]"] + c0
for (it in ub) {
  l <- th[sprintf("lambda0[%s,eta]", it)] + th[sprintf("lambda1[%s,eta,Male]", it)]
  th2[sprintf("beta[%s,Male]", it)] <- th[sprintf("beta[%s,Male]", it)] - l * c0
}
cat("logLik at estimate:", mofia_loglik(f), "  after shifting gamma by 0.7 and compensating:", mofia_loglik(f, th2), "\n")

cat("\n=== E2: moderator scale. Age in years + Age^2 vs standardized ===\n")
d <- parabelief; d$AgeY <- round(48 + 17 * d$Age)
run <- function(mods, std = FALSE) {
  w <- character(0)
  f <- withCallingHandlers(mofia(mofia_model(d, items = ub, moderators = mods, standardize = std,
                                             threshold_mod = "DevilPossession")),
                           warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  ev <- f$information_eigen
  cat(sprintf("%-32s logLik %.5f  status %-8s converged %s  eig ratio %.2e\n", deparse(mods), f$loglik, f$status,
              f$converged, ev[1] / ev[2]))
  if (length(w)) cat("   warnings:", paste(substr(w, 1, 140), collapse = " | "), "\n")
  f
}
fa <- run(~ AgeY + I(AgeY^2) + Male)
fb <- run(~ AgeY + I(AgeY^2) + Male, std = TRUE)
fc <- run(~ AgeY + Male)
fd <- run(~ AgeY + Male, std = TRUE)
# invariant quantity: the DIF Wald test of beta[DevilPossession, Male] should not depend on the scaling of Age
for (ff in list(fc, fd)) if (!is.null(ff$vcov)) {
  print(mofia_wald(ff, "beta[DevilPossession,Male]")[, 1:3])
}
