# Simulated two-factor data with loading, threshold and structural moderation.
sim_two_factor <- function(N = 400, seed = 1, continuous = FALSE) {
  set.seed(seed)
  X <- data.frame(x1 = stats::rbinom(N, 1, 0.5), x2 = round(stats::rnorm(N), 4))
  items <- paste0("y", 1:6)
  ph <- data.frame(matrix(rep(0:3, length.out = N * 6), N, 6, dimnames = list(NULL, items)), X)
  if (continuous) ph$y6 <- stats::rnorm(N)
  args <- list(loadings = list(f1 = items[1:3], f2 = items[4:6]), moderators = c("x1", "x2"),
               loading_mod = list(y2 = "x1"), threshold_mod = list(y3 = "x2"), cor_mod = TRUE,
               continuous = if (continuous) "y6" else NULL, resvar_mod = if (continuous) list(y6 = "x1") else NULL)
  m0 <- do.call(mofia_model, c(list(data = ph), args))
  th <- stats::setNames(numeric(nrow(m0$par)), m0$par$name)
  th[grepl("^tau", names(th))] <- rep(c(-1.2, 0, 1.2), length.out = sum(grepl("^tau", names(th))))
  th[grepl("^lambda0", names(th))] <- c(1.5, 1.2, 1.0, 1.4, 1.1, 1.3)
  th["lambda1[y2,f1,x1]"] <- -0.4; th["beta[y3,x2]"] <- 0.5
  th[c("gamma[f1,x1]", "gamma[f1,x2]", "gamma[f2,x1]", "gamma[f2,x2]")] <- c(0.3, -0.2, -0.2, 0.1)
  th[c("delta[f1,x1]", "delta[f1,x2]", "delta[f2,x1]", "delta[f2,x2]")] <- c(0.2, 0.1, -0.2, 0.1)
  th[c("zcor0[f1,f2]", "zcor1[f1,f2,x1]", "zcor1[f1,f2,x2]")] <- c(0.4, 0.2, -0.1)
  if (continuous) th[c("nu0[y6]", "omega0[y6]", "omega[y6,x1]")] <- c(0.5, log(0.6), 0.3)
  dat <- mofia_simulate(m0, th, newdata = X, seed = seed + 1)
  list(data = dat, theta = th, args = args)
}
