# Likelihood: marginal over the factors by (adaptive) product Gauss-Hermite
# quadrature. The C++ engine (src/mofia_core.cpp) handles any number of factors;
# the pure-R reference engine (one or two factors) is kept for testing.

.link_funs <- function(link = c("logit", "probit")) {
  link <- match.arg(link)
  if (link == "logit") list(F = stats::plogis, f = stats::dlogis, name = "logit")
  else                 list(F = stats::pnorm,  f = stats::dnorm,  name = "probit")
}

# product Gauss-Hermite rule for N(0, I_D)
.mofia_nodes <- function(nq = 31, D = 1) {
  gq <- statmod::gauss.quad.prob(nq, dist = "normal")
  if (D == 1) return(list(T = matrix(gq$nodes, ncol = 1), logw = log(gq$weights), D = 1, nq = nq))
  g <- expand.grid(rep(list(seq_len(nq)), D))
  T <- sapply(seq_len(D), function(d) gq$nodes[g[[d]]])
  logw <- rowSums(sapply(seq_len(D), function(d) log(gq$weights[g[[d]]])))
  list(T = T, logw = logw, D = D, nq = nq)
}

# default number of quadrature nodes per dimension
.default_nq <- function(D) c(31, 15, 9, 7, 5, 4)[min(D, 6)]

.mofia_threads <- function() max(1L, as.integer(getOption("mofia.threads", 1L)))

# log-likelihood dispatcher; adapt: NULL or list(m = N x D posterior means of z,
# Lc = N x D(D+1)/2 row-wise lower triangle of the posterior Cholesky factor)
.mofia_ll <- function(theta, spec, link = "logit", nodes = .mofia_nodes(31, spec$D),
                      scores = FALSE, per_person = FALSE, adapt = NULL, return_post = FALSE,
                      engine = getOption("mofia.engine", "cpp")) {
  f <- if (engine == "cpp") .mofia_ll_cpp else .mofia_ll_R
  f(theta, spec, link, nodes, scores, per_person, adapt, return_post)
}

# ---- reference engine (pure R, vectorised over persons x nodes) -----------------
.mofia_ll_R <- function(theta, spec, link = "logit", nodes = .mofia_nodes(31, spec$D),
                            scores = FALSE, per_person = FALSE, adapt = NULL, return_post = FALSE) {
  if (spec$D > 2) stop("The R reference engine supports 1 or 2 factors; use the C++ engine.")
  if (!is.null(adapt) && is.null(adapt$L11)) adapt <- .adapt_compat(adapt, spec$D)
  lf <- .link_funs(link); Fc <- lf$F; fd <- lf$f
  pr <- .unpack(theta, spec); X <- spec$X; Y <- spec$Y
  N <- spec$N; P <- spec$P; D <- spec$D; M <- nrow(nodes$T)
  for (j in seq_len(P)) if (is.unsorted(pr$tau[[j]], strictly = TRUE)) return(-Inf)
  st <- .structural(pr, spec); if (!st$ok) return(-Inf)

  Tm <- lapply(seq_len(D), function(d) matrix(nodes$T[, d], N, M, byrow = TRUE))
  LW <- matrix(nodes$logw, N, M, byrow = TRUE)
  if (is.null(adapt)) {
    Z <- Tm
  } else {
    if (D == 1) {
      Z <- list(adapt$m[, 1] + adapt$L11 * Tm[[1]])
      LW <- LW + log(adapt$L11)
    } else {
      Z <- list(adapt$m[, 1] + adapt$L11 * Tm[[1]],
                adapt$m[, 2] + adapt$L21 * Tm[[1]] + adapt$L22 * Tm[[2]])
      LW <- LW + log(adapt$L11) + log(adapt$L22)
    }
    for (d in seq_len(D)) LW <- LW + stats::dnorm(Z[[d]], log = TRUE) - stats::dnorm(Tm[[d]], log = TRUE)
  }
  E <- vector("list", D)
  E[[1]] <- st$mu[, 1] + st$s[, 1] * Z[[1]]
  if (D == 2) { sq <- sqrt(1 - st$rho[, 1]^2)
                E[[2]] <- st$mu[, 2] + st$s[, 2] * (st$rho[, 1] * Z[[1]] + sq * Z[[2]]) }

  logLq <- matrix(0, N, M); keep <- vector("list", P)
  iscont <- if (is.null(spec$cont)) rep(FALSE, P) else spec$cont
  for (j in seq_len(P)) {
    a <- drop(X %*% pr$B[j, ]) + pr$nu0[j]
    lam <- vector("list", D)
    for (d in seq_len(D)) if (spec$loadings[j, d]) {
      lam[[d]] <- pr$L0[j, d] + drop(X %*% pr$L1[j, d, ])
      a <- a + lam[[d]] * E[[d]]
    }
    if (iscont[j]) {
      th <- exp(pr$Om0[j] + drop(X %*% pr$Om[j, ]))
      y <- Y[, j]; obs <- !is.na(y); r <- ifelse(obs, y, 0) - a
      lp <- -0.5 * log(2 * pi * th) - r^2 / (2 * th); lp[!obs, ] <- 0
      logLq <- logLq + lp
      if (scores) keep[[j]] <- list(cont = TRUE, lam = lam, r = r, th = th, obs = obs)
      next
    }
    te <- c(-Inf, pr$tau[[j]], Inf)
    y <- Y[, j]; obs <- !is.na(y); yy <- ifelse(obs, y, 0L)
    lo <- te[yy + 1]; hi <- te[yy + 2]
    p <- Fc(hi - a) - Fc(lo - a); p[!obs, ] <- 1; p <- pmax(p, 1e-300)
    logLq <- logLq + log(p)
    if (scores) keep[[j]] <- list(lam = lam, a = a, lo = lo, hi = hi, p = p, obs = obs, yy = yy)
  }
  lw <- logLq + LW
  mx <- apply(lw, 1, max); lli <- mx + log(rowSums(exp(lw - mx)))
  out <- if (per_person) lli else sum(lli)
  post <- exp(lw - lli)
  if (return_post) {
    m <- sapply(seq_len(D), function(d) rowSums(post * Z[[d]])); m <- matrix(m, N, D)
    v11 <- pmax(rowSums(post * Z[[1]]^2) - m[, 1]^2, 1e-10)
    ad <- list(m = m, L11 = sqrt(v11))
    if (D == 2) {
      c21 <- rowSums(post * Z[[1]] * Z[[2]]) - m[, 1] * m[, 2]
      v22 <- rowSums(post * Z[[2]]^2) - m[, 2]^2
      ad$L21 <- c21 / ad$L11; ad$L22 <- sqrt(pmax(v22 - ad$L21^2, 1e-10))
    }
    attr(out, "adapt") <- ad
  }
  if (!scores) return(out)

  pt <- spec$par; S <- matrix(0, N, nrow(pt), dimnames = list(NULL, pt$name))
  A <- lapply(seq_len(D), function(d) matrix(0, N, M))        # d l / d eta_d at each node
  PH <- matrix(0, N, P)                                       # d l_i / d log theta_ij
  for (j in seq_len(P)) {
    k <- keep[[j]]
    if (isTRUE(k$cont)) {
      g <- k$r / k$th; g[!k$obs, ] <- 0
      h <- -0.5 + k$r^2 / (2 * k$th); h[!k$obs, ] <- 0
      pg <- rowSums(post * g); PH[, j] <- rowSums(post * h)
      for (d in seq_len(D)) if (spec$loadings[j, d]) {
        pge <- rowSums(post * g * E[[d]])
        S[, pt$block == "lambda0" & pt$item == j & pt$factor %in% d] <- pge
        for (c in which(pt$block == "lambda1" & pt$item == j & pt$factor %in% d)) S[, c] <- X[, pt$mod[c]] * pge
        A[[d]] <- A[[d]] + k$lam[[d]] * g
      }
      for (c in which(pt$block == "beta" & pt$item == j)) S[, c] <- X[, pt$mod[c]] * pg
      S[, pt$block == "nu0" & pt$item == j] <- pg
      S[, pt$block == "omega0" & pt$item == j] <- PH[, j]
      for (c in which(pt$block == "omega" & pt$item == j)) S[, c] <- X[, pt$mod[c]] * PH[, j]
      next
    }
    flo <- fd(k$lo - k$a); fhi <- fd(k$hi - k$a)
    g <- (flo - fhi) / k$p; g[!k$obs, ] <- 0
    pg <- rowSums(post * g)
    for (d in seq_len(D)) if (spec$loadings[j, d]) {
      pge <- rowSums(post * g * E[[d]])
      S[, pt$block == "lambda0" & pt$item == j & pt$factor %in% d] <- pge
      for (c in which(pt$block == "lambda1" & pt$item == j & pt$factor %in% d)) S[, c] <- X[, pt$mod[c]] * pge
      A[[d]] <- A[[d]] + k$lam[[d]] * g
    }
    for (c in which(pt$block == "beta" & pt$item == j)) S[, c] <- X[, pt$mod[c]] * pg
    dlo <- rowSums(post * (-flo / k$p)); dhi <- rowSums(post * (fhi / k$p))
    dlo[!k$obs] <- 0; dhi[!k$obs] <- 0
    tcols <- which(pt$block == "tau" & pt$item == j)
    for (kk in seq_along(tcols)) S[, tcols[kk]] <- ifelse(k$yy == kk, dlo, 0) + ifelse(k$yy + 1 == kk, dhi, 0)
  }
  dmu <- sapply(seq_len(D), function(d) rowSums(post * A[[d]])); dmu <- matrix(dmu, N, D)
  for (c in which(pt$block == "gamma")) S[, c] <- X[, pt$mod[c]] * dmu[, pt$factor[c]]
  ds <- matrix(0, N, D)
  ds[, 1] <- rowSums(post * A[[1]] * Z[[1]])
  if (D == 2) {
    rho <- st$rho[, 1]; s1 <- st$s[, 1]; s2 <- st$s[, 2]
    ds[, 2] <- rowSums(post * A[[2]] * (rho * Z[[1]] + sq * Z[[2]]))
    dr <- rowSums(post * A[[2]] * s2 * (Z[[1]] - rho / sq * Z[[2]]))   # d l / d rho
    if (spec$cov_param == "fisherz") {
      du <- dr * (1 - rho^2)
      S[, pt$block == "r0"] <- du
      for (c in which(pt$block == "r1")) S[, c] <- X[, pt$mod[c]] * du
    } else {
      ds[, 1] <- ds[, 1] - dr * rho / s1
      ds[, 2] <- ds[, 2] - dr * rho / s2
      dc <- dr / (s1 * s2)
      S[, pt$block == "r0"] <- dc
      for (c in which(pt$block == "r1")) S[, c] <- X[, pt$mod[c]] * dc
    }
  }
  for (c in which(pt$block == "delta")) {
    d <- pt$factor[c]; S[, c] <- X[, pt$mod[c]] * ds[, d] * st$s[, d] / 2
  }
  for (c in which(pt$block == "delta0")) { d <- pt$factor[c]; S[, c] <- ds[, d] * st$s[, d] / 2 }
  attr(out, "scores") <- S
  # derivatives w.r.t. the person-specific factor (residual) variances and the
  # factor covariance, each holding the others fixed (the model's person-level parameters):
  # W[, d] = d l_i / d log vf_d,i ;  W[, "cov"] = d l_i / d cov12_i
  if (D == 1) {
    W <- cbind(logvf1 = ds[, 1] * st$s[, 1] / 2)
  } else {
    ds1r <- rowSums(post * A[[1]] * Z[[1]])
    ds2r <- rowSums(post * A[[2]] * (rho * Z[[1]] + sq * Z[[2]]))
    W <- cbind(logvf1 = (ds1r - dr * rho / s1) * s1 / 2, logvf2 = (ds2r - dr * rho / s2) * s2 / 2,
               cov = dr / (s1 * s2))
  }
  attr(out, "phi_scores") <- W
  attr(out, "resvar_scores") <- PH
  out
}

# ---- C++ engine ---------------------------------------------------------------
.mofia_ll_cpp <- function(theta, spec, link = "logit", nodes = .mofia_nodes(31, spec$D),
                              scores = FALSE, per_person = FALSE, adapt = NULL, return_post = FALSE) {
  pr <- .unpack(theta, spec); X <- spec$X; N <- spec$N; P <- spec$P; D <- spec$D; Tn <- D * (D + 1) / 2
  for (j in seq_len(P)) if (is.unsorted(pr$tau[[j]], strictly = TRUE)) return(-Inf)
  st <- .structural(pr, spec); if (!st$ok) return(-Inf)
  NU <- X %*% t(pr$B) + matrix(pr$nu0, N, P, byrow = TRUE)
  cont <- if (is.null(spec$cont)) rep(FALSE, P) else spec$cont
  RV <- matrix(1, N, P); YC <- matrix(NaN, N, P)
  for (j in which(cont)) { RV[, j] <- exp(pr$Om0[j] + drop(X %*% pr$Om[j, ])); YC[, j] <- ifelse(is.na(spec$Y[, j]), NaN, spec$Y[, j]) }
  LAM <- matrix(0, N, P * D)
  for (d in seq_len(D)) for (j in seq_len(P)) if (spec$loadings[j, d])
    LAM[, j + P * (d - 1)] <- pr$L0[j, d] + drop(X %*% pr$L1[j, d, ])
  Yi <- spec$Y; Yi[, cont] <- 0; storage.mode(Yi) <- "integer"
  tau_off <- c(0L, cumsum(spec$K - 1L))
  Lr <- function(d, e) if (D == 1) rep(1, N) else st$Lr[[d]][[e]]
  CL <- matrix(0, N, Tn)                                 # Cholesky of Sigma_i = diag(s) R diag(s)
  for (d in seq_len(D)) for (e in seq_len(d)) CL[, .lt_index(d, e)] <- st$s[, d] * Lr(d, e)
  ad <- !is.null(adapt); if (ad) adapt <- .adapt_compat(adapt, D)
  res <- mofia_core_cpp(
    Yi, NU, LAM, as.integer(as.vector(spec$loadings)), unlist(pr$tau, use.names = FALSE), as.integer(tau_off),
    st$mu, CL, nodes$T, nodes$logw, ad, if (ad) adapt$m else matrix(0, 0, D), if (ad) adapt$Lc else matrix(0, 0, Tn),
    if (link == "logit") 0L else 1L, scores, return_post, .mofia_threads(),
    YC, RV, as.integer(cont))
  if (res$bad > 0) return(-Inf)
  out <- if (per_person) res$ll else sum(res$ll)
  if (return_post) attr(out, "adapt") <- .adapt_compat(list(m = res$PM, Lc = .chol_cov_batch(res$PC, D)), D)
  if (!scores) return(out)

  pt <- spec$par; S <- matrix(0, N, nrow(pt), dimnames = list(NULL, pt$name))
  for (j in seq_len(P)) {
    for (d in seq_len(D)) if (spec$loadings[j, d]) {
      pge <- res$PGE[, j + P * (d - 1)]
      S[, pt$block == "lambda0" & pt$item == j & pt$factor %in% d] <- pge
      for (c in which(pt$block == "lambda1" & pt$item == j & pt$factor %in% d)) S[, c] <- X[, pt$mod[c]] * pge
    }
    for (c in which(pt$block == "beta" & pt$item == j)) S[, c] <- X[, pt$mod[c]] * res$PG[, j]
    if (cont[j]) {
      S[, pt$block == "nu0" & pt$item == j] <- res$PG[, j]
      S[, pt$block == "omega0" & pt$item == j] <- res$PH[, j]
      for (c in which(pt$block == "omega" & pt$item == j)) S[, c] <- X[, pt$mod[c]] * res$PH[, j]
      next
    }
    obs <- !is.na(spec$Y[, j]); yy <- ifelse(obs, spec$Y[, j], 0L)
    tcols <- which(pt$block == "tau" & pt$item == j)
    for (kk in seq_along(tcols))
      S[, tcols[kk]] <- ifelse(obs & yy == kk, res$DLO[, j], 0) + ifelse(obs & yy + 1 == kk, res$DHI[, j], 0)
  }
  for (c in which(pt$block == "gamma")) S[, c] <- X[, pt$mod[c]] * res$DMU[, pt$factor[c]]
  MZ <- function(d, e) res$MZ[, (d - 1) * D + e]           # E[A_d z_e]
  # d l / d s_d holding the correlation matrix fixed
  ds <- sapply(seq_len(D), function(d) Reduce(`+`, lapply(seq_len(d), function(e) Lr(d, e) * MZ(d, e)))); ds <- matrix(ds, N, D)
  W <- matrix(0, N, D + max(0, D * (D - 1) / 2))
  if (D >= 2) {
    pairs <- st$pairs; np <- nrow(pairs); s_ <- st$s; rho <- st$rho
    drho <- matrix(0, N, np)                               # d l / d rho_p holding s fixed
    for (p in seq_len(np)) {
      dL <- .chol_dir_batch(st$Lr, pairs[p, ], D, N)
      drho[, p] <- Reduce(`+`, lapply(2:D, function(d) s_[, d] * Reduce(`+`, lapply(seq_len(d), function(e) dL[[d]][[e]] * MZ(d, e)))))
    }
    # d l / d s_d holding the covariances fixed (rho_ab = cov_ab / (s_a s_b))
    ds_cov <- ds
    for (p in seq_len(np)) for (d in pairs[p, ]) ds_cov[, d] <- ds_cov[, d] - drho[, p] * rho[, p] / s_[, d]
    dc <- drho / (s_[, pairs[, 1], drop = FALSE] * s_[, pairs[, 2], drop = FALSE])
    pp <- pt$pair
    if (spec$cov_param == "fisherz") {
      du <- drho * (1 - rho^2)
      for (c in which(pt$block == "r0")) S[, c] <- du[, pp[c]]
      for (c in which(pt$block == "r1")) S[, c] <- X[, pt$mod[c]] * du[, pp[c]]
    } else {
      ds <- ds_cov
      for (c in which(pt$block == "r0")) S[, c] <- dc[, pp[c]]
      for (c in which(pt$block == "r1")) S[, c] <- X[, pt$mod[c]] * dc[, pp[c]]
    }
    W[, seq_len(D)] <- ds_cov * s_ / 2; W[, D + seq_len(np)] <- dc
    colnames(W) <- c(paste0("logvf", seq_len(D)), if (np == 1) "cov" else paste0("cov", seq_len(np)))
  } else { W[, 1] <- ds[, 1] * st$s[, 1] / 2; colnames(W) <- "logvf1" }
  for (c in which(pt$block == "delta")) { d <- pt$factor[c]; S[, c] <- X[, pt$mod[c]] * ds[, d] * st$s[, d] / 2 }
  for (c in which(pt$block == "delta0")) { d <- pt$factor[c]; S[, c] <- ds[, d] * st$s[, d] / 2 }
  attr(out, "scores") <- S
  attr(out, "phi_scores") <- W
  attr(out, "resvar_scores") <- res$PH
  out
}

# forward-mode derivative of the batch Cholesky factor of a correlation matrix in
# the direction of correlation (a, b): dR has ones at (a, b) and (b, a)
.chol_dir_batch <- function(L, ab, D, N) {
  dL <- lapply(seq_len(D), function(i) lapply(seq_len(i), function(j) numeric(N)))
  dR <- function(i, j) if ((i == ab[1] && j == ab[2]) || (i == ab[2] && j == ab[1])) 1 else 0
  for (j in seq_len(D)) {
    acc <- if (j > 1) Reduce(`+`, lapply(seq_len(j - 1), function(k) L[[j]][[k]] * dL[[j]][[k]])) else 0
    dL[[j]][[j]] <- (dR(j, j) - 2 * acc) / (2 * L[[j]][[j]])
    if (j < D) for (i in (j + 1):D) {
      acc2 <- if (j > 1) Reduce(`+`, lapply(seq_len(j - 1), function(k) dL[[i]][[k]] * L[[j]][[k]] + L[[i]][[k]] * dL[[j]][[k]])) else 0
      dL[[i]][[j]] <- (dR(i, j) - acc2 - L[[i]][[j]] * dL[[j]][[j]]) / L[[j]][[j]]
    }
  }
  dL
}
# batch Cholesky of covariance matrices given as row-wise lower triangles (N x T)
.chol_cov_batch <- function(PC, D) {
  N <- nrow(PC); Lc <- matrix(0, N, ncol(PC))
  for (j in seq_len(D)) {
    v <- PC[, .lt_index(j, j)] - (if (j > 1) rowSums(Lc[, .lt_index(j, seq_len(j - 1)), drop = FALSE]^2) else 0)
    Lc[, .lt_index(j, j)] <- sqrt(pmax(v, 1e-10))
    if (j < D) for (i in (j + 1):D) {
      v2 <- PC[, .lt_index(i, j)] - (if (j > 1) rowSums(Lc[, .lt_index(i, seq_len(j - 1)), drop = FALSE] * Lc[, .lt_index(j, seq_len(j - 1)), drop = FALSE]) else 0)
      Lc[, .lt_index(i, j)] <- v2 / Lc[, .lt_index(j, j)]
    }
  }
  Lc
}

# ---- curvature of the nonlinear parameter transformations -----------------------
# The person-level variances exp(.), Fisher-z correlations tanh(.) * sd * sd and
# continuous residual variances exp(.) are nonlinear in the parameters. The exact
# observed information is A = -H; the "reduced" information drops
#   C = sum_i sum_phi (d l_i / d phi_i) * d2 phi_i / d theta d theta'
# (zero for transformations linear in the parameters: loadings, linear covariance),
# i.e. it uses the chain rule through the person-level parameters only. Both are
# consistent: C has expectation zero at the true parameters.
.constraint_curvature <- function(theta, spec, W) {
  pt <- spec$par; X <- spec$X; N <- spec$N; np <- nrow(pt)
  C <- matrix(0, np, np, dimnames = list(pt$name, pt$name))
  xcol <- function(c) if (pt$block[c] == "delta0") rep(1, N) else X[, pt$mod[c]]
  for (d in seq_len(spec$D)) {            # vf_d = exp(delta0_d + x' delta_d)
    idx <- which(pt$block %in% c("delta", "delta0") & pt$factor == d)
    for (a in idx) for (b in idx) C[a, b] <- C[a, b] + sum(W[, d] * xcol(a) * xcol(b))
  }
  if (!is.null(R <- attr(W, "resvar")) && any(spec$cont)) {   # theta_j = exp(omega0_j + x' omega_j)
    xo <- function(c) if (pt$block[c] == "omega0") rep(1, N) else X[, pt$mod[c]]
    for (j in which(spec$cont)) {
      idx <- which(pt$block %in% c("omega0", "omega") & pt$item == j)
      for (a in idx) for (b in idx) C[a, b] <- C[a, b] + sum(R[, j] * xo(a) * xo(b))
    }
  }
  if (spec$D >= 2 && spec$cov_param == "fisherz") {   # cov_ab = tanh(u_ab) * exp((h_a + h_b) / 2)
    pr <- .unpack(theta, spec); pairs <- if (is.null(spec$pairs)) matrix(c(1L, 2L), 1) else spec$pairs
    pp <- if (is.null(pt$pair)) rep(1L, nrow(pt)) else pt$pair
    for (p in seq_len(nrow(pairs))) {
      a <- pairs[p, 1]; b <- pairs[p, 2]
      u <- pr$r0[p] + drop(X %*% pr$r1[p, ]); t <- tanh(u)
      eh <- exp((drop(X %*% (pr$Dl[a, ] + pr$Dl[b, ])) + pr$Dl0[a] + pr$Dl0[b]) / 2)
      cuu <- -2 * t * (1 - t^2) * eh; cuh <- (1 - t^2) * eh; chh <- t * eh
      Ga <- matrix(0, N, np); Gb <- matrix(0, N, np)
      Ga[, pt$block == "r0" & pp %in% p] <- 1
      for (c in which(pt$block == "r1" & pp %in% p)) Ga[, c] <- X[, pt$mod[c]]
      for (c in which(pt$block == "delta" & pt$factor %in% c(a, b))) Gb[, c] <- X[, pt$mod[c]] / 2
      Gb[, pt$block == "delta0" & pt$factor %in% c(a, b)] <- 1 / 2
      wc <- W[, spec$D + p]
      C <- C + crossprod(Ga * (wc * cuu), Ga) + crossprod(Ga * (wc * cuh), Gb) +
               crossprod(Gb * (wc * cuh), Ga) + crossprod(Gb * (wc * chh), Gb)
    }
  }
  C
}

# adaptive-quadrature state: m (N x D posterior means of z) and Lc (N x D(D+1)/2,
# row-wise lower triangle of the posterior Cholesky factor); for D <= 2 also the
# fields L11, L21, L22 used by the R reference engine.
.adapt_compat <- function(a, D) {
  if (is.null(a$Lc)) { a$Lc <- if (D == 1) cbind(a$L11) else cbind(a$L11, a$L21, a$L22) }
  if (D <= 2 && is.null(a$L11)) { a$L11 <- a$Lc[, 1]; if (D == 2) { a$L21 <- a$Lc[, 2]; a$L22 <- a$Lc[, 3] } }
  a
}

# ---- adaptive quadrature -------------------------------------------------------
.mofia_adapt <- function(theta, spec, link = "logit", nodes = .mofia_nodes(.default_nq(spec$D), spec$D),
                        adapt = NULL, maxit = 50, tol = 1e-8) {
  if (is.null(adapt)) {   # pilot pass with a dense non-adaptive rule
    pilot <- .mofia_nodes(c(101, 41, 13, 7, 5, 4)[min(spec$D, 6)], spec$D)
    adapt <- attr(.mofia_ll(theta, spec, link, pilot, return_post = TRUE), "adapt")
  }
  ll_old <- -Inf
  for (it in seq_len(maxit)) {
    ll <- .mofia_ll(theta, spec, link, nodes, adapt = adapt, return_post = TRUE)
    if (!is.finite(ll)) break
    adapt <- attr(ll, "adapt")
    if (abs(ll - ll_old) < tol) break
    ll_old <- as.numeric(ll)
  }
  adapt
}


#' Log-likelihood at given parameter values
#'
#' Evaluates the marginal log-likelihood of a model at user-supplied parameter
#' values, e.g. to compare with estimates obtained elsewhere.
#'
#' @param object A `mofia_model` or a fitted `mofia` object.
#' @param theta Named parameter vector (names as in `object$par$name`). For a
#'   fitted object, defaults to its estimates.
#' @param link `"logit"` or `"probit"` (taken from a fitted object by default).
#' @param nq Quadrature nodes per factor dimension.
#' @return The log-likelihood (a number).
#' @examples
#' data(parabelief)
#' m <- mofia_model(parabelief, items = c("SpiritsGhosts", "GnomesElves", "Aliens"),
#'                  moderators = "Male")
#' th <- setNames(c(rep(c(1, 2, 3, 4), 3), 1.5, 1.5, 1.5, 0, 0), m$par$name)
#' mofia_loglik(m, th)
#' @export
mofia_loglik <- function(object, theta = NULL, link = NULL, nq = NULL) {
  model <- if (inherits(object, "mofia")) object$model else object
  if (is.null(theta)) { if (!inherits(object, "mofia")) stop("Give 'theta'."); theta <- object$coefficients }
  if (is.null(link)) link <- if (inherits(object, "mofia")) object$link else "logit"
  if (is.null(nq)) nq <- if (inherits(object, "mofia")) object$nq else .default_nq(model$D)
  if (is.null(names(theta))) stop("'theta' must be named.")
  miss <- setdiff(model$par$name, names(theta)); if (length(miss)) stop("Missing parameter(s): ", paste(miss, collapse = ", "))
  theta <- theta[model$par$name]
  nodes <- .mofia_nodes(nq, model$D)
  as.numeric(.mofia_ll(theta, model, link, nodes, adapt = .mofia_adapt(theta, model, link, nodes)))
}
