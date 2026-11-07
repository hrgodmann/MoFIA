# Estimation ----------------------------------------------------------------------

#' Fit a moderated factor model
#'
#' Estimates a moderated (nonlinear) factor model by marginal maximum
#' likelihood. The factors are integrated out by adaptive product
#' Gauss-Hermite quadrature; the log-likelihood and its analytic gradient are
#' computed in C++.
#'
#' @section Standard errors:
#' Three kinds are computed: `"ML"` (inverse observed information), `"MLR"`
#' (robust sandwich, the default in [summary.mofia()]) and `"MLF"` (outer
#' product of the per-person scores, valid only for a correctly specified
#' model). ML and MLR use the observed information, by default the full
#' Hessian of the log-likelihood (`information = "exact"`). This is the
#' derivative of the estimating equations, which the sandwich needs to be
#' robust to misspecification.
#'
#' `information = "reduced"` is available for comparison with software that
#' differentiates through the person-level parameters only: it omits the term
#' with the second derivatives of the nonlinear parameter transformations
#' (log-linear factor and residual variances, Fisher-z correlations). That term
#' has expectation zero when the model is correctly specified, so reduced ML
#' standard errors are consistent then. Under misspecification it does not
#' vanish, and reduced MLR standard errors and scaling factors are not robust;
#' use them only as a correct-specification check.
#'
#' @section Quadrature:
#' The factors are integrated out per person with adaptive Gauss-Hermite
#' quadrature: the nodes are centred and scaled at each person's posterior. With
#' only continuous items the posterior is normal and known in closed form, so
#' the adaptive rule is exact. Otherwise, after fitting, the log-likelihood at
#' the estimates is recomputed with a finer rule; if it changes by more than
#' `accuracy`, the model is refitted with the finer rule (up to three times, with
#' a warning if the check still fails). Very discriminating categorical items
#' can need many more nodes than the default. The nodes used and the check are
#' stored in `$nq` and `$quadrature`.
#'
#' @section Convergence:
#' A fit counts as converged when the information matrix at the solution is
#' positive definite (the solution is a maximum, not a saddle point, and the
#' model is identified by the data) and the Newton decrement, the
#' log-likelihood gain still available from one more Newton step
#' (\eqn{g'A^{-1}g/2}), is below `1e-6`. If the information matrix is not
#' positive definite, the fit is marked as not converged, `status` is
#' `"saddle"` or `"singular"`, and standard errors and tests are not reported.
#' With `se = FALSE` only the gradient is checked.
#'
#' @param model A `mofia_model` from [mofia_model()].
#' @param link `"logit"` (default) or `"probit"` for categorical items.
#' @param nq Quadrature nodes per factor dimension (default 31, 15, 9, 7, 5, 4
#'   for 1 to 6 factors; 7, 7, 5, 5, 4, 3 for models with only continuous items,
#'   whose adaptive rule is exact). See *Quadrature*.
#' @param accuracy Tolerance of the quadrature check (log-likelihood units), or
#'   `NULL` to skip the check. See *Quadrature*.
#' @param start Optional named vector of starting values.
#' @param starts Number of starts. The first uses the default starting values;
#'   the others perturb them randomly (set a seed for reproducibility). The fit
#'   with the highest log-likelihood is returned; all log-likelihoods are kept
#'   in `$starts`. Useful when a model may have several local maxima.
#' @param se Compute standard errors? (Requires one numerical differentiation of
#'   the analytic gradient.)
#' @param control List of control options passed to [stats::nlminb()].
#' @param verbose Print progress?
#' @param threads Number of threads for the likelihood (OpenMP, where
#'   available). Default: `getOption("mofia.threads", 1)`.
#'
#' @return An object of class `mofia` with components including
#'   `coefficients`, `loglik`, `npar`, `nobs`, `converged`, `status` (`"ok"`,
#'   `"saddle"` or `"singular"`), `information_eigen` (smallest and largest
#'   eigenvalue of the information matrix), `quadrature` (log-likelihood at the
#'   estimates for the node counts compared), `newton_decrement`, `vcov` (a list
#'   of covariance matrices), `scaling` (robust scaling correction factors),
#'   `model`, `link` and `nq`.
#' @seealso [mofia_model()], [summary.mofia()], [anova.mofia()],
#'   [mofia_wald()], [mofia_simulate()]
#' @references
#' Bauer, D. J. (2017). A more general model for testing measurement
#' invariance and differential item functioning. *Psychological Methods,
#' 22*(3), 507--526. \doi{10.1037/met0000077}
#'
#' Bauer, D. J., & Hussong, A. M. (2009). Psychometric approaches for
#' developing commensurate measures across independent studies: Traditional and
#' new models. *Psychological Methods, 14*(2), 101--125.
#' \doi{10.1037/a0015583}
#' @examples
#' data(parabelief)
#' ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession",
#'         "Sorcery", "HauntedPlaces")
#' m <- mofia_model(parabelief, items = ub, moderators = c("Age", "Male"),
#'                  threshold_mod = "DevilPossession")
#' fit <- mofia(m)
#' fit
#' summary(fit)
#' @export
mofia <- function(model, link = c("logit", "probit"), nq = NULL, start = NULL, se = TRUE,
                  starts = 1, control = list(), verbose = FALSE, threads = NULL, accuracy = 0.005) {
  if (!inherits(model, "mofia_model")) stop("'model' must be created with mofia_model().")
  link <- match.arg(link)
  cont_only <- all(model$cont)
  if (cont_only) link <- "logit"   # irrelevant for continuous items
  if (!is.null(threads)) { old <- options(mofia.threads = threads); on.exit(options(old), add = TRUE) }
  if (is.null(nq)) nq <- .default_nq(model$D, continuous = cont_only)
  starts <- max(1L, as.integer(starts))
  call <- match.call()
  nodes <- .mofia_nodes(nq, model$D)
  st <- .mofia_start(model, link)
  if (!is.null(start)) {
    if (is.null(names(start))) stop("'start' must be a named vector.")
    keep <- intersect(names(start), names(st)); st[keep] <- start[keep]
  }
  ctrl <- utils::modifyList(list(eval.max = 3000, iter.max = 2000, rel.tol = 1e-10), control)
  t0 <- proc.time()[[3]]
  runs <- vector("list", starts)
  runs[[1]] <- .mofia_optimize(model, link, nodes, st, ctrl, verbose)
  for (r in seq_len(starts)[-1]) {
    if (verbose) message(sprintf("random start %d", r))
    runs[[r]] <- tryCatch(.mofia_optimize(model, link, nodes, .perturb_start(st, model), ctrl, verbose),
                          error = function(e) NULL)
  }
  lls <- vapply(runs, function(r) if (is.null(r)) NA_real_ else r$loglik, 1)
  best <- which.max(lls); run <- runs[[best]]
  if (starts > 1 && lls[best] - lls[1] > 1e-3)
    message(sprintf("Random start %d reached a higher log-likelihood than the default start (+%.4f).", best, lls[best] - lls[1]))
  theta <- run$theta; adapt <- run$adapt; opt <- run$opt; outer <- run$outer

  # quadrature accuracy: recompute the log-likelihood at the estimates with a finer rule;
  # if it changes by more than 'accuracy', refit with the finer rule (up to 3 times).
  # Models with only continuous items use exact adaptation and need no check.
  quad <- NULL
  if (!is.null(accuracy) && !cont_only) {
    quad <- data.frame(nq = integer(0), logLik = numeric(0))
    for (ref in 0:3) {
      adapt <- .mofia_adapt(theta, model, link, nodes, adapt)
      ll0 <- as.numeric(.mofia_ll(theta, model, link, nodes, adapt = adapt))
      nqf <- .nq_next(nq, model$D)
      if (nqf > .nq_cap(model$D)) {
        if (ref > 0) warning(sprintf("The quadrature check reached its limit (%d nodes per factor); the log-likelihood may still be inaccurate.", nq), call. = FALSE)
        break
      }
      nodes_f <- .mofia_nodes(nqf, model$D)
      llf <- as.numeric(.mofia_ll(theta, model, link, nodes_f, adapt = .mofia_adapt(theta, model, link, nodes_f, adapt)))
      quad <- rbind(quad[quad$nq != nq, ], data.frame(nq = c(nq, nqf), logLik = c(ll0, llf)))
      if (is.finite(llf) && is.finite(ll0) && abs(llf - ll0) <= accuracy) break
      if (ref == 3) {
        warning(sprintf(paste0("The log-likelihood still changes by %.3g between %d and %d quadrature nodes per factor; ",
                               "results may be inaccurate (increase nq)."), llf - ll0, nq, nqf), call. = FALSE)
        break
      }
      if (verbose) message(sprintf("quadrature check: log-likelihood changes by %.3g with %d nodes; refitting with %d", llf - ll0, nqf, nqf))
      nq <- nqf; nodes <- nodes_f
      run <- .mofia_optimize(model, link, nodes, theta, ctrl, verbose)
      theta <- run$theta; adapt <- run$adapt; opt <- run$opt; outer <- outer + run$outer
    }
    rownames(quad) <- NULL
  }
  adapt <- .mofia_adapt(theta, model, link, nodes, adapt)
  if (!isTRUE(attr(adapt, "converged")))
    warning("The adaptive quadrature did not settle at the solution; the log-likelihood may be inaccurate ",
            "(very concentrated posteriors, e.g. near-deterministic items).", call. = FALSE)
  ll <- .mofia_ll(theta, model, link, nodes, scores = TRUE, adapt = adapt)
  if (!is.finite(ll)) stop("Estimation failed: the log-likelihood is not finite at the solution.")
  g <- colSums(attr(ll, "scores"))
  fit <- list(coefficients = theta, loglik = as.numeric(ll), npar = length(theta), nobs = model$N,
              gradient = g, max_abs_grad = max(abs(g)), converged = max(abs(g)) < 1e-3,
              newton_decrement = NA_real_, optimizer = list(convergence = opt$convergence, message = opt$message,
              iterations = opt$iterations, outer_iterations = outer),
              link = link, nq = nq, quadrature = quad, adapt = adapt, time = unname(proc.time()[[3]] - t0),
              starts = if (starts > 1) data.frame(start = seq_len(starts), logLik = lls) else NULL,
              model = model, call = call)
  class(fit) <- "mofia"
  if (se) fit <- .mofia_se(fit)
  if (!isTRUE(fit$converged) && identical(fit$status %||% "ok", "ok"))
    warning("The model may not have converged (max |gradient| = ", signif(fit$max_abs_grad, 3),
            if (is.finite(fit$newton_decrement)) paste0(", Newton decrement = ", signif(fit$newton_decrement, 3)), ").",
            call. = FALSE)
  fit
}

# one optimization run from a given start (outer loop re-adapts the quadrature)
.mofia_optimize <- function(model, link, nodes, start, ctrl, verbose) {
  raw <- .to_raw(start, model); adapt <- NULL; ll_prev <- -Inf
  for (outer in seq_len(20)) {
    adapt <- .mofia_adapt(.to_nat(raw, model), model, link, nodes, adapt)
    # objective and gradient share one likelihood evaluation
    cache <- new.env(); cache$r <- NULL
    eval_at <- function(r) {
      if (is.null(cache$r) || !identical(r, cache$r)) {
        ll <- .mofia_ll(.to_nat(r, model), model, link, nodes, scores = TRUE, adapt = adapt)
        cache$r <- r; cache$v <- as.numeric(ll)
        cache$g <- if (is.finite(ll)) -.grad_raw(r, colSums(attr(ll, "scores")), model) else rep(0, length(r))
      }
      cache
    }
    fn <- function(r) { v <- -eval_at(r)$v; if (!is.finite(v)) 1e100 else v }
    gr <- function(r) eval_at(r)$g
    opt <- stats::nlminb(raw, fn, gr, control = ctrl)
    raw <- opt$par
    if (verbose) message(sprintf("outer %d: logLik %.6f (%d iterations)", outer, -opt$objective, opt$iterations))
    if (abs(-opt$objective - ll_prev) < 1e-8) break
    ll_prev <- -opt$objective
  }
  theta <- .to_nat(raw, model); names(theta) <- model$par$name
  list(theta = theta, loglik = -opt$objective, adapt = adapt, opt = opt, outer = outer)
}

# random perturbation of the default starting values
.perturb_start <- function(theta, model) {
  pt <- model$par; th <- theta
  for (j in seq_len(model$P)) {
    s <- which(pt$block == "tau" & pt$item %in% j); if (!length(s)) next
    r <- .to_raw(th, model)[s]
    r <- r + c(stats::rnorm(1, 0, 0.5), stats::rnorm(length(s) - 1, 0, 0.25))
    th[s] <- cumsum(c(r[1], exp(r[-1])))
  }
  b <- pt$block
  th[b == "lambda0"] <- th[b == "lambda0"] * stats::runif(sum(b == "lambda0"), 0.5, 1.5)
  mod <- b %in% c("lambda1", "beta", "gamma", "delta", "omega", "r1")
  th[mod] <- stats::rnorm(sum(mod), 0, 0.1)
  th[b == "delta0"] <- stats::rnorm(sum(b == "delta0"), 0, 0.2)
  th[b == "r0"] <- stats::rnorm(sum(b == "r0"), 0.3, 0.2)
  th[b %in% c("nu0", "omega0")] <- th[b %in% c("nu0", "omega0")] + stats::rnorm(sum(b %in% c("nu0", "omega0")), 0, 0.1)
  th
}

# ---- standard errors -----------------------------------------------------------
# Relative eigenvalue threshold below which the information matrix counts as
# singular (min eigenvalue / max eigenvalue); negative beyond it: saddle point.
.pd_tol <- 1e-8

.mofia_se <- function(fit, check = TRUE) {
  model <- fit$model; theta <- fit$coefficients; link <- fit$link; adapt <- fit$adapt
  nodes <- .mofia_nodes(fit$nq, model$D)
  ll <- .mofia_ll(theta, model, link, nodes, scores = TRUE, adapt = adapt)
  S <- attr(ll, "scores")
  # observed information: central differences of the analytic gradient
  gfun <- function(th) colSums(attr(.mofia_ll(th, model, link, nodes, scores = TRUE, adapt = adapt), "scores"))
  np <- length(theta); h <- 1e-5 * pmax(1, abs(theta)); H <- matrix(0, np, np)
  for (k in seq_len(np)) {
    e <- replace(numeric(np), k, h[k])
    H[, k] <- (gfun(theta + e) - gfun(theta - e)) / (2 * h[k])
  }
  H <- (H + t(H)) / 2
  dn <- list(names(theta), names(theta))
  B <- crossprod(S); dimnames(B) <- dn
  I_exact <- -H
  W <- attr(ll, "phi_scores"); attr(W, "resvar") <- attr(ll, "resvar_scores")
  I_red <- -H + .constraint_curvature(theta, model, W)
  g <- colSums(S)

  # curvature: a maximum needs a positive-definite information matrix
  ev <- if (all(is.finite(I_exact))) eigen(I_exact, symmetric = TRUE, only.values = TRUE)$values else NA_real_
  fit$information_eigen <- c(min = min(ev), max = max(ev))
  pd <- all(is.finite(ev)) && min(ev) > .pd_tol * max(abs(ev))
  fit$status <- if (pd) "ok" else if (all(is.finite(ev)) && min(ev) < -.pd_tol * max(abs(ev))) "saddle" else "singular"
  if (check && !pd) {
    fit$converged <- FALSE; fit$vcov <- NULL; fit$scaling <- NULL; fit$newton_decrement <- NA_real_
    fit$se_note <- if (fit$status == "saddle")
      sprintf("the solution is a saddle point, not a maximum (smallest information eigenvalue %.3g)", min(ev))
    else sprintf("the information matrix is singular (smallest eigenvalue %.3g, largest %.3g): the model is not identified, or not identified by these data",
                 min(ev), max(ev))
    warning("Standard errors and tests are not available: ", fit$se_note, ". ",
            if (fit$status == "saddle") "Try other starting values (starts > 1)." else "Check the anchor items.",
            call. = FALSE)
    return(fit)
  }

  inv_pd <- function(I) {
    e <- eigen(I, symmetric = TRUE)
    if (!all(is.finite(e$values)) || min(e$values) <= .pd_tol * max(abs(e$values))) return(NULL)
    V <- e$vectors %*% (t(e$vectors) / e$values); V <- (V + t(V)) / 2; dimnames(V) <- dn; V
  }
  sandwich <- function(Ai) if (is.null(Ai)) NULL else { V <- Ai %*% B %*% Ai; (V + t(V)) / 2 }
  A_exact <- inv_pd(I_exact); A_red <- inv_pd(I_red)
  fit$vcov <- list(exact = list(ML = A_exact, MLR = sandwich(A_exact)),
                   reduced = if (is.null(A_red)) NULL else list(ML = A_red, MLR = sandwich(A_red)),
                   MLF = inv_pd(B))
  fit$scaling <- c(exact = if (is.null(A_exact)) NA_real_ else sum(diag(A_exact %*% B)) / np,
                   reduced = if (is.null(A_red)) NA_real_ else sum(diag(A_red %*% B)) / np)
  fit$information <- "exact"
  fit$newton_decrement <- if (is.null(A_exact)) NA_real_ else 0.5 * sum(g * (A_exact %*% g))
  fit$converged <- pd && is.finite(fit$newton_decrement) && fit$newton_decrement < 1e-6
  fit
}

# ---- reparameterisation (ordered thresholds) -------------------------------------
.to_raw <- function(theta, spec) {
  raw <- theta; pt <- spec$par
  for (j in seq_len(spec$P)) {
    s <- which(pt$block == "tau" & pt$item %in% j); if (!length(s)) next
    t <- theta[s]; raw[s] <- c(t[1], log(diff(t)))
  }
  raw
}
.to_nat <- function(raw, spec) {
  theta <- raw; pt <- spec$par
  for (j in seq_len(spec$P)) {
    s <- which(pt$block == "tau" & pt$item %in% j); if (!length(s)) next
    r <- raw[s]; theta[s] <- cumsum(c(r[1], exp(r[-1])))
  }
  theta
}
.grad_raw <- function(raw, g_nat, spec) {
  g <- g_nat; pt <- spec$par
  for (j in seq_len(spec$P)) {
    s <- which(pt$block == "tau" & pt$item %in% j); if (!length(s)) next
    ts <- rev(cumsum(rev(g_nat[s])))
    g[s] <- c(ts[1], exp(raw[s][-1]) * ts[-1])
  }
  g
}

.mofia_start <- function(spec, link = "logit") {
  pt <- spec$par; theta <- numeric(nrow(pt)); names(theta) <- pt$name
  q <- if (link == "logit") stats::qlogis else stats::qnorm
  for (j in seq_len(spec$P)) {
    if (isTRUE(spec$cont[j])) next
    y <- spec$Y[, j]; y <- y[!is.na(y)]
    cum <- pmin(pmax(sapply(seq_len(spec$K[j] - 1), function(k) mean(y >= k)), 0.001), 0.999)
    theta[pt$block == "tau" & pt$item %in% j] <- -q(cum) * 1.2
  }
  theta[pt$block == "lambda0"] <- if (link == "logit") 1 else 0.6
  for (j in which(spec$cont)) {
    y <- spec$Y[, j]; v <- stats::var(y, na.rm = TRUE)
    theta[pt$block == "nu0" & pt$item %in% j] <- mean(y, na.rm = TRUE)
    theta[pt$block == "omega0" & pt$item %in% j] <- log(v / 2)
    theta[pt$block == "lambda0" & pt$item %in% j] <- sqrt(v / 2)
  }
  if (spec$D >= 2) theta[pt$block == "r0"] <- 0.3
  theta
}
