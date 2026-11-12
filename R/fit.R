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
#' the adaptive rule is exact. Where the posterior is very concentrated (a
#' residual variance near zero, as in a Heywood case), it is computed in
#' square-root form, which keeps the rule exact. At the estimates, the log-likelihood
#' is also compared with its closed form (the normal density of the items); if
#' they differ by more than `1e-6`, which only rounding at residual variances
#' of about `1e-18` or below can cause, estimation fails (see *Convergence*).
#' With categorical items, after fitting, the log-likelihood at
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
#' Positive definiteness is judged on the information matrix scaled to unit
#' diagonal (smallest eigenvalue relative to the largest, tolerance `1e-8`), so
#' that it does not depend on the units of the parameters (e.g. a moderator in
#' years and its square); the covariance matrices are transformed back to the
#' parameters' units. Scaling removes differences in units, not near linear
#' dependence among the parameters. With `se = FALSE` only the gradient is
#' checked (maximum absolute gradient below `1e-3`, which does depend on units),
#' and `status` is `NULL`.
#'
#' The optimizer alternates between adapting the quadrature at the current
#' estimates and maximizing the log-likelihood with that adaptation held fixed
#' (outer iterations, at most 20). Away from the point where it was adapted a
#' fixed rule is no longer exact. So that the optimizer does not follow its
#' error, an outer iteration that lowers the re-adapted log-likelihood has its
#' step halved, up to 10 times, until it does not. If no halving helps, the
#' iteration is undone and retried with a tenth of the iterations of
#' [stats::nlminb()] (down to 2), which keeps it nearer to where the rule was
#' adapted; if that does not help either, the optimization stops there. A fall counts if it
#' exceeds `1e-6` with only continuous items (whose adaptive rule is exact) and
#' `accuracy` (default `0.005`) otherwise, the tolerance of the quadrature.
#' With only continuous items, once the loop has had to act, and while the
#' posterior of the factors is very concentrated for some person (a posterior
#' standard deviation below `1e-4`, e.g. a residual variance near zero), the
#' adaptation is instead recomputed at every evaluation, so that the optimizer
#' maximizes the exact log-likelihood.
#' These actions are recorded in `$optimizer$monotone`. A fit that stopped this
#' way is usually not at a maximum, which `status` and `converged` then show.
#'
#' Estimation fails, and `mofia()` stops with an error of class
#' `"mofia_estimation_error"` (message "Estimation failed: ..."), when the
#' log-likelihood is not finite at the solution or, with `se = TRUE`, within a
#' finite-difference step of it, so that the information matrix cannot be
#' computed, and, with only continuous items, when the log-likelihood from
#' quadrature at the solution differs from its closed form (see *Quadrature*).
#' The estimates are then at the edge of the region where the
#' likelihood can be evaluated, typically in a model that is not identified:
#' for example, a residual variance at the limit of double precision (huge
#' after a diverged optimization, or near zero), or a singular factor
#' correlation matrix. The message
#' names the parameters concerned and the optimizer's final state. No fit is
#' returned: such estimates are not a regular maximum, and after a diverged
#' optimization their log-likelihood does not describe a solution.
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
#'   available). Default: `getOption("mofia.threads", 1)`. If the package was
#'   built without OpenMP support (for example with a compiler that lacks it),
#'   more than one thread cannot be used: a message says so and the likelihood
#'   is computed serially (`$threads` of the fit records requested and used
#'   threads).
#'
#' @return An object of class `mofia` with components including
#'   `coefficients`, `loglik`, `npar`, `nobs`, `converged`, `status` (`"ok"`,
#'   `"saddle"` or `"singular"`; `NULL` with `se = FALSE`), `se_note` (why
#'   standard errors are not available, if they are not), `information_eigen`
#'   (smallest and largest eigenvalue of the information matrix, in the parameters' units),
#'   `information_eigen_scaled` (the same for the information matrix scaled to
#'   unit diagonal, which `status` is based on), `threads`, `quadrature` (log-likelihood at the
#'   estimates for the node counts compared), `newton_decrement`, `vcov` (a list
#'   of covariance matrices), `scaling` (robust scaling correction factors),
#'   `optimizer` (the final [stats::nlminb()] run's `convergence`, `message`
#'   and `iterations`, the number of `outer_iterations`, and, only if the
#'   monotone outer loop acted, `monotone`: one row per such outer iteration
#'   with the node count `nq`, its `iter_max` for [stats::nlminb()], the
#'   re-adapted log-likelihood at its start and end, the number of step
#'   `halvings` (`NA` if it was undone) and the resulting `loglik`; and, only
#'   if there were any, `exact_outer`, the number of outer iterations in which
#'   the adaptation was recomputed at every evaluation), `model`, `link` and
#'   `nq`.
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
  thr <- .thread_use()
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
  tol <- .mono_tol(cont_only, accuracy)
  runs[[1]] <- .mofia_optimize(model, link, nodes, st, ctrl, verbose, tol)
  for (r in seq_len(starts)[-1]) {
    if (verbose) message(sprintf("random start %d", r))
    runs[[r]] <- tryCatch(.mofia_optimize(model, link, nodes, .perturb_start(st, model), ctrl, verbose, tol),
                          error = function(e) NULL)
  }
  lls <- vapply(runs, function(r) if (is.null(r)) NA_real_ else r$loglik, 1)
  best <- which.max(lls); run <- runs[[best]]
  if (starts > 1 && lls[best] - lls[1] > 1e-3)
    message(sprintf("Random start %d reached a higher log-likelihood than the default start (+%.4f).", best, lls[best] - lls[1]))
  theta <- run$theta; adapt <- run$adapt; opt <- run$opt; outer <- run$outer; mono <- run$monotone; n_exact <- run$exact_outer

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
      run <- .mofia_optimize(model, link, nodes, theta, ctrl, verbose, tol)
      theta <- run$theta; adapt <- run$adapt; opt <- run$opt; outer <- outer + run$outer; mono <- rbind(mono, run$monotone)
      n_exact <- n_exact + run$exact_outer
    }
    rownames(quad) <- NULL
  }
  adapt <- .mofia_adapt(theta, model, link, nodes, adapt)
  if (!isTRUE(attr(adapt, "converged")))
    warning("The adaptive quadrature did not settle at the solution; the log-likelihood may be inaccurate ",
            "(very concentrated posteriors, e.g. near-deterministic items).", call. = FALSE)
  ll <- .mofia_ll(theta, model, link, nodes, scores = TRUE, adapt = adapt)
  if (!is.finite(ll)) .estimation_error("Estimation failed: the log-likelihood is not finite at the solution.")
  if (cont_only) .check_closed_form(as.numeric(ll), theta, model)
  g <- colSums(attr(ll, "scores"))
  fit <- list(coefficients = theta, loglik = as.numeric(ll), npar = length(theta), nobs = model$N,
              gradient = g, max_abs_grad = max(abs(g)), converged = max(abs(g)) < 1e-3,
              newton_decrement = NA_real_, optimizer = c(list(convergence = opt$convergence, message = opt$message,
              iterations = opt$iterations, outer_iterations = outer), if (!is.null(mono)) list(monotone = mono),
              if (n_exact > 0) list(exact_outer = n_exact)),
              link = link, nq = nq, quadrature = quad, adapt = adapt, time = unname(proc.time()[[3]] - t0), threads = thr,
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

# threads requested (option mofia.threads / argument threads) and threads the build
# can use: without OpenMP support the likelihood is computed serially, which is said
# rather than silently ignoring the request
.thread_use <- function() {
  req <- .mofia_threads(); omp <- mofia_openmp_threads()
  if (req > 1 && omp < 1)
    message(sprintf("This build of MoFIA has no OpenMP support: %d threads were requested, the likelihood is computed serially.", req))
  c(requested = req, openmp = omp, used = if (omp < 1) 1L else req)
}

# Tolerance of the monotone outer loop (log-likelihood units), see .mofia_optimize(): with
# only continuous items the re-adapted log-likelihood is exact up to rounding; otherwise it
# carries quadrature error, and a fall within the quadrature's own tolerance ('accuracy' of
# mofia(), by default 0.005) is not evidence that the optimizer followed that error.
.mono_tol <- function(cont_only, accuracy) if (cont_only) 1e-6 else if (is.null(accuracy)) 0.005 else accuracy

# One optimization run from a given start. Each outer iteration adapts the quadrature at
# the current point and runs nlminb with that adaptation held fixed. Away from the point
# where it was adapted the fixed rule is no longer exact, and the optimizer can follow its
# error instead of the likelihood (e.g. a residual variance collapsing towards zero, which
# a fixed rule rewards). The outer loop is therefore monotone in the re-adapted
# log-likelihood: if it is lower at the end of an outer iteration than at its start (by more
# than 'tol', see .mono_tol), the step of that iteration is halved, up to 10 times, until it is not.
# If no halving is acceptable, the iteration is undone and retried from its start with a tenth
# of nlminb's iteration budget (staying nearer to where the rule was adapted), down to 2
# iterations; after that the loop stops. Each such action is a row of $monotone (NULL if
# there was none).
#
# With only continuous items, the fixed rule is fragile once the loop has had to act, and when
# some person's posterior is concentrated (a residual variance near zero): a small step then
# moves the posterior by many of its standard deviations, and the fixed rule rewards a
# collapsing residual variance. In those outer iterations the adaptation, the exact posterior,
# is recomputed at every evaluation, so that nlminb maximizes the exact log-likelihood;
# $exact_outer counts them. Without either, the iterations, the stopping rule and the result
# are those of the plain loop.
.mofia_optimize <- function(model, link, nodes, start, ctrl, verbose, tol) {
  raw <- .to_raw(start, model); ll_prev <- -Inf; mono <- NULL; opt_ok <- NULL; ctrl_run <- ctrl; n_exact <- 0L
  ad_next <- .mofia_adapt(.to_nat(raw, model), model, link, nodes, NULL)
  ll_start <- as.numeric(.mofia_ll(.to_nat(raw, model), model, link, nodes, adapt = ad_next))
  for (outer in seq_len(20)) {
    adapt <- ad_next
    exact <- all(model$cont) && (!is.null(mono) || .concentrated(adapt, model$D)); n_exact <- n_exact + exact
    # objective and gradient share one likelihood evaluation
    cache <- new.env(); cache$r <- NULL
    eval_at <- function(r) {
      if (is.null(cache$r) || !identical(r, cache$r)) {
        ll <- .mofia_ll(.to_nat(r, model), model, link, nodes, scores = TRUE,
                        adapt = if (exact) .mofia_adapt(.to_nat(r, model), model, link, nodes, adapt) else adapt)
        cache$r <- r; cache$v <- as.numeric(ll)
        cache$g <- if (is.finite(ll)) -.grad_raw(r, colSums(attr(ll, "scores")), model) else rep(0, length(r))
      }
      cache
    }
    fn <- function(r) { v <- -eval_at(r)$v; if (!is.finite(v)) 1e100 else v }
    gr <- function(r) eval_at(r)$g
    opt <- stats::nlminb(raw, fn, gr, control = ctrl_run)
    # re-adapted at the end point (also the adaptation the next iteration starts from)
    ad_next <- .mofia_adapt(.to_nat(opt$par, model), model, link, nodes, adapt)
    ll_end <- as.numeric(.mofia_ll(.to_nat(opt$par, model), model, link, nodes, adapt = ad_next))
    if (!isTRUE(ll_end >= ll_start - tol)) {
      step <- opt$par - raw; h <- 0L; ok <- FALSE
      while (!ok && h < 10L) {
        h <- h + 1L; r <- raw + step / 2^h
        a <- .mofia_adapt(.to_nat(r, model), model, link, nodes, adapt)
        l <- as.numeric(.mofia_ll(.to_nat(r, model), model, link, nodes, adapt = a))
        ok <- isTRUE(l >= ll_start - tol)
      }
      mono <- rbind(mono, data.frame(nq = nodes$nq, outer = outer, iter_max = ctrl_run$iter.max, loglik_start = ll_start,
                                     loglik_end = ll_end, halvings = if (ok) h else NA_integer_, loglik = if (ok) l else ll_start))
      if (verbose) message(sprintf("outer %d: re-adapted logLik fell from %.6f to %.6f; %s", outer, ll_start, ll_end,
                                   if (ok) sprintf("step halved %d time(s)", h) else "iteration undone"))
      if (!ok) {                                       # undone: keep raw and its adaptation
        if (!is.null(opt_ok)) opt <- opt_ok
        ad_next <- adapt
        if (ctrl_run$iter.max <= 2) break
        ctrl_run$iter.max <- max(2, ctrl_run$iter.max %/% 10); ll_prev <- -Inf; next
      }
      opt$par <- r; opt$objective <- -l; ad_next <- a; ll_end <- l
    }
    raw <- opt$par; ll_start <- ll_end; opt_ok <- opt
    if (verbose) message(sprintf("outer %d: logLik %.6f (%d iterations)", outer, -opt$objective, opt$iterations))
    if (abs(-opt$objective - ll_prev) < 1e-8) break
    ll_prev <- -opt$objective
  }
  theta <- .to_nat(raw, model); names(theta) <- model$par$name
  list(theta = theta, loglik = if (is.null(mono) && !n_exact) -opt$objective else ll_start, adapt = adapt, opt = opt,
       outer = outer, monotone = mono, exact_outer = n_exact)
}

# Is some person's posterior concentrated, with a standard deviation (Cholesky pivot) of at
# most 1 / sqrt(.cont_qmax), where .cont_posterior() switches to the square-root form?
.concentrated <- function(adapt, D) {
  if (is.null(adapt$Lc)) return(FALSE)
  piv <- adapt$Lc[, vapply(seq_len(D), function(d) .lt_index(d, d), 1), drop = FALSE]
  !all(piv > 1 / sqrt(.cont_qmax))
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
# singular (min eigenvalue / max eigenvalue of the equilibrated matrix, see
# .info_inverse()); negative beyond it: saddle point.
.pd_tol <- 1e-8

# Positive-definiteness classification and inverse of a symmetric information-type
# matrix that do not depend on the units of the parameters. With s = sqrt(|diag(I)|)
# (1 where the diagonal is 0), the equilibrated matrix C = I / (s s') has the same
# inertia as I (Sylvester's law), so a positive diagonal change of parameter units
# changes neither the classification nor the inverse, I^-1 = C^-1 / (s s'). The
# spectrum of C is compared with the relative tolerance: "ok" (positive definite),
# "saddle" (clearly negative eigenvalue) or "singular" (otherwise, including
# non-finite entries). Equilibration removes differences in units, not near linear
# dependence: an exact or numerical rank deficiency stays visible in C.
#
# 'scale' gives the reference scale explicitly (instead of the matrix's own diagonal).
# The reduced information is I_exact + (constraint curvature): where these cancel, the
# reduced diagonal is rounding and finite-difference noise, which its own diagonal would
# blow up to 1 and present as information. It is therefore scaled by the diagonal of
# the exact information, which transforms with the parameter units in the same way
# (s -> |d| s), so the decision stays unit-invariant while a reduced matrix that has
# lost (almost) all information in some direction relative to the exact one is
# reported as singular.
.info_inverse <- function(I, tol = .pd_tol, scale = NULL) {
  out <- list(status = "singular", eigen = c(min = NA_real_, max = NA_real_), s = NULL, Cinv = NULL, inverse = NULL)
  if (!is.matrix(I) || !length(I) || !all(is.finite(I))) return(out)
  I <- (I + t(I)) / 2
  s <- if (is.null(scale)) sqrt(abs(diag(I))) else scale
  if (length(s) != nrow(I) || !all(is.finite(s))) return(out)
  s[s == 0] <- 1
  e <- eigen(I / outer(s, s), symmetric = TRUE)
  ev <- e$values; out$eigen <- c(min = min(ev), max = max(ev))
  thr <- tol * max(abs(ev))
  if (min(ev) > thr) {
    Ci <- e$vectors %*% (t(e$vectors) / ev); Ci <- (Ci + t(Ci)) / 2
    V <- Ci / outer(s, s); dimnames(V) <- dimnames(I)
    out$status <- "ok"; out$s <- s; out$Cinv <- Ci; out$inverse <- V
  } else if (min(ev) < -thr) out$status <- "saddle"
  out
}

# central-difference step for each parameter: relative to the parameter, and for
# moderation effects relative to the scale of their moderator (a coefficient of a
# moderator measured in large units, e.g. age^2 in years, lives on a small scale)
.fd_steps <- function(theta, model) {
  base <- rep(1, length(theta)); q <- model$par$mod
  if (model$Q && any(!is.na(q))) {
    rms <- sqrt(colMeans(model$X^2))
    base[!is.na(q)] <- pmin(1, 1 / rms[q[!is.na(q)]])
  }
  1e-5 * pmax(base, abs(theta))
}

# A likelihood evaluation that the information matrix can be built from: a finite
# log-likelihood that carries its score matrix. The likelihood is returned as a bare
# -Inf without scores when the parameters are invalid or a person's log-likelihood is
# not finite (e.g. a residual variance at the limit of double precision).
.ll_ok <- function(ll) length(ll) == 1 && is.finite(ll) && !is.null(attr(ll, "scores"))

# Estimation failure: an error of class "mofia_estimation_error" (also used when the
# log-likelihood is not finite at the solution, see mofia()).
.estimation_error <- function(msg) stop(structure(class = c("mofia_estimation_error", "error", "condition"),
                                                  list(message = msg, call = NULL)))

# With only continuous items the log-likelihood is also known in closed form
# (.cont_loglik). The adaptive quadrature is exact up to rounding, but that rounding grows
# as a residual variance approaches zero and the posterior of the factors collapses (below
# about 1e-18 it can exceed .cont_ll_tol). A solution at which the two disagree is not
# returned: the gradient, the information and the standard errors come from the same
# quadrature.
.cont_ll_tol <- 1e-6
.check_closed_form <- function(ll, theta, model) {
  cf <- .cont_loglik(theta, model)
  if (is.finite(cf) && abs(ll - cf) <= .cont_ll_tol) return(invisible(cf))
  rv <- .engine_inputs(theta, model)$RV; k <- which(rv == min(rv), arr.ind = TRUE)[1, 2]
  .estimation_error(sprintf(paste0(
    "Estimation failed: at the solution the log-likelihood from quadrature (%.6f) differs from its exact value (%s), ",
    "so the estimates and standard errors cannot be trusted. The smallest residual variance is %.3g (item %s); as a ",
    "residual variance approaches zero (a Heywood case), the posterior of the factors becomes too concentrated for ",
    "the quadrature. The model may not be identified (check the anchor items), or the item may be (nearly) redundant."),
    ll, if (is.finite(cf)) sprintf("%.6f", cf) else "not computable", min(rv), model$items[k]))
}

# The information matrix cannot be formed because the likelihood is not finite at, or
# within a finite-difference step of, the estimates: they are at the edge of the region
# where the likelihood can be evaluated (a variance at the limit of double precision after
# a diverged optimization, or a singular factor correlation matrix), not at a regular
# maximum. No fit is returned (as when the log-likelihood is not finite at the solution).
.se_failed <- function(fit, why) {
  op <- fit$optimizer
  .estimation_error(paste0(
    "Estimation failed: ", why, ", so the information matrix cannot be computed",
    if (!is.null(op$message)) sprintf(" (the optimizer stopped with '%s' after %d outer iteration(s); max |gradient| %.3g)",
                                      op$message, as.integer(op$outer_iterations), fit$max_abs_grad),
    ". The estimates are at the edge of the region where the likelihood can be evaluated (e.g. a variance at the limit ",
    "of double precision, or a singular factor correlation matrix). The model may not be identified (check the anchor ",
    "items), or other starting values may help (starts > 1)."))
}

.mofia_se <- function(fit, check = TRUE) {
  model <- fit$model; theta <- fit$coefficients; link <- fit$link; adapt <- fit$adapt
  nodes <- .mofia_nodes(fit$nq, model$D)
  ll <- .mofia_ll(theta, model, link, nodes, scores = TRUE, adapt = adapt)
  if (!.ll_ok(ll)) .se_failed(fit, "the log-likelihood is not finite at the estimates")
  S <- attr(ll, "scores")
  # observed information: central differences of the analytic gradient. A step at which
  # the likelihood is not finite leaves the information undefined (see .se_failed).
  gfun <- function(th) {
    v <- .mofia_ll(th, model, link, nodes, scores = TRUE, adapt = adapt)
    if (.ll_ok(v)) colSums(attr(v, "scores")) else NULL
  }
  np <- length(theta); h <- .fd_steps(theta, model); H <- matrix(0, np, np); bad <- integer(0)
  for (k in seq_len(np)) {
    e <- replace(numeric(np), k, h[k])
    gp <- gfun(theta + e); gm <- gfun(theta - e)
    if (is.null(gp) || is.null(gm)) { bad <- c(bad, k); next }
    H[, k] <- (gp - gm) / (2 * h[k])
  }
  if (length(bad))
    .se_failed(fit, sprintf("the log-likelihood is not finite within a finite-difference step of the estimates of %s",
                            paste(sprintf("%s (%.6g)", names(theta)[bad], theta[bad]), collapse = ", ")))
  H <- (H + t(H)) / 2
  dn <- list(names(theta), names(theta))
  B <- crossprod(S); dimnames(B) <- dn
  I_exact <- -H
  W <- attr(ll, "phi_scores"); attr(W, "resvar") <- attr(ll, "resvar_scores")
  I_red <- -H + .constraint_curvature(theta, model, W)
  g <- colSums(S)

  # curvature: a maximum needs a positive-definite information matrix. The decision
  # uses the equilibrated matrix (independent of parameter units, see .info_inverse);
  # the raw spectrum is kept as a diagnostic.
  dimnames(I_exact) <- dimnames(I_red) <- dn
  ex <- .info_inverse(I_exact)
  ev <- if (all(is.finite(I_exact))) eigen(I_exact, symmetric = TRUE, only.values = TRUE)$values else NA_real_
  fit$information_eigen <- c(min = min(ev), max = max(ev))
  fit$information_eigen_scaled <- ex$eigen
  pd <- identical(ex$status, "ok")
  fit$status <- ex$status
  if (check && !pd) {
    fit$converged <- FALSE; fit$vcov <- NULL; fit$scaling <- NULL; fit$newton_decrement <- NA_real_
    fit$se_note <- if (fit$status == "saddle")
      sprintf("the solution is a saddle point, not a maximum (smallest eigenvalue of the scaled information matrix %.3g)", ex$eigen[["min"]])
    else sprintf(paste0("the information matrix is singular (smallest eigenvalue of the scaled information matrix %.3g, ",
                        "largest %.3g): the model is not identified, or not identified by these data"),
                 ex$eigen[["min"]], ex$eigen[["max"]])
    warning("Standard errors and tests are not available: ", fit$se_note, ". ",
            if (fit$status == "saddle") "Try other starting values (starts > 1)." else "Check the anchor items.",
            call. = FALSE)
    return(fit)
  }

  # sandwich, scaling factor and Newton decrement in equilibrated coordinates, then
  # back to the parameter units: A^-1 B A^-1 = (C^-1 (B/ss') C^-1) / ss'
  sandwich <- function(e) if (is.null(e$inverse)) NULL else {
    Ci <- e$Cinv; V <- Ci %*% (B / outer(e$s, e$s)) %*% Ci; V <- (V + t(V)) / 2 / outer(e$s, e$s)
    dimnames(V) <- dn; V
  }
  trace_AB <- function(e) if (is.null(e$inverse)) NA_real_ else sum(e$Cinv * (B / outer(e$s, e$s))) / np
  red <- .info_inverse(I_red, scale = sqrt(abs(diag(I_exact))))
  fit$vcov <- list(exact = list(ML = ex$inverse, MLR = sandwich(ex)),
                   reduced = if (is.null(red$inverse)) NULL else list(ML = red$inverse, MLR = sandwich(red)),
                   MLF = .info_inverse(B)$inverse)
  fit$scaling <- c(exact = trace_AB(ex), reduced = trace_AB(red))
  fit$information <- "exact"
  fit$newton_decrement <- if (is.null(ex$inverse)) NA_real_ else { gs <- g / ex$s; 0.5 * sum(gs * (ex$Cinv %*% gs)) }
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
