# Factor scores and model-implied predictions ---------------------------------------

#' Factor scores
#'
#' Expected a posteriori (EAP) factor scores: the posterior mean of each
#' person's factor values given their responses and moderator values, with the
#' posterior standard deviation as its standard error. The posterior is taken
#' at the estimated parameters, so the SD does not include the uncertainty of
#' the parameter estimates. Scores are on the scale
#' of the model's factors (mean 0 and variance 1 for a person with all
#' moderators at 0, under variance identification), so they include the
#' moderated factor means and variances.
#'
#' @param object A fitted `mofia` object.
#' @return A data frame with one row per fitted case (in the order of the
#'   fitted data) and, per factor, its EAP score and posterior SD
#'   (`<factor>_sd`). Cases removed by [mofia_model()] (missing moderator
#'   values, no observed response) have no row. The row names are those of the
#'   fitted rows of the original data, and attribute `"row_index"` gives their
#'   positions in it (`model$cases`), so scores can be matched to the data by
#'   index rather than by position. Start from a full-length column of `NA`
#'   (also when reusing an existing column), so that excluded cases stay
#'   missing:
#'   ```
#'   data$score <- rep(NA_real_, nrow(data))
#'   data$score[attr(s, "row_index")] <- s$eta
#'   ```
#'   For models created with MoFIA 0.1.0, which did not record the original
#'   positions, `"row_index"` is `NULL` (unknown) and the row names are those
#'   kept in the model, if any. For several factors, the posterior covariance
#'   matrices are attached as attribute `"cov"` (persons x factors x factors,
#'   the first dimension named like the rows).
#' @examples
#' \donttest{
#' data(parabelief)
#' ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession",
#'         "Sorcery", "HauntedPlaces")
#' fit <- mofia(mofia_model(parabelief, items = ub, moderators = c("Age", "Male")))
#' head(mofia_scores(fit))
#' }
#' @export
mofia_scores <- function(object) {
  if (!inherits(object, "mofia")) stop("'object' must be a fitted mofia object.")
  model <- object$model; theta <- object$coefficients; D <- model$D; N <- model$N
  nodes <- .mofia_nodes(object$nq, D)
  ll <- .mofia_ll(theta, model, object$link, nodes, per_person = TRUE, adapt = object$adapt, return_post = TRUE)
  post <- .adapt_compat(attr(ll, "adapt"), D)               # z | y: mean m, Cholesky Lc
  pr <- .unpack(theta, model); st <- .structural(pr, model)
  Cf <- function(d, e) if (e > d) 0 else st$s[, d] * (if (D == 1) 1 else st$Lr[[d]][[e]])   # eta = mu + C z
  Lz <- function(d, e) if (e > d) 0 else post$Lc[, .lt_index(d, e)]
  mean <- matrix(0, N, D); M <- array(0, c(N, D, D))
  for (d in seq_len(D)) {
    mean[, d] <- st$mu[, d] + Reduce(`+`, lapply(seq_len(d), function(k) Cf(d, k) * post$m[, k]))
    for (e in seq_len(d)) M[, d, e] <- Reduce(`+`, lapply(e:d, function(k) Cf(d, k) * Lz(k, e)))
  }
  ids <- .case_ids(model)
  cov <- array(0, c(N, D, D), dimnames = list(ids$names, model$factors, model$factors))
  for (d in seq_len(D)) for (f in seq_len(D)) cov[, d, f] <- rowSums(M[, d, , drop = FALSE] * M[, f, , drop = FALSE])
  sds <- vapply(seq_len(D), function(d) sqrt(pmax(cov[, d, d], 0)), numeric(N))
  out <- data.frame(mean, matrix(sds, N, D))
  names(out) <- c(model$factors, paste0(model$factors, "_sd"))
  rownames(out) <- ids$names
  attr(out, "row_index") <- ids$index
  if (D > 1) attr(out, "cov") <- cov
  out
}

# identities of the fitted cases: positions in and row names of the original data
# (see below for models from MoFIA 0.1.0, which have no $cases)
.case_ids <- function(model) {
  cs <- model$cases
  if (!is.null(cs)) return(list(index = cs$index, names = cs$names))
  # models from MoFIA 0.1.0 (no $cases): the original positions of the fitted cases are
  # not known after exclusions, so none are claimed; row names kept in $Y are used
  list(index = NULL, names = rownames(model$Y))
}

# Marginal (over the factors) moments of each item's linear predictor given the
# moderators: a_ij = nu_ij + lambda_ij' eta_i with eta_i ~ N(mu_i, Sigma_i), so
# a_ij ~ N(m_ij, s_ij^2), m = nu + lambda' mu, s^2 = lambda' Sigma lambda = |C' lambda|^2.
.marginal_lp <- function(theta, model, X) {
  sp <- model; sp$X <- X; sp$N <- nrow(X)
  ei <- .engine_inputs(theta, sp)
  if (is.null(ei)) stop("The parameters imply an invalid factor correlation for some moderator values.")
  N <- nrow(X); P <- model$P; D <- model$D
  Cf <- function(d, e) if (e > d) 0 else ei$CL[, .lt_index(d, e)]
  m <- s <- matrix(0, N, P)
  for (j in seq_len(P)) {
    lam <- matrix(sapply(seq_len(D), function(d) ei$LAM[, j + P * (d - 1)]), N, D)
    m[, j] <- ei$NU[, j] + rowSums(lam * ei$st$mu)
    v2 <- 0
    for (e in seq_len(D)) v2 <- v2 + Reduce(`+`, lapply(e:D, function(d) lam[, d] * Cf(d, e)))^2
    s[, j] <- sqrt(v2)
  }
  list(m = m, s = s, pr = ei$pr)
}

# E[plogis(mu + s Z)], Z ~ N(0, 1): Gauss-Hermite over Z when s <= 4, otherwise over
# the logistic variable L (E[pnorm((mu - L) / s)]) with Gauss-Legendre on its CDF scale;
# absolute error about 1e-6 or less
.logit_normal <- function(mu, s) {
  out <- numeric(length(mu)); small <- s <= 4
  if (any(small)) {
    gh <- statmod::gauss.quad.prob(101, dist = "normal")
    out[small] <- drop(stats::plogis(outer(mu[small], rep(1, 101)) + outer(s[small], gh$nodes)) %*% gh$weights)
  }
  if (any(!small)) {
    gl <- statmod::gauss.quad(200, kind = "legendre"); u <- (gl$nodes + 1) / 2
    out[!small] <- drop(stats::pnorm(outer(mu[!small], stats::qlogis(u), "-") / s[!small]) %*% (gl$weights / 2))
  }
  out
}

# marginal category probabilities of a categorical item: P(y >= k) = E F(a - tau_k)
.marginal_probs <- function(m, s, tau, link) {
  cum <- lapply(tau, function(t)
    if (link == "probit") stats::pnorm((m - t) / sqrt(1 + s^2)) else .logit_normal(m - t, s))
  cum <- c(list(rep(1, length(m))), cum, list(rep(0, length(m))))
  lapply(seq_len(length(tau) + 1), function(k) pmax(cum[[k]] - cum[[k + 1]], 0))
}

# category probabilities of item j given a matrix of linear predictors (persons x nodes):
# P(y = k) = F(A - tau_{k-1}) - F(A - tau_k), evaluated in the tail where the
# difference does not cancel (1 - F(c) - (1 - F(b)) = F(-c) - F(-b) in the upper half)
.cat_probs <- function(A, tau, link) {
  Fc <- .link_funs(link)$F
  te <- c(-Inf, tau, Inf)
  lapply(seq_len(length(tau) + 1), function(k) {
    b <- A - te[k]; c <- A - te[k + 1]
    up <- (b + c) > 0
    out <- Fc(b) - Fc(c)
    out[up] <- Fc(-c[up]) - Fc(-b[up])
    out
  })
}

# linear predictor of item j at the given factor values (list of matrices)
.item_lp <- function(pr, j, X, eta, loadings) {
  A <- drop(X %*% pr$B[j, ])
  out <- matrix(A, nrow(eta[[1]]), ncol(eta[[1]]))
  for (d in seq_along(eta)) if (loadings[j, d]) out <- out + (pr$L0[j, d] + drop(X %*% pr$L1[j, d, ])) * eta[[d]]
  out
}

#' Model-implied predictions
#'
#' Category probabilities and expected item scores implied by a fitted model
#' for given moderator values, averaged over the factors' moderated
#' distribution (not conditional on a person's responses), or EAP factor
#' scores for the fitted data.
#'
#' @param object A fitted `mofia` object.
#' @param newdata Optional data frame with the moderator variables (on their
#'   original scale). Default: the fitted data.
#' @param type `"expected"` (expected item scores, in the items' original
#'   coding), `"prob"` (category probabilities) or `"scores"` (EAP factor
#'   scores for the fitted data; see [mofia_scores()]).
#' @param nq Not used (kept for compatibility). Each item's linear predictor is
#'   normal given the moderators, so probit probabilities are exact and logit
#'   probabilities use an accurate one-dimensional rule.
#' @param ... Not used.
#' @return For `"expected"`, a persons x items matrix. For `"prob"`, a list with
#'   one persons x categories matrix per categorical item (continuous items are
#'   omitted). For `"scores"`, a data frame. The rows are the fitted cases
#'   (without `newdata`) or the rows of `newdata`, in order; their row names
#'   are carried over, and attribute `"row_index"` gives the positions in the
#'   original data (fitted cases, see [mofia_scores()]; `NULL` for models
#'   created with MoFIA 0.1.0) or in `newdata`.
#' @examples
#' \donttest{
#' data(parabelief)
#' ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession")
#' fit <- mofia(mofia_model(parabelief, items = ub, moderators = "Male",
#'                          threshold_mod = "DevilPossession"))
#' predict(fit, newdata = data.frame(Male = c(0, 1)))
#' predict(fit, newdata = data.frame(Male = c(0, 1)), type = "prob")$DevilPossession
#' }
#' @export
predict.mofia <- function(object, newdata = NULL, type = c("expected", "prob", "scores"), nq = NULL, ...) {
  type <- match.arg(type)
  model <- object$model
  if (type == "scores") {
    if (!is.null(newdata)) stop("Factor scores are available for the fitted data only.")
    return(mofia_scores(object))
  }
  if (is.null(newdata)) {
    X <- model$X; ids <- .case_ids(model)
  } else {
    newdata <- as.data.frame(newdata)
    X <- .design_X(model, newdata); ids <- list(index = seq_len(nrow(newdata)), names = rownames(newdata))
  }
  ml <- .marginal_lp(object$coefficients, model, X); pr <- ml$pr
  probs <- list(); expd <- matrix(NA_real_, nrow(X), model$P, dimnames = list(ids$names, model$items))
  for (j in seq_len(model$P)) {
    if (model$cont[j]) { expd[, j] <- ml$m[, j]; next }
    pk <- vapply(.marginal_probs(ml$m[, j], ml$s[, j], pr$tau[[j]], object$link), identity, numeric(nrow(X)))
    pk <- matrix(pk, nrow(X), dimnames = list(ids$names, model$levels[[j]]))
    attr(pk, "row_index") <- ids$index
    probs[[model$items[j]]] <- pk
    expd[, j] <- drop(pk %*% model$levels[[j]])
  }
  attr(expd, "row_index") <- ids$index
  if (type == "prob") probs else expd
}

#' @rdname predict.mofia
#' @export
fitted.mofia <- function(object, ...) predict.mofia(object, type = "expected")

#' Confidence intervals for model parameters
#'
#' Wald intervals based on the chosen standard errors.
#'
#' @param object A fitted `mofia` object.
#' @param parm Parameter names (or indices); default all.
#' @param level Confidence level.
#' @param se,information Standard-error type and information; see [vcov.mofia()].
#' @param ... Not used.
#' @return A matrix with lower and upper limits (`NA` where the standard error
#'   is not available).
#' @export
confint.mofia <- function(object, parm, level = 0.95, se = c("MLR", "ML", "MLF"),
                          information = c("default", "reduced", "exact"), ...) {
  se <- match.arg(se)
  if (!isTRUE(object$converged))
    warning("The fit has not converged; these intervals refer to an unfinished solution.", call. = FALSE)
  est <- object$coefficients
  s <- .se_from_vcov(vcov(object, se, information))
  if (!missing(parm)) { est <- est[parm]; s <- s[parm] }
  a <- (1 - level) / 2; z <- stats::qnorm(1 - a)
  out <- cbind(est - z * s, est + z * s)
  colnames(out) <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%")
  out
}

#' Tidy and glance methods
#'
#' Parameter table and one-row model summary in the format of the broom
#' package (available when the generics or broom package is loaded).
#'
#' @param x A fitted `mofia` object.
#' @param conf.int Add confidence intervals?
#' @param conf.level Confidence level.
#' @param se,information Standard-error type and information; see [vcov.mofia()].
#' @param ... Not used.
#' @return A data frame.
#' @name mofia_tidiers
NULL

#' @rdname mofia_tidiers
#' @exportS3Method generics::tidy
tidy.mofia <- function(x, conf.int = FALSE, conf.level = 0.95, se = c("MLR", "ML", "MLF"),
                       information = c("default", "reduced", "exact"), ...) {
  se <- match.arg(se)
  s <- summary.mofia(x, se = se, information = information)$coefficients
  out <- data.frame(term = rownames(s), type = s$type, estimate = s$estimate, std.error = s$se,
                    statistic = s$z, p.value = s$p, stringsAsFactors = FALSE)
  if (conf.int) {
    ci <- confint.mofia(x, level = conf.level, se = se, information = information)
    out$conf.low <- ci[, 1]; out$conf.high <- ci[, 2]
  }
  rownames(out) <- NULL
  out
}

#' @rdname mofia_tidiers
#' @exportS3Method generics::glance
glance.mofia <- function(x, ...) {
  data.frame(logLik = x$loglik, AIC = stats::AIC(x), BIC = stats::BIC(x), npar = x$npar, nobs = x$nobs,
             converged = isTRUE(x$converged), newton_decrement = x$newton_decrement,
             scaling = if (is.null(x$scaling)) NA_real_ else x$scaling[[x$information %||% "exact"]],
             status = x$status %||% NA_character_,
             link = x$link)
}
