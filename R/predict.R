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
#' @return A data frame with one row per person (in the order of the fitted
#'   data) and, per factor, its EAP score and posterior SD (`<factor>_sd`). For
#'   several factors, the posterior covariance matrices are attached as
#'   attribute `"cov"` (persons x factors x factors).
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
  post <- attr(ll, "adapt")                                  # z | y: mean m, Cholesky Lc
  pr <- .unpack(theta, model); st <- .structural(pr, model)
  Cf <- function(d, e) if (e > d) 0 else st$s[, d] * (if (D == 1) 1 else st$Lr[[d]][[e]])   # eta = mu + C z
  Lz <- function(d, e) if (e > d) 0 else post$Lc[, .lt_index(d, e)]
  mean <- matrix(0, N, D); M <- array(0, c(N, D, D))
  for (d in seq_len(D)) {
    mean[, d] <- st$mu[, d] + Reduce(`+`, lapply(seq_len(d), function(k) Cf(d, k) * post$m[, k]))
    for (e in seq_len(d)) M[, d, e] <- Reduce(`+`, lapply(e:d, function(k) Cf(d, k) * Lz(k, e)))
  }
  cov <- array(0, c(N, D, D), dimnames = list(NULL, model$factors, model$factors))
  for (d in seq_len(D)) for (f in seq_len(D)) cov[, d, f] <- rowSums(M[, d, , drop = FALSE] * M[, f, , drop = FALSE])
  sds <- vapply(seq_len(D), function(d) sqrt(pmax(cov[, d, d], 0)), numeric(N))
  out <- data.frame(mean, matrix(sds, N, D))
  names(out) <- c(model$factors, paste0(model$factors, "_sd"))
  if (D > 1) attr(out, "cov") <- cov
  out
}

# person x node quantities under the factors' moderated (prior) distribution
.prior_eta <- function(theta, model, X, nq) {
  sp <- model; sp$X <- X; sp$N <- nrow(X)
  pr <- .unpack(theta, sp); st <- .structural(pr, sp)
  if (!st$ok) stop("The parameters imply an invalid factor correlation for some moderator values.")
  nodes <- .mofia_nodes(nq, model$D); D <- model$D
  eta <- lapply(seq_len(D), function(d) {
    Reduce(`+`, lapply(seq_len(d), function(e) outer(st$s[, d] * (if (D == 1) 1 else st$Lr[[d]][[e]]), nodes$T[, e])),
           init = matrix(st$mu[, d], nrow(X), nrow(nodes$T)))
  })
  list(pr = pr, eta = eta, w = exp(nodes$logw))
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
#' @param nq Quadrature nodes per factor for the integration.
#' @param ... Not used.
#' @return For `"expected"`, a persons x items matrix. For `"prob"`, a list with
#'   one persons x categories matrix per categorical item (continuous items are
#'   omitted). For `"scores"`, a data frame.
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
  X <- if (is.null(newdata)) model$X else .design_X(model, newdata)
  if (is.null(nq)) nq <- c(41, 21, 11, 7, 5, 4)[min(model$D, 6)]
  pe <- .prior_eta(object$coefficients, model, X, nq); pr <- pe$pr
  probs <- list(); expd <- matrix(NA_real_, nrow(X), model$P, dimnames = list(NULL, model$items))
  for (j in seq_len(model$P)) {
    if (model$cont[j]) {
      mu <- pr$nu0[j] + drop(X %*% pr$B[j, ])
      for (d in seq_len(model$D)) if (model$loadings[j, d])
        mu <- mu + (pr$L0[j, d] + drop(X %*% pr$L1[j, d, ])) * drop(pe$eta[[d]] %*% pe$w)
      expd[, j] <- mu
      next
    }
    A <- .item_lp(pr, j, X, pe$eta, model$loadings)
    pk <- vapply(.cat_probs(A, pr$tau[[j]], object$link), function(p) drop(p %*% pe$w), numeric(nrow(X)))
    pk <- matrix(pk, nrow(X)); colnames(pk) <- model$levels[[j]]
    probs[[model$items[j]]] <- pk
    expd[, j] <- drop(pk %*% model$levels[[j]])
  }
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
