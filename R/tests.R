# Hypothesis tests ------------------------------------------------------------------

#' Likelihood-ratio test between two nested models
#'
#' Compares a restricted model `fit0` with a more general model `fit1`. With
#' `scaled = TRUE` the Satorra-Bentler (2001) scaled difference test for the
#' robust (MLR) estimator is computed: \eqn{c_d = (p_0 c_0 - p_1 c_1)/(p_0 - p_1)}
#' and \eqn{T_d = -2(\ell_0 - \ell_1)/c_d}, where \eqn{p} are the numbers of
#' parameters and \eqn{c} the scaling correction factors.
#'
#' @param fit0 The restricted (nested) fitted model.
#' @param fit1 The more general fitted model.
#' @param scaled Scaled (robust) test? Otherwise the ordinary LR test.
#' @param information `"default"` (exact), `"reduced"` or `"exact"`; see [mofia()].
#' @param sb2010 Also report the strictly positive Satorra-Bentler (2010)
#'   version, which evaluates the scaling factor of `fit1` at the estimates of
#'   `fit0` (one extra information computation). Requires models nested by
#'   parameter names (the same identification).
#' @param check_nested Stop if the parameter names of `fit0` are not a subset of
#'   those of `fit1`? Set to `FALSE` only for models that are nested in a way
#'   the names cannot show (e.g. the same model under marker and variance
#'   identification).
#' @return A data frame with the test statistic, degrees of freedom, scaling
#'   factors and p-value.
#' @details Both models must be fitted to the same cases and item responses,
#'   use the same link for categorical items, and have converged (positive
#'   definite information, see [mofia()]); every moderator of `fit0` must be a
#'   moderator of `fit1` with the same values. Otherwise the test stops with the
#'   reason. A clearly negative likelihood-ratio statistic (below -0.01) means
#'   that the models are not nested or that a fit did not reach its maximum, and
#'   is an error. Scaled statistics are asymptotic mean corrections; their
#'   finite-sample calibration is not guaranteed.
#' @references
#' Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
#' statistic for moment structure analysis. *Psychometrika, 66*, 507--514.
#' \doi{10.1007/BF02296192}
#'
#' Satorra, A., & Bentler, P. M. (2010). Ensuring positiveness of the scaled
#' difference chi-square test statistic. *Psychometrika, 75*, 243--248.
#' @export
mofia_lrt <- function(fit0, fit1, scaled = TRUE, information = c("default", "reduced", "exact"),
                      sb2010 = FALSE, check_nested = TRUE) {
  if (!inherits(fit0, "mofia") || !inherits(fit1, "mofia")) stop("Both models must be fitted 'mofia' objects.")
  if (fit0$npar >= fit1$npar) stop("'fit0' must have fewer parameters than 'fit1'.")
  .lrt_check(fit0, fit1)
  if (check_nested && !all(names(fit0$coefficients) %in% names(fit1$coefficients)))
    stop("'fit0' is not nested in 'fit1' (parameter names); use check_nested = FALSE to override.")
  p0 <- fit0$npar; p1 <- fit1$npar; L0 <- fit0$loglik; L1 <- fit1$loglik
  LR <- -2 * (L0 - L1); df <- p1 - p0
  if (LR < -0.01)
    stop(sprintf(paste0("The larger model has the lower log-likelihood (LR = %.3f): the models are not nested, ",
                        "or one fit did not reach its maximum (try starts > 1)."), LR), call. = FALSE)
  LR <- max(LR, 0)
  if (!scaled)
    return(data.frame(test = "Likelihood ratio", LR = LR, df = df, statistic = LR,
                      p_value = stats::pchisq(LR, df, lower.tail = FALSE)))
  info <- .info_choice(fit1, information)
  c0 <- fit0$scaling[[info]]; c1 <- fit1$scaling[[info]]
  if (is.null(c0) || is.null(c1) || !is.finite(c0) || !is.finite(c1))
    stop("Scaling factors are not available (fit with se = TRUE; the ", info, " information must be positive definite).")
  cd <- (p0 * c0 - p1 * c1) / (p0 - p1)
  out <- data.frame(test = "Satorra-Bentler (2001)", LR = LR, df = df, c0 = c0, c1 = c1, cd = cd,
                    statistic = LR / cd)
  if (sb2010) {
    if (!all(names(fit0$coefficients) %in% names(fit1$coefficients)))
      stop("The Satorra-Bentler (2010) test needs the restricted estimates in the larger model's parameterization; ",
           "this is only implemented for models nested by parameter names. Refit both models with the same ",
           "identification (the test statistic does not depend on it).", call. = FALSE)
    th10 <- stats::setNames(numeric(p1), names(fit1$coefficients))
    th10[names(fit0$coefficients)] <- fit0$coefficients
    c10 <- .scaling_at(fit1, th10)[[info]] %||% NA_real_
    cd10 <- (p0 * c0 - p1 * c10) / (p0 - p1)
    out <- rbind(out, data.frame(test = "Satorra-Bentler (2010)", LR = LR, df = df, c0 = c0, c1 = c10,
                                 cd = cd10, statistic = LR / cd10))
  }
  out$p_value <- ifelse(is.finite(out$cd) & out$cd > 0, stats::pchisq(out$statistic, out$df, lower.tail = FALSE), NA)
  out
}

# Can two fits be compared by a likelihood-ratio test? Stops with the reason if not.
.lrt_check <- function(fit0, fit1, nm = c("fit0", "fit1")) {
  m0 <- fit0$model; m1 <- fit1$model; why <- NULL
  if (!setequal(m0$items, m1$items) || length(m0$items) != length(m1$items)) why <- "they use different items"
  else if (m0$N != m1$N) why <- "they were fitted to different numbers of cases"
  else {
    it <- m1$items
    Y0 <- m0$Y[, it, drop = FALSE]; Y1 <- m1$Y[, it, drop = FALSE]
    if (!identical(unname(m0$cont[it]), unname(m1$cont[it])))
      why <- "the items are treated differently (continuous vs categorical)"
    else if (!identical(is.na(Y0), is.na(Y1)) || any(Y0 != Y1, na.rm = TRUE) ||
             !identical(m0$levels[match(it, m0$items)], m1$levels))
      why <- "the item responses differ"
    else {
      extra <- setdiff(m0$moderators, m1$moderators)
      if (length(extra)) why <- paste0("the smaller model uses moderator(s) the larger one lacks: ", paste(extra, collapse = ", "))
      else if (length(m0$moderators)) {
        X0 <- m0$X; X1 <- m1$X[, m0$moderators, drop = FALSE]
        if (max(abs(X0 - X1)) > 1e-8 * max(1, max(abs(X1))))
          why <- "the moderator values differ (e.g. standardized in only one model)"
      }
    }
  }
  if (is.null(why) && any(!m1$cont) && !identical(fit0$link, fit1$link))
    why <- "they use different links, so their likelihoods are not nested"
  if (!is.null(why)) stop(sprintf("'%s' and '%s' cannot be compared: %s.", nm[1], nm[2], why), call. = FALSE)
  fits <- list(fit0, fit1)
  for (k in 1:2) if (!isTRUE(fits[[k]]$converged))
    stop(sprintf("'%s' has not converged%s; a likelihood-ratio test needs the maximum of both models.", nm[k],
                 if (!is.null(fits[[k]]$se_note)) paste0(" (", fits[[k]]$se_note, ")") else ""), call. = FALSE)
  if (!identical(fit0$nq, fit1$nq))
    warning(sprintf("'%s' and '%s' use different numbers of quadrature nodes; their log-likelihoods differ slightly in accuracy.",
                    nm[1], nm[2]), call. = FALSE)
  invisible(TRUE)
}

# scaling correction factors of a fitted model's specification at other parameter values
.scaling_at <- function(fit, theta) {
  nodes <- .mofia_nodes(fit$nq, fit$model$D)
  pseudo <- fit; pseudo$coefficients <- theta
  pseudo$adapt <- .mofia_adapt(theta, fit$model, fit$link, nodes)
  .mofia_se(pseudo)$scaling
}

#' Wald test that a set of parameters is zero
#'
#' @param fit A fitted `mofia` object.
#' @param params Character vector of parameter names, or (with `regex = TRUE`)
#'   a regular expression selecting parameters, e.g. `"^beta\\[DevilPossession,"`.
#' @param se Covariance type: `"MLR"` (default), `"ML"` or `"MLF"`.
#' @param information `"default"` (exact), `"reduced"` or `"exact"`; see [mofia()].
#' @param regex Treat `params` as a regular expression?
#' @return A data frame with the Wald statistic, degrees of freedom, p-value and
#'   the tested parameters. The fit must have converged; otherwise the test is
#'   refused.
#' @examples
#' \donttest{
#' data(parabelief)
#' ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession",
#'         "Sorcery", "HauntedPlaces")
#' fit <- mofia(mofia_model(parabelief, items = ub, moderators = c("Age", "Male"),
#'                          threshold_mod = "DevilPossession"))
#' mofia_wald(fit, "^beta\\[DevilPossession,", regex = TRUE)
#' }
#' @export
mofia_wald <- function(fit, params, se = c("MLR", "ML", "MLF"), information = c("default", "reduced", "exact"),
                       regex = FALSE) {
  se <- match.arg(se)
  if (!isTRUE(fit$converged))
    stop("The fit has not converged; a Wald test needs the maximum-likelihood estimates.", call. = FALSE)
  nms <- names(fit$coefficients)
  idx <- if (regex) grep(params, nms) else match(params, nms)
  if (regex && !length(idx)) stop("No parameter matches '", params, "'.")
  if (anyNA(idx)) stop("Unknown parameter(s): ", paste(params[is.na(idx)], collapse = ", "))
  V <- vcov(fit, se, information)
  th <- fit$coefficients[idx]
  Vs <- V[idx, idx, drop = FALSE]
  # checked and solved on the correlation scale, so that the decision does not depend
  # on the units of the tested parameters; the statistic is unit-free anyway
  sv <- sqrt(diag(Vs))
  if (!all(is.finite(Vs)) || !all(is.finite(sv) & sv > 0))
    stop("The covariance matrix of the tested parameters is not positive definite; no Wald test.", call. = FALSE)
  R <- Vs / outer(sv, sv); R <- (R + t(R)) / 2
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-12 * max(abs(ev)))
    stop("The covariance matrix of the tested parameters is not positive definite; no Wald test.", call. = FALSE)
  W <- drop(t(th / sv) %*% solve(R, th / sv))
  data.frame(statistic = W, df = length(idx), p_value = stats::pchisq(W, length(idx), lower.tail = FALSE),
             parameters = paste(nms[idx], collapse = ", "))
}
