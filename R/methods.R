# Methods for fitted models ----------------------------------------------------------

.info_choice <- function(object, information) {
  information <- match.arg(information, c("default", "reduced", "exact"))
  if (information == "default") object$information %||% "exact" else information
}

# standard errors from a covariance matrix: NA (not 0) where a variance is not positive
.se_from_vcov <- function(V) {
  d <- diag(V)
  ifelse(is.finite(d) & d > 0, sqrt(pmax(d, 0)), NA_real_)
}

#' Covariance matrix of the estimates
#'
#' @param object A fitted `mofia` object.
#' @param se `"MLR"` (robust sandwich, default), `"ML"` (inverse observed
#'   information) or `"MLF"` (outer product of the scores).
#' @param information `"default"` (the exact information), `"exact"` or
#'   `"reduced"` observed information; see [mofia()].
#' @param ... Not used.
#' @return A covariance matrix. An error is raised when it is not available
#'   (no standard errors computed, or an information matrix that is not
#'   positive definite).
#' @export
vcov.mofia <- function(object, se = c("MLR", "ML", "MLF"), information = c("default", "reduced", "exact"), ...) {
  se <- match.arg(se)
  if (is.null(object$vcov))
    stop("No standard errors available: ", object$se_note %||% "the model was fitted with se = FALSE", ".", call. = FALSE)
  info <- if (se == "MLF") "MLF" else .info_choice(object, information)
  V <- if (se == "MLF") object$vcov$MLF else object$vcov[[info]][[se]]
  if (is.null(V)) stop("The ", if (se == "MLF") "outer product of the scores" else paste(info, "information matrix"),
                       " is not positive definite; these standard errors are not available.", call. = FALSE)
  V
}

#' @export
coef.mofia <- function(object, ...) object$coefficients

#' @export
logLik.mofia <- function(object, ...)
  structure(object$loglik, df = object$npar, nobs = object$nobs, class = "logLik")

#' @export
nobs.mofia <- function(object, ...) object$nobs

#' @export
print.mofia <- function(x, digits = 3, ...) {
  m <- x$model; nc <- sum(m$cont)
  cat("Moderated factor model fitted by marginal maximum likelihood (MoFIA)\n")
  cat(sprintf("  %d item(s) (%s), %d factor(s), moderators: %s\n", m$P,
              paste(c(if (m$P - nc) sprintf("%d categorical", m$P - nc), if (nc) sprintf("%d continuous", nc)), collapse = ", "),
              m$D, if (m$Q) paste(m$moderators, collapse = ", ") else "none"))
  cat(sprintf("  N = %d, %d parameters%s, %d^%d adaptive quadrature nodes%s\n", x$nobs, x$npar,
              if (m$P > nc) paste0(", ", x$link, " link") else "", x$nq, m$D,
              if (nc == m$P) " (exact for continuous items)"
              else if (!is.null(x$quadrature) && nrow(x$quadrature) >= 2) {
                q <- x$quadrature; k <- nrow(q)
                sprintf(" (check: logLik changes by %.1e with %d nodes)", q$logLik[k] - q$logLik[k - 1], q$nq[k])
              } else ""))
  cat(sprintf("  logLik = %.3f, AIC = %.2f, BIC = %.2f\n", x$loglik, stats::AIC(x), stats::BIC(x)))
  status <- x$status %||% "ok"
  cat(sprintf("  %s%s\n", if (isTRUE(x$converged)) "Converged" else "NOT converged",
              if (status == "saddle") ": saddle point (information matrix not positive definite)"
              else if (status == "singular") ": information matrix singular (model not identified by these data)"
              else if (is.finite(x$newton_decrement)) sprintf(" (Newton decrement %.1e)", x$newton_decrement) else ""))
  invisible(x)
}

#' Summarize a fitted moderated factor model
#'
#' Parameter estimates with standard errors, z statistics and p-values,
#' grouped by parameter type.
#'
#' @param object A fitted `mofia` object.
#' @param se Standard-error type: `"MLR"` (default), `"ML"` or `"MLF"`.
#' @param information `"default"` (exact), `"reduced"` or `"exact"`; see [mofia()].
#' @param ... Not used.
#' @return An object of class `summary.mofia`; its `coefficients` element is a
#'   data frame with columns `estimate`, `se`, `z`, `p` and `type`. Standard
#'   errors that are not available (no positive-definite information matrix, or
#'   a non-positive variance) are `NA`, as are their z statistics and p-values.
#' @export
summary.mofia <- function(object, se = c("MLR", "ML", "MLF"), information = c("default", "reduced", "exact"), ...) {
  se <- match.arg(se)
  th <- object$coefficients
  V <- if (is.null(object$vcov)) NULL else tryCatch(vcov(object, se, information), error = function(e) NULL)
  s <- if (is.null(V)) rep(NA_real_, length(th)) else .se_from_vcov(V)
  z <- th / s
  tab <- data.frame(estimate = th, se = s, z = z, p = 2 * stats::pnorm(-abs(z)),
                    type = object$model$par$block, row.names = names(th), stringsAsFactors = FALSE)
  structure(list(coefficients = tab, fit = object, se_type = se,
                 information = if (se == "MLF") NA else .info_choice(object, information)),
            class = "summary.mofia")
}

#' @export
print.summary.mofia <- function(x, digits = 3, ...) {
  print(x$fit)
  cat(sprintf("\nStandard errors: %s%s\n", x$se_type,
              if (!is.na(x$information)) sprintf(" (%s information)", x$information) else ""))
  if (!is.null(x$fit$scaling) && !is.na(x$information))
    cat(sprintf("Robust scaling correction factor: %.4f\n", x$fit$scaling[[x$information]]))
  if (!is.null(x$fit$vcov) && !isTRUE(x$fit$converged))
    cat("WARNING: the fit has not converged. Estimates and standard errors refer to an unfinished solution;",
        "do not use them for inference.\n")
  if (is.null(x$fit$vcov)) cat("Standard errors not available: ", x$fit$se_note %||% "fitted with se = FALSE", ".\n", sep = "")
  else if (anyNA(x$coefficients$se)) cat("Some standard errors are not available (non-positive variance estimates).\n")
  groups <- list(tau = "Thresholds", nu0 = "Intercepts", lambda0 = "Loadings (baseline)",
                 lambda1 = "Loading moderation", beta = "Threshold/intercept moderation (DIF)",
                 omega0 = "Log residual variances (baseline)", omega = "Residual-variance moderation",
                 gamma = "Factor-mean moderation", delta0 = "Log factor variances (baseline)",
                 delta = "Factor-variance moderation", r0 = "Factor correlations (baseline)",
                 r1 = "Factor-correlation moderation")
  if (identical(x$fit$model$cov_param, "linear")) { groups$r0 <- "Factor covariances (baseline)"; groups$r1 <- "Factor-covariance moderation" }
  tab <- x$coefficients
  for (b in names(groups)) {
    sel <- tab$type == b; if (!any(sel)) next
    cat(sprintf("\n%s:\n", groups[[b]]))
    out <- tab[sel, c("estimate", "se", "z", "p")]
    rownames(out) <- sub("^[a-z0-9]+\\[(.*)\\]$", "\\1", rownames(out))
    out$p <- format.pval(out$p, digits = 2, eps = 1e-4)
    out[, 1:3] <- lapply(out[, 1:3], function(v) formatC(v, digits = digits, format = "f"))
    print(out, right = TRUE)
  }
  if (identical(x$fit$model$cov_param, "fisherz") && any(tab$type %in% c("r0", "r1")))
    cat("\nFactor correlations are on the Fisher-z scale: rho = tanh(zcor0 + x'zcor1).\n")
  invisible(x)
}

#' Compare nested moderated factor models
#'
#' Likelihood-ratio tests between fitted models, by default with the
#' Satorra-Bentler (2001) scaling correction for the robust (MLR) estimator.
#' Models are sorted by their number of parameters and compared in sequence.
#'
#' @param object,... Fitted `mofia` objects (at least two).
#' @param scaled Use the scaled difference test? If `FALSE`, the ordinary
#'   likelihood-ratio test is reported.
#' @param information `"default"` (exact), `"reduced"` or `"exact"` information
#'   for the scaling factors; see [mofia()].
#' @return A data frame of class `anova` with one row per model.
#' @details All models must be fitted to the same cases and item responses
#'   with the same link and must have converged; see [mofia_lrt()]. Nesting is
#'   checked by parameter names: every parameter of the smaller model must also
#'   be in the larger one. Models that fail this check are still compared, with
#'   a warning, because some models are nested in a way the names cannot show
#'   (e.g. marker vs variance identification); the test is only valid for truly
#'   nested models. A pair whose larger model has a clearly lower
#'   log-likelihood gets no test and a warning.
#'   The scaled test can be negative-definite in small samples (scaling
#'   factor below 0); the p-value is then `NA` and [mofia_lrt()] with
#'   `sb2010 = TRUE` gives a strictly positive alternative.
#' @references
#' Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
#' statistic for moment structure analysis. *Psychometrika, 66*, 507--514.
#' \doi{10.1007/BF02296192}
#' @export
anova.mofia <- function(object, ..., scaled = TRUE, information = c("default", "reduced", "exact")) {
  fits <- c(list(object), list(...))
  nm <- vapply(as.list(substitute(list(object, ...)))[-1], function(e) paste(deparse(e), collapse = ""), "")
  if (length(fits) < 2) stop("Give at least two fitted models.")
  if (!all(vapply(fits, inherits, TRUE, "mofia"))) stop("All models must be fitted 'mofia' objects.")
  ord <- order(vapply(fits, `[[`, 1, "npar")); fits <- fits[ord]; nm <- nm[ord]
  if (anyDuplicated(vapply(fits, `[[`, 1, "npar"))) stop("Two models have the same number of parameters; they cannot be nested.")
  n <- length(fits)
  out <- data.frame(npar = vapply(fits, `[[`, 1, "npar"), logLik = vapply(fits, `[[`, 1, "loglik"),
                    AIC = vapply(fits, stats::AIC, 1), BIC = vapply(fits, stats::BIC, 1),
                    Chisq = NA_real_, Df = NA_real_, Scaling = NA_real_, `Pr(>Chisq)` = NA_real_,
                    row.names = nm, check.names = FALSE)
  for (i in seq_len(n)[-1]) {
    .lrt_check(fits[[i - 1]], fits[[i]], nm[c(i - 1, i)])
    if (!all(names(fits[[i - 1]]$coefficients) %in% names(fits[[i]]$coefficients)))
      warning(sprintf("'%s' is not nested in '%s' by parameter names; the test is only valid for nested models.",
                      nm[i - 1], nm[i]), call. = FALSE)
    LR <- -2 * (fits[[i - 1]]$loglik - fits[[i]]$loglik)
    if (LR < -0.01) {
      warning(sprintf("'%s' has a lower log-likelihood than the smaller model '%s' (LR = %.3f): not nested, or a fit did not reach its maximum; no test.",
                      nm[i], nm[i - 1], LR), call. = FALSE)
      next
    }
    tst <- mofia_lrt(fits[[i - 1]], fits[[i]], scaled = scaled, information = information, check_nested = FALSE)
    out$Chisq[i] <- tst$statistic; out$Df[i] <- tst$df; out$`Pr(>Chisq)`[i] <- tst$p_value
    if (scaled) out$Scaling[i] <- tst$cd
  }
  if (!scaled) out$Scaling <- NULL
  structure(out, heading = c(if (scaled) "Scaled likelihood-ratio tests (Satorra-Bentler 2001, MLR)\n"
                             else "Likelihood-ratio tests\n"), class = c("anova", "data.frame"))
}
