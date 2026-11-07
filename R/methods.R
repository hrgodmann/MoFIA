# Methods for fitted models ----------------------------------------------------------

.info_choice <- function(object, information) {
  information <- match.arg(information, c("default", "reduced", "exact"))
  if (information == "default") object$information %||% "reduced" else information
}

#' Covariance matrix of the estimates
#'
#' @param object A fitted `mofia` object.
#' @param se `"MLR"` (robust sandwich, default), `"ML"` (inverse observed
#'   information) or `"MLF"` (outer product of the scores).
#' @param information `"default"`, `"reduced"` or `"exact"` observed
#'   information; see [mofia()].
#' @param ... Not used.
#' @return A covariance matrix.
#' @export
vcov.mofia <- function(object, se = c("MLR", "ML", "MLF"), information = c("default", "reduced", "exact"), ...) {
  se <- match.arg(se)
  if (is.null(object$vcov)) stop("No standard errors available (fit with se = TRUE).")
  if (se == "MLF") return(object$vcov$MLF)
  object$vcov[[.info_choice(object, information)]][[se]]
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
  cat(sprintf("  N = %d, %d parameters%s, %d^%d adaptive quadrature nodes\n", x$nobs, x$npar,
              if (m$P > nc) paste0(", ", x$link, " link") else "", x$nq, m$D))
  cat(sprintf("  logLik = %.3f, AIC = %.2f, BIC = %.2f\n", x$loglik, stats::AIC(x), stats::BIC(x)))
  cat(sprintf("  %s%s\n", if (isTRUE(x$converged)) "Converged" else "NOT converged",
              if (is.finite(x$newton_decrement)) sprintf(" (Newton decrement %.1e)", x$newton_decrement) else ""))
  invisible(x)
}

#' Summarize a fitted moderated factor model
#'
#' Parameter estimates with standard errors, z statistics and p-values,
#' grouped by parameter type.
#'
#' @param object A fitted `mofia` object.
#' @param se Standard-error type: `"MLR"` (default), `"ML"` or `"MLF"`.
#' @param information `"default"`, `"reduced"` or `"exact"`; see [mofia()].
#' @param ... Not used.
#' @return An object of class `summary.mofia`; its `coefficients` element is a
#'   data frame with columns `estimate`, `se`, `z`, `p` and `type`.
#' @export
summary.mofia <- function(object, se = c("MLR", "ML", "MLF"), information = c("default", "reduced", "exact"), ...) {
  se <- match.arg(se)
  th <- object$coefficients
  s <- if (is.null(object$vcov)) rep(NA_real_, length(th)) else sqrt(pmax(diag(vcov(object, se, information)), 0))
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
  if (!is.null(x$fit$scaling))
    cat(sprintf("Robust scaling correction factor: %.4f\n", x$fit$scaling[[x$information %||% "reduced"]] %||% NA))
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
#' @param information `"default"`, `"reduced"` or `"exact"` information for the
#'   scaling factors; see [mofia()].
#' @return A data frame of class `anova` with one row per model.
#' @details Nesting is checked by parameter names: every parameter of the
#'   smaller model must also be in the larger one. Models that fail this check
#'   are still compared, with a warning, because some models are nested in a
#'   way the names cannot show; the test is only valid for truly nested models.
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
  if (length(unique(vapply(fits, `[[`, 1, "nobs"))) > 1) stop("Models were fitted to different numbers of cases.")
  ord <- order(vapply(fits, `[[`, 1, "npar")); fits <- fits[ord]; nm <- nm[ord]
  n <- length(fits)
  out <- data.frame(npar = vapply(fits, `[[`, 1, "npar"), logLik = vapply(fits, `[[`, 1, "loglik"),
                    AIC = vapply(fits, stats::AIC, 1), BIC = vapply(fits, stats::BIC, 1),
                    Chisq = NA_real_, Df = NA_real_, Scaling = NA_real_, `Pr(>Chisq)` = NA_real_,
                    row.names = nm, check.names = FALSE)
  for (i in seq_len(n)[-1]) {
    tst <- mofia_lrt(fits[[i - 1]], fits[[i]], scaled = scaled, information = information, check_nested = FALSE)
    if (!all(names(fits[[i - 1]]$coefficients) %in% names(fits[[i]]$coefficients)))
      warning(sprintf("'%s' is not nested in '%s' by parameter names; the test is only valid for nested models.",
                      nm[i - 1], nm[i]), call. = FALSE)
    out$Chisq[i] <- tst$statistic; out$Df[i] <- tst$df; out$`Pr(>Chisq)`[i] <- tst$p_value
    if (scaled) out$Scaling[i] <- tst$cd
  }
  if (!scaled) out$Scaling <- NULL
  structure(out, heading = c(if (scaled) "Scaled likelihood-ratio tests (Satorra-Bentler 2001, MLR)\n"
                             else "Likelihood-ratio tests\n"), class = c("anova", "data.frame"))
}
