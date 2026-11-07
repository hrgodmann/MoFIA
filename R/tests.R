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
#' @param information `"default"`, `"reduced"` or `"exact"`; see [mofia()].
#' @param sb2010 Also report the strictly positive Satorra-Bentler (2010)
#'   version, which evaluates the scaling factor of `fit1` at the estimates of
#'   `fit0` (one extra information computation).
#' @param check_nested Stop if the parameter names of `fit0` are not a subset of
#'   those of `fit1`?
#' @return A data frame with the test statistic, degrees of freedom, scaling
#'   factors and p-value.
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
  if (fit0$npar >= fit1$npar) stop("'fit0' must have fewer parameters than 'fit1'.")
  if (check_nested && !all(names(fit0$coefficients) %in% names(fit1$coefficients)))
    stop("'fit0' is not nested in 'fit1' (parameter names); use check_nested = FALSE to override.")
  p0 <- fit0$npar; p1 <- fit1$npar; L0 <- fit0$loglik; L1 <- fit1$loglik
  LR <- -2 * (L0 - L1); df <- p1 - p0
  if (!scaled)
    return(data.frame(test = "Likelihood ratio", LR = LR, df = df, statistic = LR,
                      p_value = stats::pchisq(LR, df, lower.tail = FALSE)))
  info <- .info_choice(fit1, information)
  c0 <- fit0$scaling[[info]]; c1 <- fit1$scaling[[info]]
  if (is.null(c0) || is.null(c1)) stop("Scaling factors are missing (fit with se = TRUE).")
  cd <- (p0 * c0 - p1 * c1) / (p0 - p1)
  out <- data.frame(test = "Satorra-Bentler (2001)", LR = LR, df = df, c0 = c0, c1 = c1, cd = cd,
                    statistic = LR / cd)
  if (sb2010) {
    th10 <- stats::setNames(numeric(p1), names(fit1$coefficients))
    th10[names(fit0$coefficients)] <- fit0$coefficients
    c10 <- .scaling_at(fit1, th10)[[info]]
    cd10 <- (p0 * c0 - p1 * c10) / (p0 - p1)
    out <- rbind(out, data.frame(test = "Satorra-Bentler (2010)", LR = LR, df = df, c0 = c0, c1 = c10,
                                 cd = cd10, statistic = LR / cd10))
  }
  out$p_value <- ifelse(out$cd > 0, stats::pchisq(out$statistic, out$df, lower.tail = FALSE), NA)
  out
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
#' @param information `"default"`, `"reduced"` or `"exact"`; see [mofia()].
#' @param regex Treat `params` as a regular expression?
#' @return A data frame with the Wald statistic, degrees of freedom, p-value and
#'   the tested parameters.
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
  nms <- names(fit$coefficients)
  idx <- if (regex) grep(params, nms) else match(params, nms)
  if (regex && !length(idx)) stop("No parameter matches '", params, "'.")
  if (anyNA(idx)) stop("Unknown parameter(s): ", paste(params[is.na(idx)], collapse = ", "))
  V <- vcov(fit, se, information)
  th <- fit$coefficients[idx]
  W <- drop(t(th) %*% solve(V[idx, idx, drop = FALSE], th))
  data.frame(statistic = W, df = length(idx), p_value = stats::pchisq(W, length(idx), lower.tail = FALSE),
             parameters = paste(nms[idx], collapse = ", "))
}
