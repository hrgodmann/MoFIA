#' MoFIA: Moderated Factor and Item Analysis
#'
#' Fits moderated nonlinear factor analysis (MNLFA) models by marginal maximum
#' likelihood. Item parameters (loadings, thresholds or intercepts, residual
#' variances) and factor parameters (means, variances, correlations) can depend
#' on any number of continuous or categorical moderators. Items can be binary,
#' ordinal (logit or probit link), continuous, or a mix.
#'
#' The main functions are [mofia_model()] to specify a model, [mofia()] to fit
#' it, [summary.mofia()] and [anova.mofia()] to inspect and compare fits, and
#' [mofia_wald()] for Wald tests of sets of moderation effects.
#'
#' Computation uses adaptive product Gauss-Hermite quadrature with analytic
#' gradients in C++ (multithreaded with OpenMP where the compiler supports it;
#' set `options(mofia.threads = n)`).
#'
#' @references
#' Bauer, D. J., & Hussong, A. M. (2009). Psychometric approaches for
#' developing commensurate measures across independent studies: Traditional and
#' new models. *Psychological Methods, 14*(2), 101--125.
#' \doi{10.1037/a0015583}
#'
#' Curran, P. J., McGinley, J. S., Bauer, D. J., Hussong, A. M., Burns, A.,
#' Chassin, L., Sher, K., & Zucker, R. (2014). A moderated nonlinear factor
#' model for the development of commensurate measures in integrative data
#' analysis. *Multivariate Behavioral Research, 49*(3), 214--231.
#' \doi{10.1080/00273171.2014.889594}
#'
#' Bauer, D. J. (2017). A more general model for testing measurement
#' invariance and differential item functioning. *Psychological Methods,
#' 22*(3), 507--526. \doi{10.1037/met0000077}
#'
#' @keywords internal
#' @useDynLib MoFIA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats vcov coef logLik nobs anova
"_PACKAGE"
