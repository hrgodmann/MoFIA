# Score (Lagrange multiplier) tests ---------------------------------------------------

#' Score tests for moderation effects not in the model
#'
#' Tests whether moderation effects that are absent from a fitted model should
#' be added, without refitting: for each candidate parameter (or group of
#' parameters), a score (Lagrange multiplier) test evaluated at the fitted
#' model's estimates. This makes it cheap to screen all items for DIF
#' (threshold/intercept and loading moderation) or to check structural
#' moderation, like modification indices in SEM.
#'
#' @section Statistics:
#' Let \eqn{U_a} be the score of the candidate parameters at the restricted
#' estimates, \eqn{A} the observed information of the extended model and
#' \eqn{B} the outer product of the per-person scores. The model-based statistic
#' is \eqn{\tilde U_a' M^{-1} \tilde U_a} with \eqn{M = A_{aa} - A_{a0}A_{00}^{-1}A_{0a}}
#' and \eqn{\tilde U_a = U_a - A_{a0}A_{00}^{-1}U_0} (\eqn{U_0 \approx 0}). The robust
#' version (Boos, 1992) replaces \eqn{M} by the variance of \eqn{\tilde U_a}
#' estimated from \eqn{B}, which stays valid when the model is misspecified;
#' it is the counterpart of MLR standard errors. Both are asymptotically
#' chi-square with df equal to the number of tested parameters. The expected
#' parameter change (EPC) \eqn{M^{-1}\tilde U_a} approximates the estimate the
#' parameter would get if it were added. The exact observed information is
#' used.
#'
#' @param fit A fitted `mofia` object (the restricted model).
#' @param candidates Character vector of parameter names to test (in the
#'   naming of [mofia_model()], e.g. `"beta[Aliens,Age]"`). Default: all
#'   moderation effects the model does not contain: loading and
#'   threshold/intercept moderation for every item, residual-variance
#'   moderation for continuous items, factor mean and variance moderation, and
#'   correlation moderation (two factors). Marker loadings are never
#'   candidates.
#' @param by `"parameter"` (one 1-df test per candidate) or `"item"` (one test
#'   per item and parameter type, jointly over all moderators, e.g. all
#'   threshold moderation of one item).
#' @param se `"MLR"` (robust, default) or `"ML"` statistic for the p-values.
#' @param adjust Method for adjusted p-values across the tests, passed to
#'   [stats::p.adjust()].
#' @return A data frame sorted by the statistic, with columns `test`,
#'   `parameters`, `df`, `statistic`, `p`, `p_adj`, `statistic_ML` and (for
#'   single parameters) `epc`.
#' @references
#' Boos, D. D. (1992). On generalized score tests. *The American Statistician,
#' 46*(4), 327--333.
#' @examples
#' \donttest{
#' data(parabelief)
#' ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession",
#'         "Sorcery", "HauntedPlaces")
#' scalar <- mofia(mofia_model(parabelief, items = ub, moderators = c("Age", "Male")))
#' mofia_score_test(scalar, by = "item")
#' }
#' @details The fit must have converged (see [mofia()]); otherwise the test is
#'   refused. Each test (parameter or item group) is checked for identification
#'   given the fitted model and reported as `NA` if it is not identified.
#' @export
mofia_score_test <- function(fit, candidates = NULL, by = c("parameter", "item"), se = c("MLR", "ML"),
                             adjust = "holm") {
  by <- match.arg(by); se <- match.arg(se)
  if (!inherits(fit, "mofia")) stop("'fit' must be a fitted mofia object.")
  if (is.null(fit$vcov)) stop("No information matrix available: ", fit$se_note %||% "fit with se = TRUE", ".", call. = FALSE)
  if (!isTRUE(fit$converged))
    stop("The fit has not converged; score tests need the maximum of the restricted model.", call. = FALSE)
  model <- fit$model
  if (is.null(candidates)) candidates <- .all_candidates(model)
  candidates <- unique(candidates)
  if (!length(candidates)) stop("No candidate parameters to test.")
  inmod <- intersect(candidates, model$par$name)
  if (length(inmod)) stop("Already in the model: ", paste(inmod, collapse = ", "))
  ext <- .add_params(model, candidates)
  miss <- setdiff(candidates, ext$par$name)
  if (length(miss)) stop("Not valid candidate parameter(s): ", paste(miss, collapse = ", "))

  th <- stats::setNames(numeric(nrow(ext$par)), ext$par$name)
  th[names(fit$coefficients)] <- fit$coefficients
  nodes <- .mofia_nodes(fit$nq, ext$D)
  sc <- function(t) attr(.mofia_ll(t, ext, fit$link, nodes, scores = TRUE, adapt = fit$adapt), "scores")
  S <- sc(th); U <- colSums(S)
  p0 <- names(fit$coefficients)
  # columns of the Hessian in the candidate directions (central differences of the analytic gradient)
  H <- matrix(0, length(th), length(candidates), dimnames = list(names(th), candidates))
  for (a in candidates) {
    h <- 1e-5; e <- stats::setNames(numeric(length(th)), names(th)); e[a] <- h
    H[, a] <- (colSums(sc(th + e)) - colSums(sc(th - e))) / (2 * h)
  }
  A00inv <- fit$vcov$exact$ML[p0, p0]
  Aa0 <- -t(H[p0, , drop = FALSE]); Aaa <- -H[candidates, , drop = FALSE]; Aaa <- (Aaa + t(Aaa)) / 2
  B <- crossprod(S)

  groups <- if (by == "parameter") as.list(candidates) else split(candidates, .candidate_group(candidates))
  res <- lapply(groups, function(a) {
    K <- Aa0[a, , drop = FALSE] %*% A00inv                    # q x p0
    M <- Aaa[a, a, drop = FALSE] - K %*% t(Aa0[a, , drop = FALSE])
    Ut <- U[a] - drop(K %*% U[p0])
    Dm <- cbind(-K, diag(length(a)))
    V <- Dm %*% B[c(p0, a), c(p0, a)] %*% t(Dm)
    if (min(eigen(M, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8 * max(1, abs(diag(Aaa[a, a, drop = FALSE]))))
      return(c(df = length(a), MLR = NA_real_, ML = NA_real_, epc = NA_real_))   # not identified given the model
    T_ml <- tryCatch(drop(t(Ut) %*% solve(M, Ut)), error = function(e) NA_real_)
    T_r <- tryCatch(drop(t(Ut) %*% solve(V, Ut)), error = function(e) NA_real_)
    epc <- if (length(a) == 1) tryCatch(drop(solve(M, Ut)), error = function(e) NA_real_) else NA_real_
    c(df = length(a), MLR = unname(T_r), ML = unname(T_ml), epc = unname(epc))
  })
  lab <- if (by == "parameter") candidates else names(groups)
  out <- data.frame(test = lab, parameters = vapply(groups, paste, "", collapse = ", "),
                    df = vapply(res, `[[`, 1, "df"), statistic = vapply(res, `[[`, 1, se),
                    statistic_ML = vapply(res, `[[`, 1, "ML"), epc = vapply(res, `[[`, 1, "epc"),
                    stringsAsFactors = FALSE)
  out$p <- stats::pchisq(out$statistic, out$df, lower.tail = FALSE)
  out$p_adj <- stats::p.adjust(out$p, method = adjust)
  out <- out[order(-out$statistic), c("test", "parameters", "df", "statistic", "p", "p_adj", "statistic_ML", "epc")]
  if (by == "item") out$epc <- NULL
  rownames(out) <- NULL
  attr(out, "se") <- se
  out
}

# all moderation effects a model could have but does not
.all_candidates <- function(model) {
  lm <- model$loading_mod; lm[] <- TRUE; lm <- lm & array(model$loadings, dim(lm))
  if (identical(model$identification, "marker"))
    for (d in seq_len(model$D)) lm[model$markers[d], d, ] <- FALSE
  tm <- model$threshold_mod; tm[] <- TRUE
  rm_ <- model$resvar_mod; rm_[model$cont, ] <- TRUE
  mm <- model$mean_mod; mm[] <- TRUE; vm <- model$var_mod; vm[] <- TRUE
  cm <- model$cor_mod; if (model$D == 2 || isTRUE(model$experimental_cor)) cm[] <- TRUE
  full <- .rebuild(model, loading_mod = lm, threshold_mod = tm, mean_mod = mm, var_mod = vm, cor_mod = cm,
                   resvar_mod = rm_, check_identification = FALSE)
  setdiff(full$par$name, model$par$name)
}

# model with the given moderation parameters added
.add_params <- function(model, names) {
  lm <- model$loading_mod; tm <- model$threshold_mod; mm <- model$mean_mod
  vm <- model$var_mod; cm <- model$cor_mod; rm_ <- model$resvar_mod
  for (n in names) {
    b <- sub("\\[.*", "", n); a <- strsplit(sub(".*\\[(.*)\\]", "\\1", n), ",")[[1]]
    ok <- tryCatch({
      if (b == "lambda1") lm[a[1], a[2], a[3]] <- TRUE
      else if (b == "beta") tm[a[1], a[2]] <- TRUE
      else if (b == "gamma") mm[a[1], a[2]] <- TRUE
      else if (b == "delta") vm[a[1], a[2]] <- TRUE
      else if (b %in% c("zcor1", "cov1")) cm[paste(a[1], a[2], sep = ","), a[3]] <- TRUE
      else if (b == "omega") rm_[a[1], a[2]] <- TRUE
      else stop()
      TRUE
    }, error = function(e) FALSE)
    if (!ok) stop("Not a moderation parameter of this model: ", n)
  }
  .rebuild(model, loading_mod = lm, threshold_mod = tm, mean_mod = mm, var_mod = vm, cor_mod = cm, resvar_mod = rm_,
           check_identification = FALSE)
}

# group label: parameter type and item/factor (all moderators together)
.candidate_group <- function(names) {
  b <- sub("\\[.*", "", names); a <- strsplit(sub(".*\\[(.*)\\]", "\\1", names), ",")
  vapply(seq_along(names), function(i) {
    k <- a[[i]]
    switch(b[i], lambda1 = sprintf("lambda1[%s,%s]", k[1], k[2]), zcor1 = , cov1 = sprintf("%s[%s,%s]", b[i], k[1], k[2]),
           sprintf("%s[%s]", b[i], k[1]))
  }, "")
}
