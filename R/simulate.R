#' Simulate data from a moderated factor model
#'
#' Draws item responses from a model at given parameter values, e.g. for power
#' analyses or to check a model's estimation.
#'
#' @param object A `mofia_model` or a fitted `mofia` object (whose estimates and
#'   link are then the defaults).
#' @param theta Named parameter vector (names as in `object$par$name`).
#' @param X Moderator values: a matrix or data frame with the model's moderator
#'   columns. Default: the moderator values of the model's data.
#' @param link `"logit"` or `"probit"` for categorical items.
#' @param seed Optional random seed.
#' @return A data frame with the simulated items (categorical items in the
#'   original coding of the model's data) and the moderators. The simulated
#'   factor scores are attached as attribute `"eta"`.
#' @examples
#' \donttest{
#' data(parabelief)
#' ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession")
#' fit <- mofia(mofia_model(parabelief, items = ub, moderators = "Male"))
#' sim <- mofia_simulate(fit, seed = 1)
#' head(sim)
#' }
#' @export
mofia_simulate <- function(object, theta = NULL, X = NULL, link = NULL, seed = NULL) {
  model <- if (inherits(object, "mofia")) object$model else object
  if (is.null(theta)) {
    if (!inherits(object, "mofia")) stop("Give 'theta'.")
    theta <- object$coefficients
  }
  if (is.null(link)) link <- if (inherits(object, "mofia")) object$link else "logit"
  link <- match.arg(link, c("logit", "probit"))
  miss <- setdiff(model$par$name, names(theta)); if (length(miss)) stop("Missing parameter(s): ", paste(miss, collapse = ", "))
  theta <- theta[model$par$name]
  if (is.null(X)) X <- model$X
  X <- as.matrix(as.data.frame(X)[, model$moderators, drop = FALSE]); storage.mode(X) <- "double"
  if (!is.null(seed)) set.seed(seed)
  sp <- model; sp$X <- X; sp$N <- nrow(X); N <- nrow(X)
  pr <- .unpack(theta, sp); st <- .structural(pr, sp)
  if (!st$ok) stop("The parameters imply an invalid factor correlation for some moderator values.")
  z <- matrix(stats::rnorm(N * sp$D), N, sp$D)
  eta <- matrix(0, N, sp$D, dimnames = list(NULL, sp$factors))
  for (d in seq_len(sp$D)) {
    v <- if (sp$D == 1) z[, 1] else Reduce(`+`, lapply(seq_len(d), function(e) st$Lr[[d]][[e]] * z[, e]))
    eta[, d] <- st$mu[, d] + st$s[, d] * v
  }
  rerr <- if (link == "logit") stats::rlogis else stats::rnorm
  Y <- matrix(NA_real_, N, sp$P, dimnames = list(NULL, sp$items))
  for (j in seq_len(sp$P)) {
    a <- drop(X %*% pr$B[j, ])
    for (d in seq_len(sp$D)) if (sp$loadings[j, d]) a <- a + (pr$L0[j, d] + drop(X %*% pr$L1[j, d, ])) * eta[, d]
    if (isTRUE(sp$cont[j])) {
      Y[, j] <- a + pr$nu0[j] + stats::rnorm(N, 0, sqrt(exp(pr$Om0[j] + drop(X %*% pr$Om[j, ]))))
    } else {
      k <- findInterval(a + rerr(N), pr$tau[[j]], left.open = TRUE)
      Y[, j] <- model$levels[[j]][k + 1]
    }
  }
  out <- data.frame(Y, X, check.names = FALSE)
  attr(out, "eta") <- eta
  out
}
