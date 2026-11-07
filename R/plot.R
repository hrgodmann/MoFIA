# Plots ---------------------------------------------------------------------------------

#' Plot a fitted moderated factor model
#'
#' Three views of how parameters change with a moderator:
#' \describe{
#'   \item{`"factors"`}{factor means and standard deviations across the values
#'     of a moderator, with pointwise confidence bands;}
#'   \item{`"loadings"`}{loadings of items with moderated loadings across a
#'     moderator, with pointwise confidence bands;}
#'   \item{`"items"`}{item response functions: the expected item score given
#'     the factor, at several values of a moderator. Curves that differ show
#'     DIF.}
#' }
#' Continuous moderators are shown over their observed range on the original
#' scale; binary and categorical ones at their values. The other moderators
#' are held at their mean (continuous) or at 0 / the reference level, unless
#' set with `at`. The bands are pointwise. Under variance identification the
#' factor mean and variance of the reference person (all moderators at 0) are
#' fixed, so the bands for the factors narrow to zero width there.
#'
#' @param x A fitted `mofia` object.
#' @param type `"factors"`, `"loadings"` or `"items"`.
#' @param moderator The moderator variable to vary (a variable name, also for
#'   formula moderators). Default: the first.
#' @param items Items to show (`"loadings"`, `"items"`). Default: items with
#'   moderated loadings (`"loadings"`) or with any moderation (`"items"`), at
#'   most six.
#' @param values For `type = "items"`: moderator values at which to draw the
#'   curves. Default: mean - 1 SD, mean, mean + 1 SD for a continuous
#'   moderator, all values otherwise.
#' @param at Named list of values for the other moderator variables.
#' @param level Confidence level of the bands.
#' @param se Standard-error type for the bands; see [vcov.mofia()].
#' @param ... Further graphical parameters passed to [graphics::plot()].
#' @return Invisibly, a data frame with the plotted values.
#' @examples
#' \donttest{
#' data(parabelief)
#' ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession",
#'         "Sorcery", "HauntedPlaces")
#' fit <- mofia(mofia_model(parabelief, items = ub, moderators = c("Age", "Male"),
#'                          threshold_mod = "DevilPossession"))
#' plot(fit, "factors", moderator = "Age")
#' plot(fit, "items", items = "DevilPossession", moderator = "Male")
#' }
#' @export
plot.mofia <- function(x, type = c("factors", "loadings", "items"), moderator = NULL, items = NULL,
                       values = NULL, at = list(), level = 0.95, se = c("MLR", "ML", "MLF"), ...) {
  type <- match.arg(type); se <- match.arg(se)
  model <- x$model; des <- model$design
  if (!length(des$vars)) stop("The model has no moderators.")
  if (is.null(moderator)) moderator <- des$vars[1]
  if (!moderator %in% des$vars) stop("Unknown moderator variable: ", moderator, " (one of: ", paste(des$vars, collapse = ", "), ")")
  inf <- des$info[[moderator]]
  grid <- if (type == "items") {
    if (!is.null(values)) values
    else if (inf$type == "continuous") inf$mean + c(-1, 0, 1) * inf$sd
    else if (inf$type == "binary") inf$values else inf$levels
  } else if (inf$type == "continuous") seq(inf$range[1], inf$range[2], length.out = 60)
    else if (inf$type == "binary") inf$values else inf$levels
  nd <- .reference_data(des, length(grid), at); nd[[moderator]] <- grid
  X <- .design_X(model, nd)
  V <- vcov(x, se); th <- x$coefficients; pt <- model$par
  z <- stats::qnorm(1 - (1 - level) / 2)
  cont <- inf$type == "continuous"
  xs <- if (cont) grid else seq_along(grid)
  xlab_ <- moderator

  # linear combination X[, mods] %*% theta[idx] (+ theta[base]) with delta-method SEs
  lincomb <- function(idx_mod, base = NULL) {
    G <- matrix(0, nrow(X), length(th), dimnames = list(NULL, names(th)))
    for (r in which(!is.na(idx_mod))) G[, idx_mod[r]] <- X[, r]
    if (!is.null(base)) G[, base] <- 1
    est <- drop(G %*% th); s <- sqrt(pmax(rowSums((G %*% V) * G), 0))
    cbind(est = est, lo = est - z * s, hi = est + z * s)
  }
  band <- function(xx, m, col) {
    if (cont) {
      graphics::polygon(c(xx, rev(xx)), c(m[, "lo"], rev(m[, "hi"])), col = grDevices::adjustcolor(col, 0.2), border = NA)
      graphics::lines(xx, m[, "est"], col = col, lwd = 2)
    } else {
      ok <- m[, "hi"] - m[, "lo"] > 1e-10
      if (any(ok)) graphics::arrows(xx[ok], m[ok, "lo"], xx[ok], m[ok, "hi"], angle = 90, code = 3, length = 0.05, col = col)
      graphics::points(xx, m[, "est"], pch = 19, col = col)
    }
  }
  cols <- grDevices::hcl.colors(max(3, model$D, length(items) %||% 6), "Dark 3")
  axis_x <- function() if (!cont) graphics::axis(1, at = xs, labels = grid)

  if (type == "factors") {
    out <- list(); op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op))
    res <- lapply(seq_len(model$D), function(d) {
      im <- match(sprintf("gamma[%s,%s]", model$factors[d], model$moderators), pt$name)
      iv <- match(sprintf("delta[%s,%s]", model$factors[d], model$moderators), pt$name)
      b0 <- if (identical(model$identification, "marker")) match(sprintf("delta0[%s]", model$factors[d]), pt$name) else NULL
      m <- lincomb(im); v <- lincomb(iv, b0)
      list(mean = m, sd = exp(v / 2))
    })
    rng <- function(k) range(unlist(lapply(res, function(r) r[[k]][, c("lo", "hi")])))
    for (k in c("mean", "sd")) {
      graphics::plot(xs, xs, type = "n", ylim = rng(k), xlab = xlab_, xaxt = if (cont) "s" else "n",
                     ylab = if (k == "mean") "Factor mean" else "Factor SD", ...)
      axis_x()
      if (k == "mean") graphics::abline(h = 0, col = "grey80") else graphics::abline(h = 1, col = "grey80")
      for (d in seq_len(model$D)) band(xs + if (cont) 0 else (d - (model$D + 1) / 2) * 0.08, res[[d]][[k]], cols[d])
      if (model$D > 1) graphics::legend("topright", model$factors, col = cols[seq_len(model$D)], lwd = 2, bty = "n")
    }
    out <- do.call(rbind, lapply(seq_len(model$D), function(d)
      data.frame(factor = model$factors[d], value = grid, mean = res[[d]]$mean[, "est"], mean_lo = res[[d]]$mean[, "lo"],
                 mean_hi = res[[d]]$mean[, "hi"], sd = res[[d]]$sd[, "est"], sd_lo = res[[d]]$sd[, "lo"], sd_hi = res[[d]]$sd[, "hi"])))
    return(invisible(out))
  }

  if (type == "loadings") {
    lm <- model$loading_mod
    if (is.null(items)) items <- model$items[apply(lm, 1, any)]
    if (!length(items)) stop("No item has moderated loadings.")
    curves <- list()
    for (it in items) for (d in which(model$loadings[it, ])) {
      im <- match(sprintf("lambda1[%s,%s,%s]", it, model$factors[d], model$moderators), pt$name)
      b0 <- match(sprintf("lambda0[%s,%s]", it, model$factors[d]), pt$name)
      m <- if (is.na(b0)) { r <- lincomb(im); r + 1 } else lincomb(im, b0)     # marker loading fixed to 1
      curves[[sprintf("%s (%s)", it, model$factors[d])]] <- m
    }
    yl <- range(unlist(lapply(curves, function(m) m[, c("lo", "hi")])))
    graphics::plot(xs, xs, type = "n", ylim = yl, xlab = xlab_, ylab = "Loading", xaxt = if (cont) "s" else "n", ...)
    axis_x()
    cols <- grDevices::hcl.colors(max(3, length(curves)), "Dark 3")
    for (k in seq_along(curves)) band(xs + if (cont) 0 else (k - (length(curves) + 1) / 2) * 0.06, curves[[k]], cols[k])
    graphics::legend("topright", names(curves), col = cols[seq_along(curves)], lwd = 2, bty = "n", cex = 0.8)
    out <- do.call(rbind, lapply(names(curves), function(n) data.frame(loading = n, value = grid, curves[[n]])))
    return(invisible(out))
  }

  # item response functions
  if (is.null(items)) {
    modded <- apply(model$loading_mod, 1, any) | apply(model$threshold_mod, 1, any)
    items <- model$items[modded]
    if (!length(items)) items <- model$items
  }
  if (length(items) > 6) { message("Showing the first 6 items; choose others with 'items'."); items <- items[1:6] }
  pr <- .unpack(th, model)
  eta <- seq(-3, 3, length.out = 121)
  nr <- ceiling(length(items) / 3); nc <- min(3, length(items))
  op <- graphics::par(mfrow = c(nr, nc)); on.exit(graphics::par(op))
  cols <- grDevices::hcl.colors(max(3, length(grid)), "Dark 3")
  out <- list()
  for (it in items) {
    j <- match(it, model$items); d <- which(model$loadings[j, ])[1]
    E <- sapply(seq_along(grid), function(g) {
      Xg <- X[g, , drop = FALSE]
      etal <- lapply(seq_len(model$D), function(e) matrix(if (e == d) eta else 0, 1))
      A <- .item_lp(pr, j, Xg, etal, model$loadings)
      if (model$cont[j]) return(drop(A) + pr$nu0[j])
      pk <- .cat_probs(A, pr$tau[[j]], x$link)
      Reduce(`+`, Map(function(p, l) p * l, pk, model$levels[[j]]))[1, ]
    })
    graphics::matplot(eta, E, type = "l", lty = 1, lwd = 2, col = cols[seq_along(grid)],
                      xlab = paste("Factor", model$factors[d]), ylab = "Expected score", main = it, ...)
    lab <- if (cont) format(signif(grid, 3)) else as.character(grid)
    graphics::legend("topleft", paste(moderator, "=", lab), col = cols[seq_along(grid)], lwd = 2, bty = "n", cex = 0.8)
    out[[it]] <- data.frame(item = it, eta = rep(eta, length(grid)), value = rep(grid, each = length(eta)), expected = as.vector(E))
  }
  invisible(do.call(rbind, out))
}

# reference values for the moderator variables (original scale)
.reference_data <- function(des, n, at) {
  nd <- list()
  for (v in des$vars) {
    inf <- des$info[[v]]
    nd[[v]] <- rep(if (!is.null(at[[v]])) at[[v]]
                   else if (inf$type == "continuous") inf$mean
                   else if (inf$type == "binary") inf$values[1]
                   else factor(inf$levels[1], levels = inf$levels), n)
  }
  as.data.frame(nd, stringsAsFactors = FALSE)
}
