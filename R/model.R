# Model specification -----------------------------------------------------------

#' Specify a moderated factor model
#'
#' Sets up a moderated (nonlinear) factor model: which items load on which
#' factors, which moderators act on which item and factor parameters, and how
#' the model is identified. The result is passed to [mofia()] for estimation.
#'
#' @section Model:
#' For person \eqn{i} with moderator values \eqn{x_i} (no intercept column),
#' factors \eqn{\eta_i \mid x_i \sim N(\mu_i, \Sigma_i)} with
#' \itemize{
#'   \item factor means \eqn{\mu_{id} = x_i'\gamma_d} (baseline mean 0);
#'   \item factor variances \eqn{\sigma^2_{id} = \exp(\delta_{0d} + x_i'\delta_d)},
#'     where \eqn{\delta_{0d} = 0} under variance identification;
#'   \item factor correlations \eqn{\rho_i = \tanh(\zeta_0 + x_i'\zeta_1)}
#'     (`cov_param = "fisherz"`) or covariances \eqn{\psi_i = c_0 + x_i'c_1}
#'     (`cov_param = "linear"`).
#' }
#' Item \eqn{j} has loadings \eqn{\lambda_{ijd} = \lambda_{0jd} + x_i'\lambda_{1jd}}
#' and a moderator shift \eqn{\nu_{ij} = x_i'\beta_j}. For a categorical item
#' (binary or ordinal) with \eqn{K} categories,
#' \deqn{P(y_{ij} \ge k \mid \eta_i) = F(\nu_{ij} + \lambda_{ij}'\eta_i - \tau_{jk}),}
#' with \eqn{F} the logistic or standard normal distribution function, so a
#' positive \eqn{\beta} shifts all thresholds of the item down (uniform DIF).
#' A continuous item is normal with mean \eqn{\nu_{0j} + \nu_{ij} + \lambda_{ij}'\eta_i}
#' and residual variance \eqn{\exp(\omega_{0j} + x_i'\omega_j)}.
#'
#' @section Specifying moderation:
#' Item-level moderation (`loading_mod`, `threshold_mod`, `resvar_mod`) can be
#' given as `FALSE`/`NULL` (none), `TRUE` (all items by all moderators), a
#' character vector of item names (these items, by all moderators), a named
#' list mapping items to moderators (e.g. `list(item3 = "age")`), or a logical
#' items x moderators matrix. `loading_mod` also accepts an items x factors x
#' moderators array.
#'
#' Factor-level moderation (`mean_mod`, `var_mod`) can be `TRUE`/`FALSE`, a
#' character vector of moderators (applied to all factors), a named list
#' mapping factors to moderators, or a logical factors x moderators matrix.
#' `cor_mod` can be `TRUE`/`FALSE`, a character vector of moderators (applied
#' to all factor pairs), or a logical pairs x moderators matrix.
#'
#' Matrices and arrays must have exactly the model's dimensions. If they carry
#' dimension names (items, factors, moderators, factor pairs such as
#' `"f1,f2"`), these must match the model's names and are used to align the
#' entries; without names the entries are taken in the model's order. Nothing is
#' recycled.
#'
#' @section Identification:
#' Fixing the baseline factor mean and variance (or a marker loading) does not
#' identify every moderation pattern. If a factor's mean and the
#' thresholds/intercepts of *all* its items depend on the same moderator, a
#' shift of the factor mean cannot be told apart from equal shifts of all
#' items; likewise for the factor variance and *all* loadings. Keep at least
#' one anchor item per factor and moderator whose parameters are not moderated.
#' Patterns that are exactly unidentified are refused with an error; patterns
#' identified only through the functional form of the moderation (for example,
#' exponential variance against linear loadings for a continuous moderator)
#' give a warning. [mofia()] additionally checks that the information matrix at
#' the solution is positive definite. Which items are invariant is an
#' assumption: choosing anchors from the same data (for example with
#' [mofia_score_test()]) and then testing on those data is a form of
#' post-selection inference.
#'
#' Cases without any observed item response are removed with a warning (they
#' do not contribute to the likelihood). Cases with some observed responses
#' contribute the likelihood of their observed items (missing at random).
#'
#' @section Parameter names:
#' `tau[item,k]` thresholds; `nu0[item]` intercepts of continuous items;
#' `lambda0[item,factor]` baseline loadings; `lambda1[item,factor,moderator]`
#' loading moderation; `beta[item,moderator]` threshold/intercept moderation
#' (DIF); `omega0[item]`, `omega[item,moderator]` log residual variances of
#' continuous items; `gamma[factor,moderator]` factor-mean moderation;
#' `delta0[factor]` baseline log factor variance (marker identification only);
#' `delta[factor,moderator]` factor-variance moderation; `zcor0[f1,f2]`,
#' `zcor1[f1,f2,moderator]` Fisher-z factor correlation (or `cov0`, `cov1` for
#' the linear covariance).
#'
#' @param data A data frame containing the items and moderators.
#' @param items Character vector of item names. Can be omitted when `loadings`
#'   is a named list.
#' @param moderators The moderators: a character vector of numeric columns, or
#'   a one-sided formula such as `~ age + gender + age:gender + I(age^2)`. With a
#'   formula, factor and character variables are dummy-coded (treatment
#'   contrasts) and interactions and transformations are allowed; the resulting
#'   model-matrix columns are the moderators (a comma in a column name, as in
#'   `poly(age, 2)1`, becomes a semicolon). Data-dependent transformations such
#'   as `scale()` or `poly()` keep the values learned from the fitted data when
#'   new data are predicted. Cases with missing moderator values are removed
#'   with a warning; moderators must be finite, non-constant and not linearly
#'   dependent.
#' @param loadings Which items load on which factors: a measurement string in
#'   lavaan style (`"f1 =~ y1 + y2 + y3; f2 =~ y4 + y5"`), a named list mapping
#'   factor names to item names (e.g. `list(f1 = c("y1", "y2"), f2 = ...)`), or
#'   a logical items x factors matrix. `NULL` gives one factor measured by all
#'   items.
#' @param factors Factor names; taken from `loadings` when it is a list.
#' @param loading_mod,threshold_mod Moderation of loadings and of
#'   thresholds/intercepts; see *Specifying moderation*. Default: none.
#' @param mean_mod,var_mod Moderation of factor means and variances. Default:
#'   all moderators.
#' @param cor_mod Moderation of the factor correlations (two or more factors).
#'   Default: none.
#' @param cov_param `"fisherz"` (moderated correlation on the Fisher-z scale,
#'   recommended) or `"linear"` (linearly moderated covariance).
#' @param identification `"variance"`: baseline factor variances fixed to 1 and
#'   all loadings free. `"marker"`: one loading per factor fixed to 1 and the
#'   baseline log variance `delta0` estimated.
#' @param markers With `identification = "marker"`: one item per factor whose
#'   loading is fixed to 1 (it must not be moderated). Default: the first item
#'   with an unmoderated loading.
#' @param continuous Names of items treated as continuous (normal). All other
#'   items are categorical (binary or ordinal; any integer coding).
#' @param resvar_mod Moderation of the residual variances of continuous items.
#' @param standardize Standardize continuous moderators (z-scores) before the
#'   model is built? `TRUE` standardizes all continuous moderator variables
#'   (those with more than two distinct values), a character vector only the
#'   named ones. Baseline parameters then refer to a person at the mean, and
#'   moderation effects are per standard deviation. New data for [predict()] and
#'   plots are given on the original scale.
#' @param sparse Warn when a category of a categorical item has fewer than this
#'   many responses (0 = never). Categories never chosen are always reported.
#' @param experimental_cor With three or more factors, moderated correlations
#'   are refused unless `TRUE`: independent pairwise Fisher-z correlations do
#'   not guarantee a positive-definite correlation matrix for all moderator
#'   values, and this case is not methodologically established.
#'
#' @return An object of class `mofia_model`.
#' @seealso [mofia()], [mofia_drop()]
#' @examples
#' data(parabelief)
#' m <- mofia_model(parabelief,
#'   loadings = list(EHA = c("HealingHands", "LifeCoursePossessions", "Graphology",
#'                           "Dowsing", "MindReading", "Willpower",
#'                           "SeeingFuture", "SeeingEvents"),
#'                   UB = c("SpiritsGhosts", "GnomesElves", "Aliens",
#'                          "DevilPossession", "Sorcery", "HauntedPlaces")),
#'   moderators = c("Age", "Male"),
#'   threshold_mod = c("Dowsing", "DevilPossession"), cor_mod = TRUE)
#' m
#'
#' # the same, with measurement syntax and a moderator formula
#' m2 <- mofia_model(parabelief,
#'   loadings = "EHA =~ HealingHands + LifeCoursePossessions + Graphology + Dowsing +
#'                      MindReading + Willpower + SeeingFuture + SeeingEvents
#'               UB =~ SpiritsGhosts + GnomesElves + Aliens + DevilPossession +
#'                     Sorcery + HauntedPlaces",
#'   moderators = ~ Age * Male, threshold_mod = c("Dowsing", "DevilPossession"),
#'   cor_mod = TRUE)
#' m2$moderators
#'
#' # not identified: the factor means and all thresholds depend on Age and Male
#' try(mofia_model(parabelief, loadings = list(UB = c("SpiritsGhosts", "Aliens",
#'                 "Sorcery")), moderators = c("Age", "Male"), threshold_mod = TRUE))
#' @export
mofia_model <- function(data, items = NULL, moderators = NULL, loadings = NULL, factors = NULL,
                        loading_mod = NULL, threshold_mod = NULL,
                        mean_mod = TRUE, var_mod = TRUE, cor_mod = FALSE,
                        cov_param = c("fisherz", "linear"),
                        identification = c("variance", "marker"), markers = NULL,
                        continuous = NULL, resvar_mod = NULL, standardize = FALSE,
                        sparse = 5, experimental_cor = FALSE) {
  cov_param <- match.arg(cov_param); identification <- match.arg(identification)
  if (!is.data.frame(data)) data <- as.data.frame(data)

  # ---- moderators: columns or a formula; missing values; standardization ----
  des <- .moderator_design(data, moderators, standardize)
  data <- des$data; moderators <- des$moderators; moddata <- des$raw

  # ---- items, factors and the loading pattern ----
  if (is.character(loadings) && length(loadings) == 1 && grepl("=~", loadings, fixed = TRUE))
    loadings <- .parse_loadings(loadings)
  if (is.list(loadings) && !is.data.frame(loadings)) {
    if (is.null(names(loadings)) || any(names(loadings) == "")) stop("'loadings' must be a named list (factor = items).")
    if (anyDuplicated(names(loadings))) stop("Duplicated factor names in 'loadings'.")
    if (is.null(factors)) factors <- names(loadings)
    bad <- setdiff(names(loadings), factors)
    if (length(bad) || length(factors) != length(names(loadings)))
      stop("'factors' must contain exactly the factor names of 'loadings'.")
    if (is.null(items)) items <- unique(unlist(loadings, use.names = FALSE))
    lp <- matrix(FALSE, length(items), length(factors), dimnames = list(items, factors))
    for (f in names(loadings)) {
      bad <- setdiff(loadings[[f]], items); if (length(bad)) stop("Unknown item(s) in loadings$", f, ": ", paste(bad, collapse = ", "))
      lp[loadings[[f]], f] <- TRUE
    }
    loadings <- lp
  }
  if (is.null(items) && is.matrix(loadings) && !is.null(rownames(loadings))) items <- rownames(loadings)
  if (is.null(items)) stop("Give 'items', or 'loadings' as a named list.")
  items <- as.character(items)
  if (is.null(factors)) factors <- if (is.null(loadings)) "eta" else (colnames(loadings) %||% paste0("f", seq_len(ncol(as.matrix(loadings)))))
  P <- length(items); Q <- length(moderators); D <- length(factors)
  miss <- setdiff(c(items, moderators), names(data))
  if (length(miss)) stop("Not in the data: ", paste(miss, collapse = ", "))
  if (anyDuplicated(items)) stop("Duplicated item names.")
  if (anyDuplicated(factors)) stop("Duplicated factor names.")
  if (is.null(loadings)) {
    if (D > 1) stop("Give the loading pattern ('loadings') for more than one factor.")
    loadings <- matrix(TRUE, P, 1)
  }
  loadings <- .named_logical(loadings, list(items, factors), c("items", "factors"), "loadings")
  if (any(rowSums(loadings) == 0)) stop("Item(s) without a loading: ", paste(items[rowSums(loadings) == 0], collapse = ", "))

  if (any(moderators %in% items)) stop("Items cannot also be moderators: ", paste(intersect(moderators, items), collapse = ", "))

  # ---- cases without any observed item response carry no information ----
  empty <- rowSums(!is.na(as.matrix(data[, items, drop = FALSE]))) == 0
  if (any(empty)) {
    warning(sum(empty), " case(s) without any observed item response removed (they carry no information).", call. = FALSE)
    data <- data[!empty, , drop = FALSE]; moddata <- moddata[!empty, , drop = FALSE]
  }
  if (!nrow(data)) stop("No cases with observed item responses.")

  # ---- moderation patterns ----
  lm <- array(FALSE, c(P, D, Q), dimnames = list(items, factors, moderators))
  if (length(dim(loading_mod)) == 3) {
    lm[] <- .named_logical(loading_mod, list(items, factors, moderators), c("items", "factors", "moderators"),
                           "loading_mod")
  } else {
    m2 <- .item_mod(loading_mod, items, moderators, "loading_mod")
    for (d in seq_len(D)) lm[, d, ] <- m2 & loadings[, d]
  }
  lm <- lm & array(loadings, c(P, D, Q))
  tm <- .item_mod(threshold_mod, items, moderators, "threshold_mod")
  mean_mod <- .factor_mod(mean_mod, factors, moderators, "mean_mod")
  var_mod <- .factor_mod(var_mod, factors, moderators, "var_mod")

  # identification
  mk_idx <- NULL
  if (identification == "marker") {
    mk_idx <- integer(D)
    for (d in seq_len(D)) {
      if (!is.null(markers)) {
        mk_idx[d] <- match(markers[[d]], items)
        if (is.na(mk_idx[d]) || !loadings[mk_idx[d], d]) stop("Marker ", markers[[d]], " does not load on ", factors[d], ".")
      } else {
        cand <- which(loadings[, d] & !apply(lm[, d, , drop = FALSE], 1, any))
        if (!length(cand)) stop("No item with an unmoderated loading on ", factors[d], " to use as marker.")
        mk_idx[d] <- cand[1]
      }
      if (any(lm[mk_idx[d], d, ])) stop("The marker loading (", items[mk_idx[d]], " on ", factors[d], ") must not be moderated.")
    }
    names(mk_idx) <- factors
  }

  # factor pairs (1,2), (1,3), ..., (D-1,D)
  pairs <- if (D >= 2) t(utils::combn(D, 2)) else matrix(integer(0), 0, 2)
  npair <- nrow(pairs)
  cmat <- matrix(FALSE, npair, Q, dimnames = list(if (npair) paste(factors[pairs[, 1]], factors[pairs[, 2]], sep = ",") else NULL, moderators))
  if (npair && Q && !is.null(cor_mod) && !identical(cor_mod, FALSE)) {
    if (isTRUE(cor_mod)) cmat[] <- TRUE
    else if (is.character(cor_mod)) {
      bad <- setdiff(cor_mod, moderators); if (length(bad)) stop("Unknown moderator(s) in cor_mod: ", paste(bad, collapse = ", "))
      cmat[, cor_mod] <- TRUE
    } else if (is.null(dim(cor_mod))) {
      cmat[] <- matrix(.per_moderator(cor_mod, moderators, "cor_mod"), npair, Q, byrow = TRUE)
    } else cmat[] <- .named_logical(cor_mod, dimnames(cmat), c("factor pairs", "moderators"), "cor_mod")
  } else if (!npair && is.character(cor_mod) && length(cor_mod))
    stop("'cor_mod' needs at least two factors.")
  cor_mod <- cmat
  if (D >= 3 && any(cor_mod) && !isTRUE(experimental_cor))
    stop("Moderated factor correlations with 3 or more factors are not methodologically established ",
         "(pairwise Fisher-z correlations need not form a positive-definite matrix). Use unmoderated ",
         "correlations (cor_mod = FALSE), or set experimental_cor = TRUE to fit them anyway.")

  # ---- items: continuous or categorical ----
  bad <- setdiff(continuous, items); if (length(bad)) stop("Unknown continuous item(s): ", paste(bad, collapse = ", "))
  cont <- stats::setNames(items %in% (if (is.null(continuous)) character(0) else continuous), items)
  rm_ <- .item_mod(resvar_mod, items, moderators, "resvar_mod")
  if (any(rm_[!cont, ])) stop("Residual-variance moderation is only defined for continuous items.")
  Y <- as.matrix(data[, items, drop = FALSE])
  if (!is.numeric(Y) && !is.logical(Y)) stop("Items must be numeric (integer-coded categories, or continuous).")
  storage.mode(Y) <- "double"
  K <- integer(P); levels_list <- vector("list", P); gaps <- sparse_items <- character(0)
  for (j in seq_len(P)) {
    if (cont[j]) { K[j] <- 1L; next }
    yj <- Y[!is.na(Y[, j]), j]
    if (any(yj != round(yj))) stop("Categorical item ", items[j], " must be integer-coded (or declare it continuous).")
    lv <- sort(unique(yj)); levels_list[[j]] <- lv
    Y[, j] <- match(Y[, j], lv) - 1L; K[j] <- length(lv)
    if (K[j] < 2) stop("Item ", items[j], " has fewer than 2 observed categories.")
    if (length(lv) > 1 && any(diff(lv) > 1))
      gaps <- c(gaps, sprintf("%s (not observed: %s)", items[j],
                              paste(setdiff(seq(min(lv), max(lv)), lv), collapse = ", ")))
    cnt <- tabulate(Y[!is.na(Y[, j]), j] + 1L, K[j])
    if (sparse > 0 && any(cnt < sparse))
      sparse_items <- c(sparse_items, sprintf("%s (%s)", items[j],
                        paste(sprintf("%s: %d", lv[cnt < sparse], cnt[cnt < sparse]), collapse = ", ")))
  }
  if (length(gaps))
    warning("Categories never chosen were dropped and the remaining ones recoded consecutively: ",
            paste(gaps, collapse = "; "), call. = FALSE)
  if (length(sparse_items))
    warning("Categories with fewer than ", sparse, " responses: ", paste(sparse_items, collapse = "; "),
            ". Estimates of their thresholds can be unstable; see mofia_collapse().", call. = FALSE)
  X <- matrix(as.numeric(as.matrix(data[, moderators, drop = FALSE])), nrow(data), Q, dimnames = list(NULL, moderators))
  .check_moderators(X)
  .check_identification(loadings, lm, tm, mean_mod, var_mod, X, identification, cov_param)

  # ---- parameter table ----
  rows <- list()
  add <- function(n, b, j = NA, d = NA, q = NA, k = NA, pr = NA)
    rows[[length(rows) + 1]] <<- list(name = n, block = b, item = j, factor = d, mod = q, k = k, pair = pr)
  for (j in seq_len(P)) {
    if (cont[j]) add(sprintf("nu0[%s]", items[j]), "nu0", j)
    else for (k in seq_len(K[j] - 1)) add(sprintf("tau[%s,%d]", items[j], k), "tau", j, k = k)
  }
  for (j in seq_len(P)) for (d in seq_len(D)) if (loadings[j, d] && !(identification == "marker" && mk_idx[[d]] == j))
    add(sprintf("lambda0[%s,%s]", items[j], factors[d]), "lambda0", j, d)
  for (j in seq_len(P)) for (d in seq_len(D)) for (q in seq_len(Q)) if (lm[j, d, q])
    add(sprintf("lambda1[%s,%s,%s]", items[j], factors[d], moderators[q]), "lambda1", j, d, q)
  for (j in seq_len(P)) for (q in seq_len(Q)) if (tm[j, q])
    add(sprintf("beta[%s,%s]", items[j], moderators[q]), "beta", j, q = q)
  for (j in which(cont)) {
    add(sprintf("omega0[%s]", items[j]), "omega0", j)
    for (q in which(rm_[j, ])) add(sprintf("omega[%s,%s]", items[j], moderators[q]), "omega", j, q = q)
  }
  for (d in seq_len(D)) for (q in seq_len(Q)) if (mean_mod[d, q])
    add(sprintf("gamma[%s,%s]", factors[d], moderators[q]), "gamma", d = d, q = q)
  for (d in seq_len(D)) for (q in seq_len(Q)) if (var_mod[d, q])
    add(sprintf("delta[%s,%s]", factors[d], moderators[q]), "delta", d = d, q = q)
  if (identification == "marker") for (d in seq_len(D)) add(sprintf("delta0[%s]", factors[d]), "delta0", d = d)
  pre <- if (cov_param == "fisherz") "zcor" else "cov"
  for (p in seq_len(npair)) {
    f1 <- factors[pairs[p, 1]]; f2 <- factors[pairs[p, 2]]
    add(sprintf("%s0[%s,%s]", pre, f1, f2), "r0", pr = p)
    for (q in which(cor_mod[p, ])) add(sprintf("%s1[%s,%s,%s]", pre, f1, f2, moderators[q]), "r1", q = q, pr = p)
  }
  pt <- data.frame(name = vapply(rows, `[[`, "", "name"), block = vapply(rows, `[[`, "", "block"),
                   item = vapply(rows, function(r) as.integer(r$item), 1L),
                   factor = vapply(rows, function(r) as.integer(r$factor), 1L),
                   mod = vapply(rows, function(r) as.integer(r$mod), 1L),
                   k = vapply(rows, function(r) as.integer(r$k), 1L),
                   pair = vapply(rows, function(r) as.integer(r$pair), 1L), stringsAsFactors = FALSE)
  structure(list(Y = Y, X = X, K = K, levels = levels_list, items = items, moderators = moderators,
                 factors = factors, P = P, Q = Q, D = D, N = nrow(Y), loadings = loadings,
                 loading_mod = lm, threshold_mod = tm, mean_mod = mean_mod, var_mod = var_mod,
                 cor_mod = cor_mod, pairs = pairs, npair = npair, cov_param = cov_param,
                 identification = identification, markers = mk_idx, cont = cont, resvar_mod = rm_,
                 experimental_cor = isTRUE(experimental_cor), design = des$design, moddata = moddata, par = pt),
            class = "mofia_model")
}

# ---- input validation ------------------------------------------------------------
# A logical array with the given dimension names. Named inputs must carry exactly
# these names on every axis (in any order) and are reordered; unnamed inputs must
# have exactly the right dimensions. Nothing is recycled.
.named_logical <- function(x, dn, axes, what) {
  x <- if (is.null(dim(x)) && length(dn) == 2) as.matrix(x) else x
  if (is.data.frame(x)) x <- as.matrix(x)
  d <- dim(x); target <- lengths(dn)
  if (length(d) != length(target) || any(d != target))
    stop(sprintf("'%s' must be a %s array/matrix (%s); got dimensions %s.", what, paste(axes, collapse = " x "),
                 paste(target, collapse = " x "), paste(d, collapse = " x ")), call. = FALSE)
  xdn <- dimnames(x)
  if (!is.null(xdn)) for (k in seq_along(dn)) {
    nk <- xdn[[k]]
    if (is.null(nk)) next
    if (anyDuplicated(nk) || !setequal(nk, dn[[k]]) || length(nk) != length(dn[[k]]))
      stop(sprintf("The %s names of '%s' do not match the model's %s: %s.", sub("s$", "", axes[k]), what, axes[k],
                   paste(dn[[k]], collapse = ", ")), call. = FALSE)
  }
  if (!is.null(xdn)) {
    idx <- lapply(seq_along(dn), function(k) if (is.null(xdn[[k]])) seq_along(dn[[k]]) else match(dn[[k]], xdn[[k]]))
    x <- do.call(`[`, c(list(x), idx, list(drop = FALSE)))
  }
  v <- as.vector(x)
  if (is.numeric(v) && all(v %in% c(0, 1))) v <- v == 1
  if (!is.logical(v) || anyNA(v)) stop("'", what, "' must contain only TRUE/FALSE (or 0/1).", call. = FALSE)
  array(v, target, dimnames = dn)
}

# a per-moderator logical vector: length 1, or one value per moderator (by name if named)
.per_moderator <- function(x, moderators, what) {
  if (!is.logical(x) || anyNA(x)) stop("'", what, "' must be TRUE/FALSE, moderator names, a named list or a matrix.", call. = FALSE)
  if (length(x) == 1) return(rep(x, length(moderators)))
  if (length(x) != length(moderators))
    stop("A logical '", what, "' vector needs one value per moderator (", length(moderators), ").", call. = FALSE)
  if (!is.null(names(x))) {
    if (!setequal(names(x), moderators)) stop("The names of '", what, "' do not match the moderators.", call. = FALSE)
    x <- x[moderators]
  }
  unname(x)
}

.check_moderators <- function(X) {
  if (!ncol(X)) return(invisible())
  nf <- colnames(X)[colSums(!is.finite(X)) > 0]
  if (length(nf)) stop("Moderators must be finite numbers; not finite in: ", paste(nf, collapse = ", "), call. = FALSE)
  const <- colnames(X)[apply(X, 2, function(v) max(v) - min(v) == 0)]
  if (length(const)) stop("Constant moderator(s) (no variation in the data): ", paste(const, collapse = ", "), call. = FALSE)
  q <- qr(cbind(1, X))
  if (q$rank < ncol(X) + 1) {
    dep <- colnames(X)[q$pivot[(q$rank + 1):(ncol(X) + 1)] - 1]
    stop("The moderators are linearly dependent (together with the intercept); remove: ",
         paste(dep[dep != ""], collapse = ", "), call. = FALSE)
  }
  invisible()
}

# Identification of the moderation effects. Fixing the baseline factor mean and
# variance (or a marker loading) does not identify every moderation pattern:
#  * mean / DIF: if the factor mean and the thresholds/intercepts of all items of
#    the factor depend on the same moderator, a shift of the factor mean equals
#    equal shifts of all items (exactly so when no loading of the factor is moderated);
#  * scale: if the factor variance and all loadings of the factor depend on the same
#    moderator, a larger factor SD equals smaller loadings (exactly so for a binary
#    moderator when nothing else on the factor is moderated; for a continuous moderator
#    the model is identified only through the different functional forms, exp() vs
#    linear).
# Exactly unidentified patterns are errors, functional-form-only ones warnings. At
# least one anchor item per factor and moderator avoids both.
.mofia_state <- new.env(parent = emptyenv())   # internal switches (not user options)

.check_identification <- function(loadings, lm, tm, mean_mod, var_mod, X, identification, cov_param) {
  D <- ncol(loadings); Q <- ncol(X); if (!Q) return(invisible())
  if (isFALSE(.mofia_state$check_identification)) return(invisible())
  items <- rownames(loadings); factors <- colnames(loadings); mods <- colnames(X)
  binary <- apply(X, 2, function(v) length(unique(v)) == 2)
  errs <- warns <- character(0)
  for (d in seq_len(D)) {
    on_d <- which(loadings[, d])
    any_lm_d <- any(lm[on_d, d, ])
    for (q in seq_len(Q)) {
      if (mean_mod[d, q] && all(tm[on_d, q])) {
        msg <- sprintf("the mean of factor '%s' and the thresholds/intercepts of all its items depend on '%s'", factors[d], mods[q])
        if (!any_lm_d) errs <- c(errs, msg) else warns <- c(warns, msg)
      }
      if (var_mod[d, q] && all(lm[on_d, d, q])) {
        msg <- sprintf("the variance of factor '%s' and the loadings of all its items depend on '%s'", factors[d], mods[q])
        exact <- binary[q] && !any(lm[on_d, d, -q]) && !any(mean_mod[d, -q]) && (D == 1 || cov_param == "fisherz")
        if (exact) errs <- c(errs, msg) else warns <- c(warns, msg)
      }
    }
  }
  if (length(errs))
    stop("The model is not identified: ", paste(errs, collapse = "; "), ". Keep at least one anchor item per ",
         "factor whose parameters do not depend on that moderator, or drop the factor-level moderation ",
         "(mean_mod / var_mod).", call. = FALSE)
  if (length(warns))
    warning("Possibly not identified: ", paste(warns, collapse = "; "), ". Without an anchor item the model is ",
            "identified only through the functional form of the moderation (if at all); estimates can be unstable ",
            "and depend on that assumption. Consider keeping an anchor item.", call. = FALSE)
  invisible()
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# items x moderators logical matrix from the flexible input forms
.item_mod <- function(x, items, moderators, what) {
  m <- matrix(FALSE, length(items), length(moderators), dimnames = list(items, moderators))
  if (is.null(x) || identical(x, FALSE)) return(m)
  if (isTRUE(x)) { m[] <- TRUE; return(m) }
  if (is.character(x) && is.null(dim(x))) {
    bad <- setdiff(x, items); if (length(bad)) stop("Unknown item(s) in ", what, ": ", paste(bad, collapse = ", "))
    m[x, ] <- TRUE; return(m)
  }
  if (is.list(x) && !is.data.frame(x)) {
    if (length(x) && (is.null(names(x)) || any(names(x) == "")))
      stop("'", what, "' given as a list must be named (item = moderators).", call. = FALSE)
    for (r in names(x)) {
      if (!r %in% items) stop("Unknown item in ", what, ": ", r)
      if (is.null(x[[r]]) || identical(x[[r]], FALSE)) next
      mods <- if (isTRUE(x[[r]])) moderators else x[[r]]
      if (!is.character(mods)) stop("Entries of '", what, "' must be moderator names or TRUE: ", r, call. = FALSE)
      bad <- setdiff(mods, moderators); if (length(bad)) stop("Unknown moderator(s) in ", what, ": ", paste(bad, collapse = ", "))
      m[r, mods] <- TRUE
    }
    return(m)
  }
  if (is.logical(x) && is.null(dim(x)) && length(x) == 1) { m[] <- x; return(m) }
  .named_logical(x, dimnames(m), c("items", "moderators"), what)
}

# factors x moderators logical matrix
.factor_mod <- function(x, factors, moderators, what) {
  m <- matrix(FALSE, length(factors), length(moderators), dimnames = list(factors, moderators))
  if (is.null(x) || identical(x, FALSE) || !length(moderators)) return(m)
  if (isTRUE(x)) { m[] <- TRUE; return(m) }
  if (is.character(x) && is.null(dim(x))) {
    bad <- setdiff(x, moderators); if (length(bad)) stop("Unknown moderator(s) in ", what, ": ", paste(bad, collapse = ", "))
    m[, x] <- TRUE; return(m)
  }
  if (is.list(x) && !is.data.frame(x)) {
    if (length(x) && (is.null(names(x)) || any(names(x) == "")))
      stop("'", what, "' given as a list must be named (factor = moderators).", call. = FALSE)
    for (f in names(x)) {
      if (!f %in% factors) stop("Unknown factor in ", what, ": ", f)
      if (is.null(x[[f]]) || identical(x[[f]], FALSE)) next
      mods <- if (isTRUE(x[[f]])) moderators else x[[f]]
      if (!is.character(mods)) stop("Entries of '", what, "' must be moderator names or TRUE: ", f, call. = FALSE)
      bad <- setdiff(mods, moderators); if (length(bad)) stop("Unknown moderator(s) in ", what, ": ", paste(bad, collapse = ", "))
      m[f, mods] <- TRUE
    }
    return(m)
  }
  if (is.null(dim(x))) {
    m[] <- matrix(.per_moderator(x, moderators, what), length(factors), length(moderators), byrow = TRUE)
    return(m)
  }
  .named_logical(x, dimnames(m), c("factors", "moderators"), what)
}

# ---- moderators: columns or a one-sided formula -------------------------------------
# Returns the data with (i) cases with missing moderator variables removed, (ii) the
# chosen variables standardized, and (iii) for a formula, the model-matrix columns
# (dummies, interactions, transformations) added; plus the design needed to rebuild
# the moderator matrix for new data.
.moderator_design <- function(data, moderators, standardize) {
  if (is.null(moderators)) moderators <- character(0)
  is_f <- inherits(moderators, "formula")
  if (is_f && length(moderators) != 2) stop("The moderator formula must be one-sided, e.g. ~ age + male.")
  vars <- if (is_f) all.vars(moderators) else as.character(moderators)
  miss <- setdiff(vars, names(data))
  if (length(miss)) stop("Not in the data: ", paste(miss, collapse = ", "))
  if (length(vars)) {
    nax <- !stats::complete.cases(data[, vars, drop = FALSE])
    if (any(nax)) {
      warning(sum(nax), " case(s) with missing moderator values removed.", call. = FALSE)
      data <- data[!nax, , drop = FALSE]
    }
    nf <- vars[vapply(vars, function(v) is.numeric(data[[v]]) && any(!is.finite(data[[v]])), TRUE)]
    if (length(nf)) stop("Moderators must be finite numbers; not finite in: ", paste(nf, collapse = ", "), call. = FALSE)
  }
  raw <- data[, vars, drop = FALSE]                    # moderator variables on their original scale
  rownames(raw) <- NULL
  # description of each raw variable (on its original scale), used for plots and new data
  info <- lapply(vars, function(v) {
    x <- data[[v]]
    if (is.numeric(x) || is.logical(x)) {
      ux <- sort(unique(as.numeric(x)))
      if (length(ux) <= 2) list(type = "binary", values = ux)
      else list(type = "continuous", mean = mean(x), sd = stats::sd(x), range = range(x))
    } else list(type = "factor", levels = levels(factor(x)))
  })
  names(info) <- vars
  std <- if (isTRUE(standardize)) vars[vapply(info, function(i) i$type == "continuous", TRUE)]
         else if (is.character(standardize)) standardize else character(0)
  bad <- setdiff(std, vars[vapply(info, function(i) i$type == "continuous", TRUE)])
  if (length(bad)) stop("Only continuous moderators can be standardized: ", paste(bad, collapse = ", "))
  scaling <- list()
  for (v in std) {
    scaling[[v]] <- c(mean = info[[v]]$mean, sd = info[[v]]$sd)
    data[[v]] <- (data[[v]] - info[[v]]$mean) / info[[v]]$sd
  }
  if (!is_f) {
    for (q in vars) if (!is.numeric(data[[q]]) && !is.logical(data[[q]]))
      stop("Moderator '", q, "' is not numeric. Give the moderators as a formula (e.g. ~ ", q,
           ") to dummy-code it automatically.")
    for (q in vars) data[[q]] <- as.numeric(data[[q]])
    return(list(data = data, moderators = vars, raw = raw,
                design = list(type = "columns", vars = vars, info = info, scaling = scaling)))
  }
  tt <- stats::terms(moderators)
  if (attr(tt, "intercept") == 0) stop("The moderator formula must keep the intercept (do not use '- 1' or '+ 0').")
  if (!length(attr(tt, "term.labels"))) stop("The moderator formula has no terms.")
  mf <- stats::model.frame(tt, data, na.action = stats::na.pass)
  # the model frame's terms carry 'predvars': data-dependent transformations such as
  # scale() or poly() are re-applied to new data with the training values
  tt <- stats::terms(mf)
  mm <- stats::model.matrix(tt, mf)
  xlev <- stats::.getXlevels(tt, mf); contr <- attr(mm, "contrasts")
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  # commas separate the indices in parameter names; terms such as poly(x, 2) get ';'
  mods <- gsub(",[[:space:]]*", "; ", colnames(mm))
  for (k in seq_along(mods)) data[[mods[k]]] <- unname(mm[, k])
  list(data = data, moderators = mods, raw = raw,
       design = list(type = "formula", formula = moderators, terms = tt, xlevels = xlev, contrasts = contr,
                     vars = vars, info = info, scaling = scaling))
}

# moderator matrix for new data (raw variables on their original scale)
.design_X <- function(model, newdata) {
  d <- model$design
  newdata <- as.data.frame(newdata)
  miss <- setdiff(d$vars, names(newdata))
  if (length(miss)) stop("'newdata' lacks the moderator variable(s): ", paste(miss, collapse = ", "))
  for (v in names(d$scaling)) newdata[[v]] <- (newdata[[v]] - d$scaling[[v]][["mean"]]) / d$scaling[[v]][["sd"]]
  if (d$type == "columns") {
    X <- as.matrix(newdata[, model$moderators, drop = FALSE])
  } else {
    tt <- stats::delete.response(d$terms)
    mf <- stats::model.frame(tt, newdata, xlev = d$xlevels, na.action = stats::na.pass)
    X <- stats::model.matrix(tt, mf, contrasts.arg = d$contrasts)
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }
  X <- matrix(as.numeric(X), nrow(newdata), length(model$moderators), dimnames = list(NULL, model$moderators))
  if (anyNA(X)) stop("Missing moderator values in 'newdata'.")
  X
}

# lavaan-style measurement syntax: "f1 =~ y1 + y2 + y3; f2 =~ y4 + y5" (statements may
# span several lines; a new statement starts at "<name> =~")
.parse_loadings <- function(txt) {
  txt <- gsub("#[^\n]*", "", txt)
  txt <- gsub("[\r\n\t]", " ", txt)
  txt <- gsub("([A-Za-z.][A-Za-z0-9._]*[[:space:]]*=~)", ";\\1", txt, perl = TRUE)
  stm <- trimws(unlist(strsplit(txt, ";", fixed = TRUE)))
  stm <- stm[nzchar(stm)]
  if (!length(stm)) stop("No measurement statements ('factor =~ items') found.")
  out <- list()
  for (l in stm) {
    parts <- strsplit(l, "=~", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("Cannot read the measurement statement: ", l)
    f <- trimws(parts[1]); its <- trimws(strsplit(parts[2], "+", fixed = TRUE)[[1]]); its <- its[nzchar(its)]
    if (!nzchar(f) || !length(its)) stop("Cannot read the measurement statement: ", l)
    out[[f]] <- unique(c(out[[f]], its))
  }
  out
}

#' Remove moderation effects from a model
#'
#' Returns the model with selected moderation parameters removed, e.g. to build
#' the restricted model for a likelihood-ratio test.
#'
#' @param model A `mofia_model` (or a fitted `mofia` object, whose model is used).
#' @param drop Character vector of parameter names to remove. Only moderation
#'   effects can be removed: `lambda1`, `beta`, `gamma`, `delta`, `zcor1`/`cov1`
#'   and `omega` parameters.
#' @return A `mofia_model`.
#' @examples
#' data(parabelief)
#' m1 <- mofia_model(parabelief, items = c("SpiritsGhosts", "GnomesElves", "Aliens",
#'                   "DevilPossession", "Sorcery", "HauntedPlaces"),
#'                   moderators = c("Age", "Male"), threshold_mod = "DevilPossession")
#' m0 <- mofia_drop(m1, c("beta[DevilPossession,Age]", "beta[DevilPossession,Male]"))
#' @export
mofia_drop <- function(model, drop) {
  if (inherits(model, "mofia")) model <- model$model
  lm <- model$loading_mod; tm <- model$threshold_mod; mm <- model$mean_mod
  vm <- model$var_mod; cm <- model$cor_mod; rm_ <- model$resvar_mod
  for (n in drop) {
    if (!n %in% model$par$name) stop("Not a parameter of this model: ", n)
    b <- sub("\\[.*", "", n); a <- strsplit(sub(".*\\[(.*)\\]", "\\1", n), ",")[[1]]
    if (b == "lambda1") lm[a[1], a[2], a[3]] <- FALSE
    else if (b == "beta") tm[a[1], a[2]] <- FALSE
    else if (b == "gamma") mm[a[1], a[2]] <- FALSE
    else if (b == "delta") vm[a[1], a[2]] <- FALSE
    else if (b %in% c("zcor1", "cov1")) { pr <- model$par$pair[model$par$name == n]; cm[pr, a[3]] <- FALSE }
    else if (b == "omega") rm_[a[1], a[2]] <- FALSE
    else stop("Only moderation effects can be dropped (lambda1, beta, gamma, delta, zcor1/cov1, omega): ", n)
  }
  .rebuild(model, loading_mod = lm, threshold_mod = tm, mean_mod = mm, var_mod = vm, cor_mod = cm, resvar_mod = rm_)
}

# rebuild a model from its own data with changed moderation patterns. With
# check_identification = FALSE the model is only used to evaluate scores (e.g. the
# union of all score-test candidates, which need not be identified jointly).
.rebuild <- function(model, ..., check_identification = TRUE) {
  if (!check_identification) {
    old <- .mofia_state$check_identification; .mofia_state$check_identification <- FALSE
    on.exit(.mofia_state$check_identification <- old)
  }
  Y <- model$Y
  for (j in seq_len(model$P)) if (!model$cont[j]) Y[, j] <- model$levels[[j]][Y[, j] + 1]
  dat <- data.frame(Y, model$X, check.names = FALSE)
  args <- list(data = dat, items = model$items, moderators = model$moderators, loadings = model$loadings,
               factors = model$factors, loading_mod = model$loading_mod, threshold_mod = model$threshold_mod,
               mean_mod = model$mean_mod, var_mod = model$var_mod, cor_mod = model$cor_mod,
               cov_param = model$cov_param, identification = model$identification,
               markers = if (identical(model$identification, "marker")) model$items[model$markers] else NULL,
               continuous = model$items[model$cont], resvar_mod = model$resvar_mod,
               experimental_cor = model$experimental_cor, sparse = 0)
  new <- suppressWarnings(do.call(mofia_model, utils::modifyList(args, list(...))))
  new$design <- model$design; new$moddata <- model$moddata
  new
}

#' @export
print.mofia_model <- function(x, ...) {
  nc <- sum(x$cont); cat("Moderated factor model\n")
  cat(sprintf("  N = %d; %d item(s) (%s); %d factor(s): %s\n", x$N, x$P,
              paste(c(if (x$P - nc) sprintf("%d categorical", x$P - nc), if (nc) sprintf("%d continuous", nc)), collapse = ", "),
              x$D, paste(x$factors, collapse = ", ")))
  cat(sprintf("  moderators: %s\n", if (x$Q) paste(x$moderators, collapse = ", ") else "none"))
  if (length(x$design$scaling)) cat(sprintf("  standardized: %s\n", paste(names(x$design$scaling), collapse = ", ")))
  cat(sprintf("  identification: %s; %d free parameters\n", x$identification, nrow(x$par)))
  tab <- table(factor(x$par$block, levels = c("tau", "nu0", "lambda0", "lambda1", "beta", "omega0", "omega",
                                              "gamma", "delta0", "delta", "r0", "r1")))
  lab <- c(tau = "thresholds", nu0 = "intercepts", lambda0 = "loadings", lambda1 = "loading moderation",
           beta = "threshold/intercept moderation", omega0 = "residual variances", omega = "residual-variance moderation",
           gamma = "factor-mean moderation", delta0 = "baseline factor variances", delta = "factor-variance moderation",
           r0 = "factor correlations", r1 = "correlation moderation")
  tab <- tab[tab > 0]
  for (b in names(tab)) cat(sprintf("    %-32s %d\n", lab[[b]], tab[[b]]))
  invisible(x)
}

# ---- internal: unpack the parameter vector -------------------------------------
.unpack <- function(theta, spec) {
  pt <- spec$par; P <- spec$P; Q <- spec$Q; D <- spec$D
  tau <- lapply(seq_len(P), function(j) theta[pt$block == "tau" & pt$item %in% j])
  L0 <- matrix(0, P, D); L1 <- array(0, c(P, D, Q)); B <- matrix(0, P, Q)
  np <- spec$npair
  G <- matrix(0, D, Q); Dl <- matrix(0, D, Q); r1 <- matrix(0, np, Q); r0 <- numeric(np)
  s <- pt$block == "lambda0"; L0[cbind(pt$item[s], pt$factor[s])] <- theta[s]
  s <- pt$block == "lambda1"; L1[cbind(pt$item[s], pt$factor[s], pt$mod[s])] <- theta[s]
  s <- pt$block == "beta";    B[cbind(pt$item[s], pt$mod[s])] <- theta[s]
  s <- pt$block == "gamma";   G[cbind(pt$factor[s], pt$mod[s])] <- theta[s]
  s <- pt$block == "delta";   Dl[cbind(pt$factor[s], pt$mod[s])] <- theta[s]
  s <- pt$block == "r0";      r0[pt$pair[s]] <- theta[s]
  s <- pt$block == "r1";      r1[cbind(pt$pair[s], pt$mod[s])] <- theta[s]
  Dl0 <- numeric(D); s <- pt$block == "delta0"; Dl0[pt$factor[s]] <- theta[s]
  nu0 <- numeric(P); s <- pt$block == "nu0"; nu0[pt$item[s]] <- theta[s]
  Om0 <- numeric(P); s <- pt$block == "omega0"; Om0[pt$item[s]] <- theta[s]
  Om <- matrix(0, P, Q); s <- pt$block == "omega"; Om[cbind(pt$item[s], pt$mod[s])] <- theta[s]
  if (identical(spec$identification, "marker")) L0[cbind(spec$markers, seq_len(D))] <- 1
  list(tau = tau, L0 = L0, L1 = L1, B = B, G = G, Dl = Dl, Dl0 = Dl0, r0 = r0, r1 = r1,
       nu0 = nu0, Om0 = Om0, Om = Om)
}

# person-level structural quantities: mu (N x D), s (N x D) factor means / SDs,
# rho (N x npair) correlations (cv: covariances for cov_param = "linear"), Lr:
# lower-triangular Cholesky factor of the correlation matrix ([[i]][[j]] N-vectors)
.lt_index <- function(i, j) i * (i - 1) / 2 + j                     # 1-based, i >= j
.chol_batch <- function(Rf, D, N) {                                   # Rf(i, j) -> N-vector, i > j
  L <- lapply(seq_len(D), function(i) vector("list", i)); ok <- rep(TRUE, N)
  for (j in seq_len(D)) {
    d <- 1 - (if (j > 1) Reduce(`+`, lapply(seq_len(j - 1), function(k) L[[j]][[k]]^2)) else 0)
    ok <- ok & (d > 1e-12); L[[j]][[j]] <- sqrt(pmax(d, 1e-12))
    if (j < D) for (i in (j + 1):D) {
      v <- Rf(i, j) - (if (j > 1) Reduce(`+`, lapply(seq_len(j - 1), function(k) L[[i]][[k]] * L[[j]][[k]])) else 0)
      L[[i]][[j]] <- v / L[[j]][[j]]
    }
  }
  list(L = L, ok = ok)
}
.structural <- function(pr, spec) {
  X <- spec$X; D <- spec$D; N <- nrow(X)
  mu <- X %*% t(pr$G)
  s <- exp((X %*% t(pr$Dl) + matrix(pr$Dl0, N, D, byrow = TRUE)) / 2)
  out <- list(mu = mu, s = s, rho = NULL, ok = TRUE)
  if (D >= 2) {
    pairs <- spec$pairs
    U <- sapply(seq_len(nrow(pairs)), function(p) pr$r0[p] + drop(X %*% pr$r1[p, ])); U <- matrix(U, N)
    if (spec$cov_param == "fisherz") out$rho <- tanh(U)
    else { out$cv <- U; out$rho <- U / (s[, pairs[, 1], drop = FALSE] * s[, pairs[, 2], drop = FALSE]) }
    if (any(!is.finite(out$rho)) || any(abs(out$rho) >= 1 - 1e-10)) { out$ok <- FALSE; return(out) }
    pidx <- matrix(0L, D, D); pidx[pairs] <- seq_len(nrow(pairs)); pidx <- pidx + t(pidx)
    ch <- .chol_batch(function(i, j) out$rho[, pidx[i, j]], D, N)
    if (!all(ch$ok)) { out$ok <- FALSE; return(out) }
    out$Lr <- ch$L; out$pidx <- pidx; out$pairs <- pairs
  }
  out
}
