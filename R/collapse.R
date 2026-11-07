#' Collapse sparse categories of ordinal items
#'
#' Merges categories with fewer than `min_count` responses into an adjacent
#' category, one at a time, until all categories are large enough (or two
#' remain). A sparse category is merged with whichever neighbour has fewer
#' responses (with ties, the one towards the middle of the scale). The
#' remaining categories are recoded as consecutive integers starting at the
#' item's lowest observed value.
#'
#' Sparse categories make threshold estimates unstable and can slow
#' convergence; collapsing them changes the measurement model, so report it.
#'
#' @param data A data frame.
#' @param items Names of the ordinal items to check (integer-coded).
#' @param min_count Minimum number of responses per category.
#' @return The data frame with the collapsed items. The recoding is attached as
#'   attribute `"collapsed"`: a data frame with columns `item`, `from` and `to`
#'   for the items that changed.
#' @examples
#' data(parabelief)
#' d <- mofia_collapse(parabelief, c("GnomesElves", "Aliens"), min_count = 10)
#' attr(d, "collapsed")
#' table(d$GnomesElves)
#' @export
mofia_collapse <- function(data, items, min_count = 10) {
  miss <- setdiff(items, names(data)); if (length(miss)) stop("Not in the data: ", paste(miss, collapse = ", "))
  log <- list()
  for (it in items) {
    x <- data[[it]]
    if (!is.numeric(x) || any(x[!is.na(x)] != round(x[!is.na(x)]))) stop("Item ", it, " must be integer-coded.")
    lv <- sort(unique(x[!is.na(x)]))
    groups <- as.list(lv); counts <- as.numeric(table(factor(x, levels = lv)))
    while (length(counts) > 2 && min(counts) < min_count) {
      i <- which.min(counts); n <- length(counts)
      j <- if (i == 1) 2 else if (i == n) n - 1 else {
        l <- counts[i - 1]; r <- counts[i + 1]
        if (l < r) i - 1 else if (r < l) i + 1 else if (i - 1 >= n - i) i - 1 else i + 1
      }
      keep <- min(i, j); drop <- max(i, j)
      groups[[keep]] <- sort(c(groups[[keep]], groups[[drop]])); counts[keep] <- counts[i] + counts[j]
      groups <- groups[-drop]; counts <- counts[-drop]
    }
    if (length(groups) == length(lv)) next
    new <- min(lv) + seq_along(groups) - 1
    map <- unlist(lapply(seq_along(groups), function(g) stats::setNames(rep(new[g], length(groups[[g]])), groups[[g]])))
    data[[it]] <- unname(map[as.character(x)])
    log[[it]] <- data.frame(item = it, from = as.numeric(names(map)), to = unname(map))
  }
  attr(data, "collapsed") <- if (length(log)) do.call(rbind, log) else
    data.frame(item = character(0), from = numeric(0), to = numeric(0))
  rownames(attr(data, "collapsed")) <- NULL
  data
}
