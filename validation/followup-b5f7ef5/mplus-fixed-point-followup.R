# Targeted regression on retained B08 data/estimates. No optimization or simulation.
options(digits = 17)
# Follow-up copy of the review's mplus-fixed-point.R: base = the review's retained 08b9acb R sources,
# subject = the MoFIA working tree (run on a clean checkout of the commit under test).
# Outputs go to this folder; the review's historical outputs are not touched.
root <- "/Users/henrikgodmann/Desktop/GitHub/jaspSem/tempAI/ReworkMNLFA-engine-fix-review-20261111"
subject <- "/Users/henrikgodmann/Desktop/GitHub/MoFIA"
out_dir <- "/Users/henrikgodmann/Desktop/GitHub/MoFIA/validation/followup-b5f7ef5"
cat("Subject commit:", system2("git", c("-C", subject, "rev-parse", "HEAD"), stdout = TRUE),
    "| R/ clean:", !length(system2("git", c("-C", subject, "status", "--porcelain", "--", "R"), stdout = TRUE)), "\n")
stopifnot(!length(system2("git", c("-C", subject, "diff", "--name-only", "08b9acb", "HEAD", "--", "src"), stdout = TRUE)))
old <- "/Users/henrikgodmann/Desktop/GitHub/jaspSem/tempAI/ReworkMNLFA-B08-mplus-check-20261110"
lib <- "/Users/henrikgodmann/Library/Caches/org.R-project.R/R/renv/library/jaspSem-a86347ba/macos/R-4.5/aarch64-apple-darwin20"
.libPaths(c(lib, .libPaths()))
ns <- asNamespace("MoFIA")
cat("R:", R.version.string, "\nCompiled backend:", find.package("MoFIA"), "\n")
print(packageDescription("MoFIA")[c("Version", "RemoteSha", "Built")])
cat("C++ source unchanged between review base and subject; R code sourced by revision.\n")
load_source <- function(which) {
  en <- new.env(parent = ns)
  for (f in list.files(if (which == "base") file.path(root, "source", "base", "R") else file.path(subject, "R"), pattern = "\\.R$", full.names = TRUE)) sys.source(f, en)
  en$mofia <- function(...) stop("No refits permitted in fixed-point regression")
  en
}
envs <- lapply(c("base", "fixed"), load_source); names(envs) <- c("base", "fixed")
ref <- readRDS(file.path(old, "reference.rds"))
d <- read.table(file.path(old, "b08_rep40001.dat")); names(d) <- ref$columns
mods <- c("x1", "groupg1", "groupg2", "x1:groupg1", "x1:groupg2")
d[mods] <- d[c("x1", "g1", "g2", "xg1", "xg2")]
its <- paste0("y", 1:6); fref <- ref$fits$full
fits <- lapply(envs, function(en) {
  m <- en$mofia_model(d, loadings = list(f1 = its), continuous = its,
                     moderators = mods, standardize = FALSE,
                     loading_mod = list(y3 = "x1"), threshold_mod = list(y4 = mods),
                     resvar_mod = list(y5 = "x1"))
  stopifnot(identical(m$par$name, names(fref$est)), m$N == 500L)
  nodes <- en$.mofia_nodes(fref$nq, m$D)
  adapt <- en$.mofia_adapt(fref$est, m, "logit", nodes)
  ll <- as.numeric(en$.mofia_ll(fref$est, m, "logit", nodes, adapt = adapt))
  stopifnot(abs(ll - fref$loglik) < 1e-9)
  fit <- list(model = m, coefficients = fref$est, link = "logit", nq = fref$nq, adapt = adapt,
              loglik = ll, npar = length(fref$est))
  en$.mofia_se(fit)
})
summary <- list()
for (which in names(fits)) {
  f <- fits[[which]]
  cat(which, "status", f$status, "converged", f$converged, "loglik", f$loglik,
      "scaling", f$scaling, "\n")
  stopifnot(identical(f$status, "ok"), isTRUE(f$converged))
}
for (kind in c("exact", "reduced")) for (se in c("ML", "MLR")) {
  A <- fits$base$vcov[[kind]][[se]]; B <- fits$fixed$vcov[[kind]][[se]]
  stopifnot(identical(is.null(A), is.null(B)))
  if (is.null(A)) next
  label <- paste(kind, se, sep = "_")
  err <- max(abs(B - A)); serr <- max(abs(sqrt(diag(B)) / sqrt(diag(A)) - 1))
  cat("New vs base", label, "covariance max absolute", err, "SE max relative", serr, "\n")
  stopifnot(serr < 1e-7)
  if (!is.null(fref$vcov[[label]])) {
    referr <- max(abs(A - fref$vcov[[label]]))
    cat("Reconstructed base vs retained", label, "max covariance difference", referr, "\n")
    stopifnot(referr < 1e-7)
  }
  summary[[label]] <- data.frame(method=label, max_abs_covariance_change=err, max_relative_SE_change=serr)
}
tab <- read.csv(file.path(old, "results", "mplus_covariances.csv"))
labels <- unname(ref$mapping[names(fref$est)])
T <- rbind(diag(length(labels)), matrix(0, 3, length(labels)))
dimnames(T) <- list(c(labels, "sg3", "sg1", "sg2"), labels)
T[c("sg3", "sg1", "sg2"), "dx"] <- 1
T["sg1", "dxg1"] <- 1; T["sg2", "dxg2"] <- 1
out <- list()
for (se in c("ML", "MLR")) {
  V <- T %*% fits$fixed$vcov$exact[[se]] %*% t(T)
  z <- tab[tab$model == "full" & tab$estimator == se, ]
  z$fixed_covariance <- V[cbind(match(z$label1, rownames(V)), match(z$label2, colnames(V)))]
  z$fixed_difference <- z$mplus_covariance - z$fixed_covariance
  diagz <- z[z$label1 == z$label2, ]
  seerr <- max(abs(sqrt(diagz$mplus_covariance) - sqrt(diagz$fixed_covariance)))
  cat("Saved Mplus vs fixed", se, "unique covariance entries", nrow(z),
      "max covariance difference", max(abs(z$fixed_difference)), "max SE difference", seerr, "\n")
  stopifnot(!anyNA(z), nrow(z)==741L, seerr < 0.001)
  out[[se]] <- z
}
cat("MLR scaling difference from Mplus printed 0.9752:", fits$fixed$scaling[["exact"]] - 0.9752, "\n")
stopifnot(abs(fits$fixed$scaling[["exact"]] - 0.9752) < 0.001)
write.csv(do.call(rbind, summary), file.path(out_dir, "fixed-point-changes.csv"), row.names=FALSE)
write.csv(do.call(rbind, out), file.path(out_dir, "mplus-fixed-point-covariances.csv"), row.names=FALSE)
saveRDS(fits, file.path(out_dir, "fixed-point-results.rds"))
cat("PASS: fixed-point B08 exact/reduced regression and saved-Mplus exact ML/MLR comparison.\n")
cat("Not a refit, new Mplus execution, or JASP whole-chain Check 3.\n")
