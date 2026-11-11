args <- commandArgs(TRUE); .libPaths(c(args[1], args[2], .libPaths()))
suppressMessages(library(testthat))
cat("MoFIA", as.character(packageVersion("MoFIA")), "from", find.package("MoFIA"), "\n")
source("helper-sim.R")
for (f in c("test-review3.R", "test-comparisons.R")) {
  r <- as.data.frame(test_file(f, package = "MoFIA", load_package = "installed", reporter = "silent", stop_on_failure = FALSE))
  r$result <- ifelse(r$error, "ERROR", ifelse(r$failed > 0, "FAIL", ifelse(r$skipped, "SKIP", "pass")))
  cat("\n==", f, "==\n"); print(r[, c("test", "nb", "failed", "result")], right = FALSE, row.names = FALSE)
}
