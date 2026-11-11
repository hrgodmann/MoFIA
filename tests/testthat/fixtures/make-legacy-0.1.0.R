# Creates legacy-0.1.0-fit.rds: a model and fit made with MoFIA 0.1.0 (commit 08b9acb),
# which did not record $cases, from data with interleaved exclusions and custom row
# names; plus the scores and predictions 0.1.0 computed for it. Run with 0.1.0 installed
# in the library given as first argument:
#   Rscript make-legacy-0.1.0.R <library with MoFIA 0.1.0>
.libPaths(c(commandArgs(TRUE)[1], .libPaths()))
library(MoFIA); stopifnot(packageVersion("MoFIA") == "0.1.0")
set.seed(11); N <- 60; eta <- rnorm(N)
d <- data.frame(y1 = as.integer(eta + rlogis(N) > 0), y2 = as.integer(1.2 * eta + rlogis(N) > 0),
                y3 = as.integer(0.8 * eta + rlogis(N) > 0), y4 = as.integer(eta + rlogis(N) > 0.3),
                age = round(runif(N, 20, 70)), g = rep(0:1, N / 2))
rownames(d) <- sprintf("id%03d", seq_len(N) * 7)
d$age[c(2, 7)] <- NA
d[c(4, 9, 60), c("y1", "y2", "y3", "y4")] <- NA          # includes a trailing exclusion
d$y2[5] <- NA
m <- suppressWarnings(mofia_model(d, items = c("y1", "y2", "y3", "y4"), moderators = c("age", "g"),
                                  standardize = TRUE, var_mod = FALSE))
f <- mofia(m)
saveRDS(list(fit = f, data = d, kept = setdiff(seq_len(N), c(2, 7, 4, 9, 60)),
             scores = mofia_scores(f), expected = predict(f), version = "0.1.0", commit = "08b9acb"),
        "legacy-0.1.0-fit.rds", version = 2)
