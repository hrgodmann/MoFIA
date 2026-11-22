################################################################################
# Compare the MoFIA fits with the Mplus fits
#
# What this script does
#   For each of the 12 models (4 model specifications x configural / metric /
#   scalar) it puts the MoFIA result next to the Mplus result and reports:
#     - the number of cases (N) and of free parameters,
#     - the maximized log-likelihood,
#     - the largest difference between the parameter estimates,
#     - the range of the ratio of the standard errors (Mplus / MoFIA),
#   and, per model specification, the scaled likelihood-ratio tests between
#   the nested invariance models.
#
# What you need
#   1. MoFIA/MoFIA_Comp_results.RData  : saved at the end of MoFIA/MoFIA_Comp.R
#   2. MPLUS/Model*_*.out              : the Mplus outputs of MPLUS/Model*_*.inp
#
# Note on precision
#   Mplus prints results with 3 decimals. Differences of about 0.0005 in the
#   estimates and standard-error ratios of about 0.99 to 1.01 are therefore
#   rounding, not real differences.
################################################################################

library(MoFIA)   # needed for vcov(), mofia_lrt() on the saved fits

############### Paths (adjust to your copy of the repository)
comparisons_dir <- "/Users/henrikgodmann/Desktop/GitHub/MoFIA/inst/examples/comparisons"
mofia_results   <- file.path(comparisons_dir, "MoFIA", "MoFIA_Comp_results.RData")
mplus_dir       <- file.path(comparisons_dir, "MPLUS")


############### Load the MoFIA fits
# The .RData file contains the whole workspace of MoFIA_Comp.R. We load it into
# its own environment, so that it does not overwrite anything here.
mofia <- new.env()
load(mofia_results, envir = mofia)


############### The four model specifications, as written in the Mplus inputs
# Mplus allows names of at most 8 characters, so some items and one factor
# have a shorter name in Mplus.
mplus_name <- c(HealingHands = "HealHand", LifeCoursePossessions = "LifePoss",
                Graphology = "Graphol", Dowsing = "Dowsing", MindReading = "MindRead",
                Willpower = "Willpow", SeeingFuture = "SeeFut", SeeingEvents = "SeeEvent",
                SpiritsGhosts = "Spirits", GnomesElves = "Gnomes", Aliens = "Aliens",
                DevilPossession = "DevilPos", Sorcery = "Sorcery", HauntedPlaces = "Haunted",
                Depression = "Depress")
to_mplus <- function(x) ifelse(x %in% names(mplus_name), mplus_name[x], x)

# In the Mplus inputs, the parameters defined in MODEL CONSTRAINT are named by
# the position of the item ("j") and of the factor ("f") in the lists below,
# and by one letter per moderator, e.g. L2_1_A = effect of Age on the loading
# of item 2 on factor 1. The item order below is the order of USEVARIABLES.
models <- list(
  Model1 = list(items     = c("HealHand", "LifePoss", "Graphol", "Dowsing", "MindRead", "Willpow",
                              "SeeFut", "SeeEvent", "Spirits", "Gnomes", "Aliens", "DevilPos",
                              "Sorcery", "Haunted"),
                factors   = c("EHB", "UB"),
                anchors   = c("HealHand", "Spirits"),
                ordinal   = TRUE),
  Model2 = list(items     = c("x1", "x2", "x3", "x9", "x7", "x8"),
                factors   = c("Visual", "Speed"),
                anchors   = c("x1", "x7"),
                ordinal   = FALSE),
  Model3 = list(items     = c("x1", "x2", "x3", "x9", "x4", "x5", "x6", "x7", "x8"),
                factors   = c("Visual", "Textual", "Speed"),
                anchors   = c("x1", "x4", "x7"),
                ordinal   = FALSE),
  Model4 = list(items     = c("i4", "i19", "i34", "i49", "i9", "i24", "i39", "i54"),
                factors   = c("Anxiety", "Depress"),
                anchors   = c("i4", "i9"),
                ordinal   = TRUE)
)
moderator_letter <- c(Age = "A", Male = "M", school = "S", Gender = "G")
types <- c("conf", "metric", "scalar")


############### Read an Mplus output file

# The last number on the first line that contains `label`
mplus_value <- function(lines, label) {
  line <- grep(label, lines, fixed = TRUE, value = TRUE)[1]
  as.numeric(tail(strsplit(trimws(line), "\\s+")[[1]], 1))
}

# The MODEL RESULTS table as a data frame with one row per parameter:
# "where" is the section and row name as Mplus prints them, e.g. "VISUAL BY X1"
read_mplus_results <- function(lines) {
  first <- grep("^MODEL RESULTS", lines) + 1
  # The table ends at the next heading in capital letters, e.g. "QUALITY OF NUMERICAL RESULTS"
  headings <- grep("^[A-Z][A-Z ]+$", lines)
  last  <- min(headings[headings > first]) - 1
  section <- ""
  rows <- data.frame()
  for (line in lines[first:last]) {
    if (trimws(line) == "") next
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    is_parameter <- startsWith(line, "    ") && length(parts) == 5
    if (is_parameter) {
      rows <- rbind(rows, data.frame(where    = toupper(paste(section, parts[1])),
                                     estimate = as.numeric(parts[2]),
                                     se       = as.numeric(parts[3])))
    } else {
      section <- trimws(line)   # e.g. "VISUAL   BY", "Thresholds", "New/Additional Parameters"
      section <- gsub("\\s+", " ", section)
    }
  }
  rows
}

read_mplus <- function(file) {
  lines <- readLines(file)
  list(N       = mplus_value(lines, "Number of observations"),
       npar    = mplus_value(lines, "Number of Free Parameters"),
       logLik  = mplus_value(lines, "H0 Value"),
       scaling = mplus_value(lines, "H0 Scaling Correction Factor"),
       results = read_mplus_results(lines))
}


############### Where Mplus prints a MoFIA parameter
# MoFIA names its parameters like "lambda1[x2,Visual,school]". This function
# returns the matching Mplus section and row, e.g. "NEW/ADDITIONAL PARAMETERS L2_1_S".
mplus_row <- function(mofia_name, model, type) {
  block <- sub("\\[.*", "", mofia_name)                               # e.g. "lambda1"
  args  <- strsplit(gsub(".*\\[|\\]", "", mofia_name), ",")[[1]]      # e.g. c("x2", "Visual", "school")
  args  <- to_mplus(args)
  j <- function(item)   match(item, model$items)                      # item position
  f <- function(factor) match(factor, model$factors)                  # factor position
  m <- function(mod)    moderator_letter[[mod]]                       # moderator letter
  new <- "New/Additional Parameters"
  moderated_loading <- type == "conf" && !(args[1] %in% model$anchors)

  where <- switch(block,
    lambda0 = if (moderated_loading) sprintf("%s L%d_%d_0", new, j(args[1]), f(args[2]))
              else                   sprintf("%s BY %s", args[2], args[1]),
    lambda1 = sprintf("%s L%d_%d_%s", new, j(args[1]), f(args[2]), m(args[3])),
    tau     = sprintf("Thresholds %s$%s", args[1], args[2]),
    nu0     = sprintf("Intercepts %s", args[1]),
    beta    = sprintf("%s ON %s", args[1], args[2]),
    omega0  = sprintf("%s W%d_0", new, j(args[1])),
    omega   = sprintf("%s W%d_%s", new, j(args[1]), m(args[2])),
    gamma   = sprintf("%s ON %s", args[1], args[2]),
    delta   = sprintf("%s D%d_%s", new, f(args[1]), m(args[2])),
    zcor0   = sprintf("%s Z%d%d_0", new, f(args[1]), f(args[2])),
    zcor1   = sprintf("%s Z%d%d_%s", new, f(args[1]), f(args[2]), m(args[3])))
  toupper(where)
}


############### Compare the 12 models
# Standard errors: MLR in both programs. For MoFIA, the exact information matrix
# for continuous items and the reduced information matrix for ordinal items,
# which reproduces the standard errors Mplus reports for the ordinal models.
comparison <- data.frame()
for (model_name in names(models)) {
  model <- models[[model_name]]
  information <- if (model$ordinal) "reduced" else "exact"
  for (type in types) {
    fit   <- get(paste0("fit_", model_name, "_", type), envir = mofia)
    mplus <- read_mplus(file.path(mplus_dir, paste0(model_name, "_", type, ".out")))

    mofia_estimate <- coef(fit)
    mofia_se       <- sqrt(diag(vcov(fit, se = "MLR", information = information)))
    rows           <- match(sapply(names(mofia_estimate), mplus_row, model = model, type = type),
                            mplus$results$where)
    if (anyNA(rows)) stop("Not found in Mplus output: ", paste(names(mofia_estimate)[is.na(rows)], collapse = ", "))

    estimate_difference <- abs(mofia_estimate - mplus$results$estimate[rows])
    se_ratio            <- mplus$results$se[rows] / mofia_se

    comparison <- rbind(comparison, data.frame(
      model = model_name, type = type,
      N_MoFIA = fit$nobs, N_Mplus = mplus$N,
      npar_MoFIA = fit$npar, npar_Mplus = mplus$npar,
      logLik_MoFIA = round(fit$loglik, 3), logLik_Mplus = mplus$logLik,
      max_abs_estimate_difference = round(max(estimate_difference), 4),
      se_ratio_min = round(min(se_ratio), 3), se_ratio_max = round(max(se_ratio), 3),
      information_MoFIA = information))
  }
}
print(comparison)


############### Scaled likelihood-ratio tests between the nested models
# MoFIA: mofia_lrt() (Satorra-Bentler 2001 scaled difference test).
# Mplus does not print this test; it is computed from the log-likelihoods,
# numbers of parameters and MLR scaling correction factors of the two models:
#   cd = (p0 * c0 - p1 * c1) / (p0 - p1),   TRd = -2 * (L0 - L1) / cd
# with model 0 the smaller (nested) model.
satorra_bentler <- function(small, large) {
  cd <- (small$npar * small$scaling - large$npar * large$scaling) / (small$npar - large$npar)
  -2 * (small$logLik - large$logLik) / cd
}

lr_tests <- data.frame()
for (model_name in names(models)) {
  information <- if (models[[model_name]]$ordinal) "reduced" else "exact"
  for (pair in list(c("metric", "conf"), c("scalar", "metric"))) {          # c(smaller, larger)
    fit_small <- get(paste0("fit_", model_name, "_", pair[1]), envir = mofia)
    fit_large <- get(paste0("fit_", model_name, "_", pair[2]), envir = mofia)
    mofia_test <- mofia_lrt(fit_small, fit_large, scaled = TRUE, information = information)
    mplus_small <- read_mplus(file.path(mplus_dir, paste0(model_name, "_", pair[1], ".out")))
    mplus_large <- read_mplus(file.path(mplus_dir, paste0(model_name, "_", pair[2], ".out")))

    lr_tests <- rbind(lr_tests, data.frame(
      model = model_name, test = paste(pair[1], "vs", pair[2]),
      df = mofia_test$df,
      chisq_MoFIA = round(mofia_test$statistic, 3),
      chisq_Mplus = round(satorra_bentler(mplus_small, mplus_large), 3)))
  }
}
print(lr_tests)

# To save both tables, uncomment:
# write.csv(comparison, file.path(comparisons_dir, "comparison_models.csv"), row.names = FALSE)
# write.csv(lr_tests,   file.path(comparisons_dir, "comparison_lr_tests.csv"), row.names = FALSE)
