# Engine audit of 11 November 2026: scripts, logs and results

This folder keeps the evidence for the audit of MoFIA 0.1.0 and for the corrections that
followed it (NEWS, 0.1.0.9000). It is not part of the built package (`.Rbuildignore`).

- **Audited source:** `main` at commit `08b9acb364570683baf63a98e761f061bd3aabac` (version
  0.1.0), clean tree.
- **Auditor:** Claude (Anthropic) as a task-agent for Henrik Godmann, in the same session
  that later implemented the corrections. This is therefore not an independent review of
  the corrections.
- **Environment:** R 4.5.3 (aarch64-apple-darwin20), Apple clang without OpenMP (MoFIA
  ran serially), lavaan 0.6.21, mirt 1.48, numDeriv 2016.8.1.1, testthat 3.3.2.

## Scripts (`scripts/`)

Each script simulates its own data with a fixed seed or uses `parabelief`. Run with
`Rscript scripts/<name>.R` and the package to test installed. `f1_age2.R` and
`run_tests_table.R` take the library paths as arguments (see their first lines).

| Script | What it checks |
|---|---|
| `e1_ident_scale.R` | E1: binary moderator with mean, all thresholds and all loadings moderated (exact alias); E2: age in years plus its square against the standardized model |
| `e3_external.R` | E3a: continuous items vs lavaan multigroup (logLik, estimates, ML/MLR SEs); E3b: graded items vs mirt `multipleGroup` |
| `e4_closedform.R` | Two factors, continuous items, continuous moderators: closed-form likelihood, an independent `optim` fit and numerical-Hessian SEs |
| `e5_integrate.R` | Ordinal items, moderated loadings/thresholds/mean/variance: per-person likelihood by `integrate()`, logit and probit |
| `e6_calibration.R` | 1000 replications at N = 250 and N = 1000: bias, CI coverage (ML, MLR), type-I error of Wald, scaled/unscaled LRT, robust/ML score tests |
| `f1_age2.R` | After the fix: age in years plus its square against the standardized model, coefficients and covariances compared after the exact reparameterization |
| `run_tests_table.R` | Runs `tests/testthat/test-review3.R` and `test-comparisons.R` against a given installation and tabulates the results |

## Logs (`logs/`)

- `original-08b9acb/`: output at the audited commit. `e1`–`e5` and `R-CMD-check.log` were
  recorded during the audit; `e6_calibration.log` is the calibration run whose replicate
  results are in `results/calib_N*_R1000.rds`. `new-regression-tests.log` runs the new
  regression tests against the original commit: 13 of 15 fail as expected, and the two
  controls pass (a genuinely unidentified model is still reported as singular; the
  refused alias is an exact likelihood invariance).
- `fixed/`: the same scripts and the new tests run against the corrected package. The
  calibration (`e6`) was not rerun, and no new calibration campaign was requested: for
  well-scaled models the corrections leave the likelihood and estimates unchanged and the
  standard errors equal up to numerical-differentiation error (`cmp_old_new.log`: at most
  2.5e-10 relative, from the moderator-scaled difference step). Data with cases that have
  no item response can change, because those cases no longer enter the standardization
  and classification of the moderators.

## Limits

- The comparisons cover the configurations listed above, not every model MoFIA fits
  (e.g. not three or more factors, sparse categories or moderated correlations beyond two
  factors).
- The reference log-likelihood `-8394.7028` in `tests/testthat/test-fit.R` was recorded by
  MoFIA itself when that test was added; its original software and script are not known.
  It is a regression baseline, not an independent comparison.
- OpenMP: this machine's build is serial, so only the serial fallback was tested; the
  multithreaded comparison in `test-review3.R` (R5) runs only on an OpenMP build.
