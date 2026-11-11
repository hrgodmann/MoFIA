# Follow-up to b5f7ef5 (review findings F1–F4): evidence

- **Subject:** commit `2ce8b51a5739640f1cd6e828246361648dd007fa` (tree
  `27e0fdb7330b563f75f10d5ce3bcdb47fad63542`), branch `engine-review-fixes`, on top of
  `b5f7ef5` and `08b9acb`. All scripts below were run on a clean `R/` at that commit.
  `src/` is unchanged since `08b9acb` (asserted in `mplus-fixed-point-followup.R`).
- **Review being answered:** `jaspSem/tempAI/ReworkMNLFA-engine-fix-review-20261111/REVIEW.md`
  and `FOLLOWUP-FOR-ENGINE-AGENT.md`. That folder was only read; nothing in it was changed.
- **Author:** Claude (Anthropic), the same session that wrote `b5f7ef5`; not an
  independent check.
- **Environment:** R 4.5.3 (aarch64-apple-darwin20), Apple clang without OpenMP (serial
  build), Rcpp 1.1.1.1.1, statmod 1.5.2, testthat 3.3.2, numDeriv 2016.8.1.1, lavaan
  0.6.21, mirt 1.48. Package and mirt installed in a session scratch library; lavaan
  and the other packages from the user library.

## Build and check

| What | Result |
|---|---|
| `R CMD build` of the commit | `MoFIA_0.1.0.9000.tar.gz`, sha256 `bf11cb935b5772467123f8ad46e8b3bafdc2e86a6e582b9bdc99d37b362eee26` (`*.sha256`; the tarball is kept in this folder but ignored by git) |
| `NOT_CRAN=true _R_CHECK_FORCE_SUGGESTS_=true R CMD check --as-cran --run-donttest` | 0 errors, 0 warnings, 2 NOTEs (repository URL 404 because private / development version; HTML Tidy and V8 unavailable) — `R-CMD-check.log` |
| testthat in the check | 368 passed, 0 failed, 1 skipped: the multithreaded comparison (no OpenMP in this build; only the serial fallback is covered) — `testthat.Rout` |

## Findings

| Finding | Change | Evidence |
|---|---|---|
| F1 reduced-information cancellation | Reduced information is classified and inverted on the exact information's diagonal scale (`.info_inverse(I_red, scale = sqrt(abs(diag(I_exact))))`); exact information and MLF unchanged. No ridge, clipping, pseudoinverse or fallback. | `reduced-boundary-followup.R/.log`: the review's boundary probe now gives reduced covariance `NULL` and scaling `NA`, as in 08b9acb. `info-probe.R/.log` (the review's `information/probe.R`, unchanged): 54 passed, 0 failed (52/2 on b5f7ef5). Test F1 covers singular, clearly positive and indefinite reduced matrices under three unit/sign transformations. |
| F2 origin-dependent false alias | `.mean_alias` uses centred/scaled moderators with the products expanded exactly, snaps residuals of moderators lying in the span to zero, and takes the rank from singular values (no Gram matrix). `.check_moderators` ranks centred/scaled columns. | `identification-final-checks-followup.R/.log`: the review's seven contract checks pass and raw `x = 10000 + w` now gets the same warning status as centred `w`. `f2_origin_probe.R/.log`: origins 0, 1e4, 1e9 agree; binary (codings 0/1, 1/2, 0/1000, −5e4, 1e9), offset interactions and `x + I(x^2)` remain errors; anchor/no-mean controls pass. |
| F3 legacy row metadata | Models without `$cases` (0.1.0) keep `rownames(model$Y)`; `"row_index"` is `NULL`. | Test F3 with `tests/testthat/fixtures/legacy-0.1.0-fit.rds`, a fit made by 08b9acb itself (`make-legacy-0.1.0.R`) with interleaved and trailing exclusions and custom row names; scores/predictions also equal those 0.1.0 computed. |
| F4 documented assignment | `mofia_scores()` docs initialize a full-length `NA` column first. | Test F4 (trailing exclusion, reused column with stale values). It passes on b5f7ef5 too: the defect was the documentation, not the code. |

`b5f7ef5_vs_followup_tests.log` / `2ce8b51_vs_followup_tests.log`: F1 and F2 error and F3
fails on b5f7ef5; all pass at 2ce8b51. Earlier R1–R5 regressions pass on both.

## Mplus reference (B08, exact information)

`mplus-fixed-point-followup.R/.log` is the review's `mplus-fixed-point.R` with the
subject R sources taken from the MoFIA working tree at 2ce8b51 (base: the review's
retained 08b9acb sources) and outputs written here. At the saved B08 estimates, no refit:

- log-likelihood base and subject −3691.0095275777958;
- subject vs base, exact and reduced ML/MLR SEs: max relative change 3.6e−14;
- saved Mplus vs subject, exact ML SE max difference 4.39e−7, MLR 1.052e−6 (741
  covariance entries each, max covariance difference 3.85e−7); MLR scaling
  0.975184 vs printed 0.9752.

Same values as the review's run on b5f7ef5. This is the exact-information B08 reference,
not the reduced-information categorical harness of the planned JASP Check 3, which was
not run. No Mplus was executed and no model was refitted for it.

## Not done

No calibration, enabled-OpenMP run, Desktop or JASP whole-chain check. The jaspSem pin,
frozen pack, S2a records and earlier Mplus outputs were not touched.
