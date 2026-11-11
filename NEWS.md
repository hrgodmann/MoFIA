# MoFIA 0.1.0.9000

## Numerical stability and compatibility

* The reduced information matrix is judged and inverted on the scale of the
  exact information's diagonal, not its own: where the constraint curvature
  cancels the exact information, the remaining finite-difference noise is no
  longer presented as information (reduced standard errors and scaling factors
  are then unavailable, as in 0.1.0). Exact information, MLF and well-scaled
  reduced information are unchanged.
* The mean/DIF alias check works on centred and scaled moderators with the
  products expanded exactly, so a continuous moderator far from its origin
  (e.g. 10000 + small variation) is classified like the centred one (identified
  only through the loading moderation) instead of being refused as exactly
  unidentified. Binary, interaction and quadratic aliases remain errors under
  any coding. The linear-dependence check of the moderators is conditioned the
  same way.
* For models created with 0.1.0 (no `$cases`), scores and predictions keep the
  row names stored in the model and no longer claim original positions
  (`"row_index"` is `NULL`).
* The documented way to attach scores to the data starts from a full-length
  `NA` column, so excluded cases (also trailing ones) stay missing.

## Inference and data handling

* Positive definiteness of the information matrix (status `"ok"`, `"saddle"`
  or `"singular"`) is judged on the matrix scaled to unit diagonal, and all
  inverses (exact and reduced ML/MLR, MLF), sandwich covariances, scaling
  factors and the Newton decrement are computed in these scaled coordinates
  and transformed back. A moderator in large units (e.g. age in years and its
  square) was wrongly reported as `"singular"`, with no standard errors or
  tests. Genuine singular and saddle points are still reported. The raw
  spectrum stays in `$information_eigen`; the scaled one is in
  `$information_eigen_scaled`. Wald tests check the covariance of the tested
  parameters on the correlation scale; score tests compare the information
  left after conditioning with the candidates' own information, so that
  neither depends on parameter units. Numerical second derivatives use steps
  relative to the scale of the moderator for moderation effects.
* Factor mean moderation together with threshold moderation of all items of
  the factor is refused also when loadings are moderated, if the products of
  the loading and mean moderators are absorbed by the thresholds (always for a
  binary moderator, z * z = z; for binary moderators and all their
  interactions; for `x` and `I(x^2)` both in the threshold moderation). Such
  models were only warned about, with a message suggesting identification
  through the functional form; they are exactly unidentified. Patterns
  identified only through moderated loadings now say so. Identification
  errors and warnings have classes `"mofia_identification_error"` and
  `"mofia_identification_warning"` and carry a data frame `$issues`; the
  warnings are also stored in `model$identification_issues`.
* Cases with missing moderator values and cases without any observed item
  response are both removed before moderator types, standardization and
  data-dependent formula terms are determined, so that these describe exactly
  the fitted cases (a case without responses could change the standardization
  or the classification of a moderator).
* `mofia_model()` records which input rows were fitted (`$cases`: positions,
  row names, number of input rows and the excluded rows of both kinds).
  `mofia_scores()` and `predict()` keep the number and type of their rows, but
  now carry the row names of the data and attribute `"row_index"` (positions
  in the data or in `newdata`); posterior covariance arrays are named likewise.
* `mofia(threads = n)` with `n > 1` on a build without OpenMP support says that
  it computes serially; the fit records requested and used threads in
  `$threads`.
* The citation year comes from the release date in DESCRIPTION, no longer from
  the current date. The documentation states that threshold moderation is one
  shift per item and moderator (no category-specific threshold moderation),
  and how to install from the private repository.
* New tests: the corrections above, and independent comparisons with lavaan
  and mirt (when installed), a closed-form likelihood and `integrate()`. The
  reference log-likelihood of the parabelief test is documented as a
  regression baseline.

# MoFIA 0.1.0

## Numerical accuracy and prediction

* Quadrature accuracy is checked after every fit: the log-likelihood at the
  estimates is recomputed with a finer rule, and the model is refitted with the
  finer rule when it changes by more than `accuracy` (default 0.005). Very
  discriminating items can need many more nodes than the default; the check is
  stored in `$quadrature` and printed.
* With only continuous items the posterior of the factors is normal, and the
  adaptive rule now uses it in closed form: the integration is exact also for
  very precise items (fewer nodes by default, faster fits). In mixed models it
  is the starting point of the adaptation. Adaptation that does not settle now
  gives a warning instead of a silently wrong likelihood.
* `predict()` computes marginal probabilities exactly for probit items and
  with an accurate one-dimensional rule for logit items (it used a coarse
  fixed grid that could be off by several percentage points).
* `mofia_score_test()` works again when all candidate effects together would
  not be identified (e.g. one binary moderator); each test is still checked on
  its own.
* `mofia_wald()` and `mofia_score_test()` refuse fits that have not converged;
  summaries and `confint()` of such fits carry a warning.
* `mofia_lrt(sb2010 = TRUE)` refuses models that are not nested by parameter
  names, instead of failing.
* Item plots work for items on the second or later factor and for
  cross-loading items, and no longer print a note about standard errors.
* `mofia_scores()` works with the R reference engine.

## Inference and input validation

* Robust (MLR) standard errors, scaled likelihood-ratio tests and Wald tests now
  use the exact observed information by default. The "reduced" information
  (which drops the curvature of the nonlinear parameter transformations and was
  the default for categorical items) is consistent only under correct
  specification; it remains available as `information = "reduced"`.
* `mofia_model()` refuses exactly unidentified moderation patterns (factor mean
  and the thresholds/intercepts of all its items moderated by the same
  moderator; factor variance and all loadings moderated by the same binary
  moderator) and warns about patterns identified only through the functional
  form. Constant, non-finite and linearly dependent moderators are refused.
* `mofia()` declares convergence only if the information matrix is positive
  definite. Saddle points and singular solutions get `status = "saddle"` or
  `"singular"`, a warning, and no standard errors or tests.
* Non-positive variance estimates give `NA` standard errors, intervals and
  plot bands instead of zero.
* Category probabilities and their derivatives are computed without
  cancellation or overflow in the tails (both links; C++ engine, R reference
  engine and predictions). The probability floor is gone.
* `mofia_lrt()` and `anova()` check that the models were fitted to the same
  cases and responses, with the same link and compatible moderator values, and
  that both converged; a clearly negative likelihood-ratio statistic is an
  error (no test in `anova()`).
* Cases without any observed item response are removed with a warning, so they
  no longer count in N and BIC.
* Specification matrices and arrays must have the exact dimensions; named ones
  are aligned by name, unnamed lists are refused, nothing is recycled.
* Data-dependent formula transformations (`scale()`, `poly()`) keep their
  training values for new data; commas in moderator names become semicolons.
* `print(summary(fit, se = "MLF"))` works.
* `mofia_simulate()` returns the moderator variables on their original scale
  (new argument `newdata`; `X` for design-scale values).

## Features

* Model syntax for the loading pattern (`"f1 =~ y1 + y2"`) and moderator
  formulas (`~ age * gender + group`) with automatic dummy coding,
  interactions, transformations and optional standardization (`standardize`).
* `mofia_score_test()`: model-based and robust score (Lagrange multiplier)
  tests with expected parameter changes, per parameter or per item.
* `mofia_scores()`: EAP factor scores with posterior SDs.
* `predict()` / `fitted()`: model-implied category probabilities and expected
  item scores for new moderator values; `confint()`; `tidy()` and `glance()`.
* `plot()`: factor means and SDs, moderated loadings and item response
  functions across a moderator.
* `mofia(starts = )`: random starts; warnings for sparse and unused categories;
  `mofia_collapse()` to merge sparse categories.

* First version: moderated nonlinear factor analysis by marginal maximum
  likelihood with adaptive Gauss-Hermite quadrature and analytic gradients in
  C++.
* Binary, ordinal (logit or probit), continuous and mixed items; moderation of
  loadings, thresholds/intercepts, residual variances, factor means, variances
  and (for two factors) correlations; any number of factors in principle (the
  quadrature cost grows exponentially; accuracy checked for up to three).
* ML, MLR and MLF standard errors; scaled likelihood-ratio tests (`anova()`,
  `mofia_lrt()`) and Wald tests (`mofia_wald()`).
* Simulation from fitted or specified models (`mofia_simulate()`).
* Example data `parabelief` (Hoogeveen et al., 2024).
