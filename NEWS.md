# MoFIA 0.1.0

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
  and (for two factors) correlations; any number of factors.
* ML, MLR and MLF standard errors; scaled likelihood-ratio tests (`anova()`,
  `mofia_lrt()`) and Wald tests (`mofia_wald()`).
* Simulation from fitted or specified models (`mofia_simulate()`).
* Example data `parabelief` (Hoogeveen et al., 2024).
