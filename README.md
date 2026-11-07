# MoFIA: Moderated Factor and Item Analysis

MoFIA is an open-source R package for **moderated nonlinear factor analysis
(MNLFA)**. In these factor models, item parameters (loadings, thresholds or
intercepts, residual variances) and factor parameters (means, variances,
correlations) can depend on any number of moderators, continuous or
categorical. Typical uses are testing measurement invariance and differential
item functioning (DIF) across age, gender or other characteristics without
forming groups, and comparing latent means and variances across moderators on
a common scale.

## Features

* Binary, ordinal (logit or probit link), continuous and mixed items.
* Moderation of loadings, thresholds/intercepts, residual variances, factor
  means, factor variances and factor correlations (Fisher-z or linear
  covariance).
* One or more factors (the cost of the quadrature grows exponentially with the
  number of factors; accuracy has been checked for up to three); variance or
  marker-item identification; missing item responses (full information,
  assuming missingness at random).
* Marginal maximum likelihood with adaptive Gauss-Hermite quadrature (exact for
  continuous items; checked against a finer rule and refined automatically
  otherwise) and analytic gradients in C++ (multithreaded with OpenMP where
  available).
* ML, robust (MLR, sandwich with the exact observed information) and
  score-based (MLF) standard errors.
* Scaled likelihood-ratio tests (Satorra-Bentler), Wald tests, and score
  (Lagrange multiplier) tests that screen all items for DIF from one fit.
* Model syntax (`"f1 =~ y1 + y2"`) and moderator formulas with factors,
  interactions and optional standardization.
* EAP factor scores, model-implied predictions, confidence intervals, broom-style
  `tidy()`/`glance()`, and plots of moderated parameters and item response
  functions.
* Safeguards: unidentified moderation patterns are refused, saddle points and
  singular information matrices are reported instead of standard errors,
  likelihood-ratio tests check that the models are comparable; random starts,
  warnings for sparse categories and a helper to collapse them.
* Simulation of data from any model.

## Installation

```r
# install.packages("remotes")
remotes::install_github("hrgodmann/MoFIA")
```

A C++ compiler is needed (Rtools on Windows, the Xcode command line tools on
macOS).

## Example

```r
library(MoFIA)
data(parabelief)

ub <- c("SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession",
        "Sorcery", "HauntedPlaces")
others <- setdiff(ub, "SpiritsGhosts")          # SpiritsGhosts as anchor item

scalar <- mofia(mofia_model(parabelief, items = ub, moderators = c("Age", "Male")))
metric <- mofia(mofia_model(parabelief, items = ub, moderators = c("Age", "Male"),
                            threshold_mod = others))

anova(scalar, metric)                            # threshold DIF?
mofia_score_test(scalar, by = "item")            # which items? (one fit)
mofia_wald(metric, "^beta\\[DevilPossession,", regex = TRUE)
summary(metric)
plot(metric, "items", moderator = "Male", items = "DevilPossession")
```

See `vignette("MoFIA")` for a walk-through.

## Data

The example data `parabelief` are a subset of the data of Hoogeveen, S.,
Borsboom, D., Kucharský, Š., Marsman, M., Molenaar, D., de Ron, J., Sekulovski,
N., Visser, I., van Elk, M., & Wagenmakers, E.-J. (2024). Prevalence, patterns
and predictors of paranormal beliefs in The Netherlands: A several-analysts
approach. *Royal Society Open Science, 11*(9), 240049.
<https://doi.org/10.1098/rsos.240049>. Please cite this article when you use
the data.

## References

* Bauer, D. J., & Hussong, A. M. (2009). Psychometric approaches for developing
  commensurate measures across independent studies: Traditional and new
  models. *Psychological Methods, 14*(2), 101–125.
  <https://doi.org/10.1037/a0015583>
* Bauer, D. J. (2017). A more general model for testing measurement invariance
  and differential item functioning. *Psychological Methods, 22*(3), 507–526.
  <https://doi.org/10.1037/met0000077>

## License

GPL (>= 3)
