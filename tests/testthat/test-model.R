test_that("loading patterns and moderation inputs are equivalent", {
  data(parabelief)
  eha <- c("HealingHands", "LifeCoursePossessions", "Graphology")
  ub <- c("SpiritsGhosts", "GnomesElves", "Aliens")
  m1 <- mofia_model(parabelief, loadings = list(EHA = eha, UB = ub), moderators = c("Age", "Male"),
                    threshold_mod = "Aliens", loading_mod = list(Graphology = "Male"))
  L <- cbind(EHA = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), UB = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
  tm <- matrix(FALSE, 6, 2); tm[6, ] <- TRUE
  lm <- matrix(FALSE, 6, 2); lm[3, 2] <- TRUE
  m2 <- mofia_model(parabelief, items = c(eha, ub), moderators = c("Age", "Male"), loadings = L,
                    threshold_mod = tm, loading_mod = lm)
  expect_identical(m1$par, m2$par)
  expect_true(all(c("beta[Aliens,Age]", "beta[Aliens,Male]", "lambda1[Graphology,EHA,Male]",
                    "zcor0[EHA,UB]", "delta[UB,Male]") %in% m1$par$name))
  expect_false(any(grepl("zcor1", m1$par$name)))
  expect_output(print(m1), "Moderated factor model")
})

test_that("invalid specifications are refused", {
  data(parabelief)
  expect_error(mofia_model(parabelief, items = c("Aliens", "Nope"), moderators = "Age"), "Not in the data")
  expect_error(mofia_model(parabelief, items = c("Aliens", "Sorcery"), moderators = "Age",
                           threshold_mod = "Nope"), "Unknown item")
  three <- list(a = c("HealingHands", "Dowsing"), b = c("Aliens", "Sorcery"), c = c("Willpower", "MindReading"))
  expect_error(mofia_model(parabelief, loadings = three, moderators = "Age", cor_mod = TRUE), "not methodologically")
  expect_silent(mofia_model(parabelief, loadings = three, moderators = "Age", cor_mod = TRUE, experimental_cor = TRUE))
  expect_error(mofia_model(parabelief, items = c("Aliens", "Sorcery"), moderators = "Age",
                           identification = "marker", loading_mod = TRUE), "unmoderated loading")
})

test_that("cases with missing moderators are removed with a warning", {
  data(parabelief)
  d <- parabelief[1:50, ]; d$Age[1:3] <- NA
  expect_warning(m <- mofia_model(d, items = c("Aliens", "Sorcery", "HauntedPlaces"), moderators = "Age"), "3 case")
  expect_equal(m$N, 47)
})

test_that("mofia_drop removes moderation effects", {
  data(parabelief)
  m <- mofia_model(parabelief, items = c("Aliens", "Sorcery", "HauntedPlaces"), moderators = c("Age", "Male"),
                   threshold_mod = "Sorcery")
  m0 <- mofia_drop(m, c("beta[Sorcery,Age]", "gamma[eta,Male]"))
  expect_setequal(setdiff(m$par$name, m0$par$name), c("beta[Sorcery,Age]", "gamma[eta,Male]"))
  expect_error(mofia_drop(m, "tau[Sorcery,1]"), "Only moderation")
})
