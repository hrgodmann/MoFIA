engine_ll <- function(theta, model, link, engine, nq = 7) {
  nodes <- MoFIA:::.mofia_nodes(nq, model$D)
  ad <- MoFIA:::.mofia_adapt(theta, model, link, nodes)
  MoFIA:::.mofia_ll(theta, model, link, nodes, scores = TRUE, per_person = TRUE, adapt = ad, engine = engine)
}

test_that("C++ and R engines agree (ordinal and mixed, logit and probit)", {
  for (cont in c(FALSE, TRUE)) {
    s <- sim_two_factor(N = 150, continuous = cont)
    m <- do.call(mofia_model, c(list(data = s$data), s$args))
    th <- s$theta[m$par$name]
    for (link in c("logit", "probit")) {
      a <- engine_ll(th, m, link, "R"); b <- engine_ll(th, m, link, "cpp")
      expect_equal(as.numeric(a), as.numeric(b), tolerance = 1e-10)
      expect_equal(attr(a, "scores"), attr(b, "scores"), tolerance = 1e-8)
    }
  }
})

test_that("analytic gradient matches numerical derivatives", {
  skip_if_not_installed("numDeriv")
  s <- sim_two_factor(N = 120, continuous = TRUE)
  m <- do.call(mofia_model, c(list(data = s$data), s$args))
  th <- s$theta[m$par$name]
  for (link in c("logit", "probit")) {
    nodes <- MoFIA:::.mofia_nodes(7, 2)
    ad <- MoFIA:::.mofia_adapt(th, m, link, nodes)
    g <- colSums(attr(MoFIA:::.mofia_ll(th, m, link, nodes, scores = TRUE, adapt = ad), "scores"))
    gn <- numDeriv::grad(function(t) MoFIA:::.mofia_ll(t, m, link, nodes, adapt = ad), th)
    expect_equal(unname(g), gn, tolerance = 1e-5)
  }
})

test_that("binary items and models without moderators work", {
  data(parabelief)
  d <- parabelief; d$Aliens <- as.integer(d$Aliens > 1); d$Sorcery <- as.integer(d$Sorcery > 1)
  m <- mofia_model(d, items = c("Aliens", "Sorcery", "HauntedPlaces", "SpiritsGhosts"))
  expect_equal(m$K, c(2L, 2L, 5L, 5L))
  fit <- mofia(m)
  expect_true(fit$converged)
  expect_equal(fit$npar, 1 + 1 + 4 + 4 + 4)
})
