#' Paranormal beliefs in the Netherlands (subset)
#'
#' Fourteen paranormal-belief items and two covariates for 600 respondents, a
#' subset of the nationally representative Dutch sample (N = 2,534) of
#' Hoogeveen et al. (2024). The subset (300 men, 300 women) was selected to
#' illustrate moderated factor analysis; results do not generalize to the
#' Dutch population.
#'
#' All items ask how much the respondent believes in a phenomenon, rated from
#' 1 (*definitely not*) to 5 (*definitely*). The first eight items measure
#' belief in *extraordinary human abilities* (EHA), the last six belief in
#' *unearthly beings* (UB).
#'
#' @format A data frame with 600 rows and 16 variables:
#' \describe{
#'   \item{HealingHands}{Healing by laying on of hands (EHA).}
#'   \item{LifeCoursePossessions}{Life course description by means of possessions (EHA).}
#'   \item{Graphology}{Graphology, i.e. handwriting analysis (EHA).}
#'   \item{Dowsing}{Dowsing, i.e. locating water or minerals with a rod (EHA).}
#'   \item{MindReading}{Mind-reading (EHA).}
#'   \item{Willpower}{Causing events by willpower (EHA).}
#'   \item{SeeingFuture}{Seeing into the future (EHA).}
#'   \item{SeeingEvents}{Seeing events happening elsewhere (EHA).}
#'   \item{SpiritsGhosts}{Spirits and ghosts (UB).}
#'   \item{GnomesElves}{Gnomes and elves (UB).}
#'   \item{Aliens}{Extraterrestrial beings (UB).}
#'   \item{DevilPossession}{Devil possession (UB).}
#'   \item{Sorcery}{Sorcery or black magic (UB).}
#'   \item{HauntedPlaces}{Haunted places (UB).}
#'   \item{Age}{Age, standardized (z-scores).}
#'   \item{Male}{Gender: 1 = male, 0 = female.}
#' }
#' @source Hoogeveen, S., Borsboom, D., Kucharský, Š., Marsman, M., Molenaar,
#'   D., de Ron, J., Sekulovski, N., Visser, I., van Elk, M., & Wagenmakers,
#'   E.-J. (2024). Prevalence, patterns and predictors of paranormal beliefs in
#'   The Netherlands: A several-analysts approach. *Royal Society Open Science,
#'   11*(9), 240049. \doi{10.1098/rsos.240049}
#'
#'   The article is published under the Creative Commons Attribution License
#'   (CC BY 4.0). Please cite it when you use these data.
#' @examples
#' data(parabelief)
#' str(parabelief)
#' table(parabelief$GnomesElves)
"parabelief"
