# Builds data/parabelief.rda from the subset used in the MoFIA examples.
# Source: Hoogeveen et al. (2024), Royal Society Open Science, 11(9), 240049.
raw <- read.csv("data-raw/ParaBeliefReduced.csv")
items <- c("HealingHands", "LifeCoursePossessions", "Graphology", "Dowsing", "MindReading", "Willpower",
           "SeeingFuture", "SeeingEvents", "SpiritsGhosts", "GnomesElves", "Aliens", "DevilPossession",
           "Sorcery", "HauntedPlaces")
parabelief <- raw[, c(items, "Age", "Male")]
for (v in c(items, "Male")) parabelief[[v]] <- as.integer(parabelief[[v]])
parabelief$Age <- as.numeric(parabelief$Age)
save(parabelief, file = "data/parabelief.rda", compress = "xz")
