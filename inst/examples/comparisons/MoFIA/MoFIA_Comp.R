
############### MoFIA runs on real data - read in data sets

############### Prepare package
# remotes::install_github("hrgodmann/MoFIA@2553b5b")

library(MoFIA)

# Data Set 1 -- Two factors
ParaNormal <- read.csv("/Users/henrikgodmann/Desktop/GitHub/MoFIA/inst/examples/comparisons/datasets/paranormal_beliefs_all_cleanedFull.csv")
head(ParaNormal)
# remove cases without any observed item response
ParaNormal <- ParaNormal[rowSums(is.na(ParaNormal[,c(1:14)])) < 14,]

# Factor Structure
# Factor EHB
# Item 1: HealingHands
# Item 2: LifeCoursePossessions
# Item 3: Graphology
# Item 4: Dowsing
# Item 5: MindReading
# Item 6: Willpower
# Item 7: SeeingFuture
# Item 8: SeeingEvents
# Factor UB
# Item 9: SpiritsGhosts
# Item 10: GnomesElves
# Item 11: Aliens
# Item 12: DevilPossession
# Item 13: Sorcery
# Item 14: HauntedPlaces

# Fit -- Moderators: Age (continuous), Male (binary); Items: Ordinal
# Anchor Items: HealingHands and SpiritsGhosts

EHB_items <- c("HealingHands", "LifeCoursePossessions", 
            "Graphology", "Dowsing", "MindReading", 
            "Willpower", "SeeingFuture", "SeeingEvents")
UB_items <- c("SpiritsGhosts", "GnomesElves", "Aliens", 
            "DevilPossession", "Sorcery", "HauntedPlaces")
anchor_items <- c("HealingHands", "SpiritsGhosts")

Model1_conf <- mofia_model(ParaNormal,
                      loadings = list(EHB = EHB_items,
                      UB = UB_items),
                      moderators = c("Age","Male"),
                      loading_mod = c(setdiff(EHB_items, anchor_items), 
                                      setdiff(UB_items, anchor_items)),
                      threshold_mod = c(setdiff(EHB_items, anchor_items), 
                                        setdiff(UB_items, anchor_items)),
                      # ordinal items do not have residual variances:
                      # resvar_mod = c(setdiff(EHB_items, anchor_items), 
                      #                setdiff(UB_items, anchor_items)) 
                      mean_mod = TRUE, var_mod = TRUE, cor_mod = TRUE)
set.seed(42); fit_Model1_conf <- mofia(Model1_conf, starts = 10)



Model1_metric <- mofia_model(ParaNormal,
                           loadings = list(EHB = EHB_items,
                                           UB = UB_items),
                           moderators = c("Age","Male"),
                           # loading_mod = c(setdiff(EHB_items, anchor_items), 
                                           # setdiff(UB_items, anchor_items)),
                           threshold_mod = c(setdiff(EHB_items, anchor_items), 
                                             setdiff(UB_items, anchor_items)),
                           # ordinal items do not have residual variances:
                           # resvar_mod = c(setdiff(EHB_items, anchor_items), 
                           #                setdiff(UB_items, anchor_items)) 
                           mean_mod = TRUE, var_mod = TRUE, cor_mod = TRUE)
set.seed(42); fit_Model1_metric <- mofia(Model1_metric, starts = 10)



Model1_scalar <- mofia_model(ParaNormal,
                             loadings = list(EHB = EHB_items,
                                             UB = UB_items),
                             moderators = c("Age","Male"),
                             # loading_mod = c(setdiff(EHB_items, anchor_items), 
                             # setdiff(UB_items, anchor_items)),
                             # threshold_mod = c(setdiff(EHB_items, anchor_items), 
                                               # setdiff(UB_items, anchor_items)),
                             # ordinal items do not have residual variances:
                             # resvar_mod = c(setdiff(EHB_items, anchor_items), 
                             #                setdiff(UB_items, anchor_items)) 
                             mean_mod = TRUE, var_mod = TRUE, cor_mod = TRUE)
set.seed(42); fit_Model1_scalar <- mofia(Model1_scalar, starts = 10)

# ParaNormal (ordinal): matches Mplus MLR + ALGORITHM = INTEGRATION
summary(fit_Model1_conf,   se = "MLR", information = "reduced")
summary(fit_Model1_metric, se = "MLR", information = "reduced")
summary(fit_Model1_scalar, se = "MLR", information = "reduced")

# scaled LR tests (Satorra–Bentler, like Mplus MLR difference testing)
anova(fit_Model1_conf, fit_Model1_metric, fit_Model1_scalar, information = "reduced")

################################################################################
# Data Set 2 -- Two factors and a cross loading AND separately Three factors
HS39 <- lavaan::HolzingerSwineford1939
HS39 <- HS39[c(5,7:15)]
head(HS39)

# Factor Structure
# Factor Visual	(Fit 1 and 2)
# Item 1: x1 Visual perception
# Item 2: x2 Cubes
# Item 3: x3 Lozenges
# Item 9: x9 cross-loading (see Speed)
# Factor Textual	(Fit 2 only)	
# Item 4: x4 Paragraph comprehension
# Item 5: x5 Sentence completion
# Item 6: x6 Word meaning
# Factor Speed (Fit 1 and 2)
# Item 7: x7 Speeded addition
# Item 8: x8 Speeded counting of dots
# Item 9: x9 Speeded discrimination of straight vs. curved capitals

# Fit 1 -- Moderators: School (binary); Items: Continuous
# Anchor Items: x1 and x7

Vis_items <- c("x1", "x2", "x3","x9")
Spe_items <- c("x7", "x8", "x9")
anchor_items <- c("x1", "x7")

HS39$school <- ifelse(HS39$school == "Pasteur", 1, 0)

Model2_conf <- mofia_model(HS39,
                           loadings = list(Visual = Vis_items,
                                           Speed = Spe_items),
                           moderators = c("school"),
                           loading_mod = c(setdiff(Vis_items, anchor_items), 
                                           setdiff(Spe_items, anchor_items)),
                           threshold_mod = c(setdiff(Vis_items, anchor_items), 
                                             setdiff(Spe_items, anchor_items)),
                           resvar_mod = unique(c(Vis_items, anchor_items,
                                          Spe_items, anchor_items)),
                           continuous = unique(c(Vis_items, Spe_items)),
                           mean_mod = TRUE, var_mod = TRUE, cor_mod = TRUE)
set.seed(42); fit_Model2_conf <- mofia(Model2_conf, starts = 10)

Model2_metric <- mofia_model(HS39,
                             loadings = list(Visual = Vis_items,
                                             Speed = Spe_items),
                             moderators = c("school"),
                             # loading_mod = c(setdiff(Vis_items, anchor_items), 
                                             # setdiff(Spe_items, anchor_items)),
                             threshold_mod = c(setdiff(Vis_items, anchor_items), 
                                               setdiff(Spe_items, anchor_items)),
                             resvar_mod = unique(c(Vis_items, anchor_items,
                                            Spe_items, anchor_items)),
                             continuous = c(Vis_items, Spe_items),
                             mean_mod = TRUE, var_mod = TRUE, cor_mod = TRUE)
set.seed(42); fit_Model2_metric <- mofia(Model2_metric, starts = 10)

Model2_scalar <- mofia_model(HS39,
                             loadings = list(Visual = Vis_items,
                                             Speed = Spe_items),
                             moderators = c("school"),
                             # loading_mod = c(setdiff(Vis_items, anchor_items), 
                                             # setdiff(Spe_items, anchor_items)),
                             # threshold_mod = c(setdiff(Vis_items, anchor_items), 
                                               # setdiff(Spe_items, anchor_items)),
                             resvar_mod = unique(c(Vis_items, anchor_items,
                                            Spe_items, anchor_items)),
                             continuous = c(Vis_items, Spe_items),
                             mean_mod = TRUE, var_mod = TRUE, cor_mod = TRUE)
set.seed(42); fit_Model2_scalar <- mofia(Model2_scalar, starts = 10)

# HS39 (continuous): matches Mplus MLR + INFORMATION = OBSERVED
summary(fit_Model2_conf,   se = "MLR", information = "exact")
summary(fit_Model2_metric, se = "MLR", information = "exact")
summary(fit_Model2_scalar, se = "MLR", information = "exact")
anova(fit_Model2_conf, fit_Model2_metric, fit_Model2_scalar, information = "exact")


# Fit 2 -- Moderators: School (binary); Items: Continuous
# Anchor Items: x1 and x4 and x7

Tex_items <- c("x4", "x5", "x6")
anchor_items <- c("x1", "x4", "x7")

Model3_conf <- mofia_model(HS39,
                           loadings = list(Visual = Vis_items,
                                           Textual = Tex_items,
                                           Speed = Spe_items),
                           moderators = c("school"),
                           loading_mod = c(setdiff(Vis_items, anchor_items), 
                                           setdiff(Tex_items, anchor_items),
                                           setdiff(Spe_items, anchor_items)),
                           threshold_mod = c(setdiff(Vis_items, anchor_items), 
                                             setdiff(Tex_items, anchor_items),
                                             setdiff(Spe_items, anchor_items)),
                           resvar_mod = unique(c(Vis_items, anchor_items,
                                          Tex_items, anchor_items,
                                          Spe_items, anchor_items)),
                           continuous = c(Vis_items, Tex_items, Spe_items),
                           mean_mod = TRUE, var_mod = TRUE, cor_mod = FALSE)
set.seed(42); fit_Model3_conf <- mofia(Model3_conf, starts = 10)

Model3_metric <- mofia_model(HS39,
                             loadings = list(Visual = Vis_items,
                                             Textual = Tex_items,
                                             Speed = Spe_items),
                             moderators = c("school"),
                             # loading_mod = c(setdiff(Vis_items, anchor_items), 
                                             # setdiff(Tex_items, anchor_items),
                                             # setdiff(Spe_items, anchor_items)),
                             threshold_mod = c(setdiff(Vis_items, anchor_items), 
                                               setdiff(Tex_items, anchor_items),
                                               setdiff(Spe_items, anchor_items)),
                             resvar_mod = unique(c(Vis_items, anchor_items,
                                            Tex_items, anchor_items,
                                            Spe_items, anchor_items)),
                             continuous = c(Vis_items, Tex_items, Spe_items),
                             mean_mod = TRUE, var_mod = TRUE, cor_mod = FALSE)
set.seed(42); fit_Model3_metric <- mofia(Model3_metric, starts = 10)

Model3_scalar <- mofia_model(HS39,
                             loadings = list(Visual = Vis_items,
                                             Textual = Tex_items,
                                             Speed = Spe_items),
                             moderators = c("school"),
                             # loading_mod = c(setdiff(Vis_items, anchor_items), 
                                             # setdiff(Tex_items, anchor_items),
                                             # setdiff(Spe_items, anchor_items)),
                             # threshold_mod = c(setdiff(Vis_items, anchor_items), 
                                               # setdiff(Tex_items, anchor_items),
                                               # setdiff(Spe_items, anchor_items)),
                             resvar_mod = unique(c(Vis_items, anchor_items,
                                            Tex_items, anchor_items,
                                            Spe_items, anchor_items)),
                             continuous = c(Vis_items, Tex_items, Spe_items),
                             mean_mod = TRUE, var_mod = TRUE, cor_mod = FALSE)
set.seed(42); fit_Model3_scalar <- mofia(Model3_scalar, starts = 10)

# HS39 (continuous): matches Mplus MLR + INFORMATION = OBSERVED
summary(fit_Model3_conf,   se = "MLR", information = "exact")
summary(fit_Model3_metric, se = "MLR", information = "exact")
summary(fit_Model3_scalar, se = "MLR", information = "exact")
anova(fit_Model3_conf, fit_Model3_metric, fit_Model3_scalar, information = "exact")


################################################################################
# Data Set 3 -- Two factors

data(BFI2, package = "ShinyItemAnalysis")
BFI2 <- BFI2[c(1:62)]
head(BFI2)

Anx_items <- c("i4", "i19", "i34", "i49")
Dep_items <- c("i9", "i24", "i39", "i54")
BFI2 <- BFI2[c(Anx_items, Dep_items, "Gender", "Age")]
BFI2$Age <- as.numeric(scale(BFI2$Age))
anchor_items <- c("i4", "i9")

# Factor Structure
# Factor Anxiety
# Item 1: i4	is relaxed, handles stress well (R)
# Item 2: i19	can be tense
# Item 3: i34	worries a lot
# Item 4: i49	rarely feels anxious or afraid (R)
# Factor Depression
# Item 5: i9	stays optimistic after experiencing a setback (R)
# Item 6: i24	feels secure, comfortable with self (R)
# Item 7: i39	often feels sad
# Item 8: i54	tends to feel depressed, blue

# Fit 1 -- Moderators: Age (continuous), Gender (binary); Items: Ordinal
# Anchor Items: i4 and i9

Model4_conf <- mofia_model(BFI2,
                           loadings = list(Anxiety = Anx_items,
                                           Depression = Dep_items),
                           moderators = c("Age","Gender"),
                           loading_mod = c(setdiff(Anx_items, anchor_items), 
                                           setdiff(Dep_items, anchor_items)),
                           threshold_mod = c(setdiff(Anx_items, anchor_items), 
                                             setdiff(Dep_items, anchor_items)),
                           # ordinal items do not have residual variances:
                           # resvar_mod = c(setdiff(Anx_items, anchor_items), 
                           #                setdiff(Dep_items, anchor_items)) 
                           mean_mod = TRUE, var_mod = TRUE, cor_mod = TRUE)
set.seed(42); fit_Model4_conf <- mofia(Model4_conf, starts = 10)

Model4_metric <- mofia_model(BFI2,
                             loadings = list(Anxiety = Anx_items,
                                             Depression = Dep_items),
                             moderators = c("Age","Gender"),
                             # loading_mod = c(setdiff(Anx_items, anchor_items), 
                                             # setdiff(Dep_items, anchor_items)),
                             threshold_mod = c(setdiff(Anx_items, anchor_items), 
                                               setdiff(Dep_items, anchor_items)),
                             # ordinal items do not have residual variances:
                             # resvar_mod = c(setdiff(Anx_items, anchor_items), 
                             #                setdiff(Dep_items, anchor_items)) 
                             mean_mod = TRUE, var_mod = TRUE, cor_mod = TRUE)
set.seed(42); fit_Model4_metric <- mofia(Model4_metric, starts = 10)

Model4_scalar <- mofia_model(BFI2,
                             loadings = list(Anxiety = Anx_items,
                                             Depression = Dep_items),
                             moderators = c("Age","Gender"),
                             # loading_mod = c(setdiff(Anx_items, anchor_items), 
                                             # setdiff(Dep_items, anchor_items)),
                             # threshold_mod = c(setdiff(Anx_items, anchor_items), 
                                               # setdiff(Dep_items, anchor_items)),
                             # ordinal items do not have residual variances:
                             # resvar_mod = c(setdiff(Anx_items, anchor_items), 
                             #                setdiff(Dep_items, anchor_items)) 
                             mean_mod = TRUE, var_mod = TRUE, cor_mod = TRUE)
set.seed(42); fit_Model4_scalar <- mofia(Model4_scalar, starts = 10)

# BFI-2 (ordinal): matches Mplus MLR + ALGORITHM = INTEGRATION
summary(fit_Model4_conf,   se = "MLR", information = "reduced")
summary(fit_Model4_metric, se = "MLR", information = "reduced")
summary(fit_Model4_scalar, se = "MLR", information = "reduced")

# scaled LR tests (Satorra–Bentler, like Mplus MLR difference testing)
anova(fit_Model4_conf, fit_Model4_metric, fit_Model4_scalar, information = "reduced")


# save the R workspace with all results (read by ../compare_MoFIA_Mplus.R)
save.image(file = "/Users/henrikgodmann/Desktop/GitHub/MoFIA/inst/examples/comparisons/MoFIA/MoFIA_Comp_results.RData")


