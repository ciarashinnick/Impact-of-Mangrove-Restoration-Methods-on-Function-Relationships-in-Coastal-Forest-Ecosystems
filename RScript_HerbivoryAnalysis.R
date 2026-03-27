################################################################################
# Analyses of rates of herbivory between three mangrove forest types
#
# Context:
#   - Herbivory is measured from leaf scans taken within quadrats at two sites.
#   - Multiple scans/leaves belong to a single “leaf basket”, and baskets map ~1:1
#     to TREE_NUMBER (tree-level unit of replication).
#   - Goal: fit a Beta GLMM to herbivory proportion, testing for
#     differences among QUADRAT_TYPE while accounting for PLANT_SPECIES and QUADRAT.
#
# Workflow outline:
#   A) Housekeeping + package setup
#   B) Load raw scan-level data
#   C) Aggregate to tree/basket-level (one row per tree)
#   D) Subset Likupang, explore distribution
#   E) Handle zeros for Beta regression (boundary adjustment)
#   F) Spatial autocorrelation test (Moran’s I on quadrat-mean residuals)
#   G) Fit candidate Beta GLMMs, improve dispersion model, select via AIC
#   H) Final model refit with clear reference levels
#   I) Export results (tables, diagnostics) + publication figures
#   J) Subset Tiwoho, explore distribution
#   K) Handle zeros for Beta regression (boundary adjustment)
#   L) Spatial autocorrelation test (Moran’s I on quadrat-mean residuals)
#   M) Fit candidate Beta GLMMs, improve dispersion model, select via AIC
#   N) Final model refit with clear reference levels
#   O) Export results (tables, diagnostics) + publication figures
################################################################################



################################################################################
# A) Housekeeping + package setup
################################################################################

# remove everything currently held in the R memory
rm(list=ls())

# close all open graphics windows 
graphics.off() 

install.packages("spdep")        # spatial dependence tools incl. Moran’s I
install.packages("DHARMa")       # simulation-based residual diagnostics for GLMMs
install.packages("dplyr")        # data manipulation
install.packages("ggplot2")      # plotting
install.packages("glmmTMB")      # GLMM fitting incl. beta + dispersion models
install.packages("broom.mixed")  # tidy extraction of mixed model coefficients
install.packages("emmeans")      # estimated marginal means (model-adjusted means)
install.packages("sjPlot")       # model summary tables for reporting
install.packages("performance")  # model performance metrics incl. R2 for GLMMs



#################################################################################
# B) Load raw dataset (scan-level)
#################################################################################

# Read in the raw CSV. file.choose() opens a dialog to select the file.
dat0 <- read.csv(file.choose())  



################################################################################
# C) Aggregate scan-level measurements to tree/basket level (the scans are not 
  #independent if multiple scans/leaves come from the same basket/tree. The “unit
  #of replication” is one basket per tree. Aggregating to TREE_NUMBER gives one 
  #row per independent sampling unit)
################################################################################

library(dplyr)

dat_tree <- dat0 %>%
  mutate(
    # factors = categorical predictors for modelling / plotting
    SITE         = factor(SITE),
    QUADRAT      = factor(QUADRAT),
    QUADRAT_TYPE = factor(QUADRAT_TYPE),
    PLANT_SPECIES = factor(PLANT_SPECIES),
    TREE_NUMBER  = factor(TREE_NUMBER)
  ) %>%
  group_by(TREE_NUMBER) %>%
  summarise(
  #Aggregate scan-level areas to tree-level totals
  #na.rm=TRUE ensures missing scans don't break sums
    TOTAL_RESTORED_LEAF_AREA_mm = sum(TOTAL_RESTORED_LEAF_AREA_mm, na.rm = TRUE),
    TOTAL_HERBIVORY_AREA_mm     = sum(TOTAL_HERBIVORY_AREA_mm,     na.rm = TRUE),
    
    
#Carry forward variables expected to be constant within treewith first() takeing
#the first observed value in the group
    SITE         = first(SITE),
    QUADRAT      = first(QUADRAT),
    QUADRAT_TYPE = first(QUADRAT_TYPE),
    longitude_x  = first(longitude_x),
    latitude_y   = first(latitude_y),
    PLANT_SPECIES = first(PLANT_SPECIES),
    DBH_cm       = first(DBH_cm),
    CaCo_INDEX   = first(CaCo_INDEX),
    

#QA checks: how many distinct values occurred within each tree group. If any of 
#these >1, then that “constant” variable varied within tree, suggesting a data 
#integrity issue or grouping issue.
    n_SITE         = n_distinct(SITE),
    n_QUADRAT      = n_distinct(QUADRAT),
    n_QUADRAT_TYPE = n_distinct(QUADRAT_TYPE),
    n_lon          = n_distinct(longitude_x),
    n_lat          = n_distinct(latitude_y),
    n_SPECIES      = n_distinct(PLANT_SPECIES),
    n_DBH          = n_distinct(DBH_cm),
    n_CaCo         = n_distinct(CaCo_INDEX),
    

#drop grouping after summarise so downstream code works on a normal data frame 
#and derive herbivory metrics at tree level using the formulae:   
  #proportion = herbivory area / restored leaf area
  #percentage = 100 * proportion
    .groups = "drop"
  ) %>%
  mutate(
    HERBIVORY_PROPORTION = TOTAL_HERBIVORY_AREA_mm / TOTAL_RESTORED_LEAF_AREA_mm,
    HERBIVORY_PERCENTAGE = 100 * HERBIVORY_PROPORTION
  )


# Identify any trees where a “constant” variable was not actually constant.
# Expectation: zero rows. Any rows returned should be investigated.
problems <- dat_tree %>%
  filter(
    n_SITE != 1 | n_QUADRAT != 1 | n_QUADRAT_TYPE != 1 |
      n_lon  != 1 | n_lat    != 1 | n_SPECIES     != 1 |
      n_DBH  != 1 | n_CaCo   != 1
  )

problems   #if 0 rows -> constants behaved as expected during aggregation


# Export aggregated tree-level dataset (useful for manual QA and for archiving)
write.csv(
  dat_tree,
  file = "Mangrove_Thesis_TreeLevel_2026.02.04.csv",
  row.names = FALSE
)



###############################################################################
# D) Likupang analyses (site-specific)
###############################################################################

# Subset to Likupang (your analysis focus is site = Likupang)
dat_likupang <- dat_tree %>%
  filter(SITE == "Likupang")
nrow(dat_likupang)


###Explore distribution of herbivory (skew check). If distribution is 
  #strongly right-skewed and bounded, Beta regression is generally appropriate 
  #for the proportion scale.

library(ggplot2)

ggplot(
  dat_likupang,
  aes(x = HERBIVORY_PERCENTAGE)
) +
  geom_histogram(bins = 100) +
  labs(x = "Herbivory (%)",
       y = "Frequency",
       title = "Distribution of leaf herbivory — Likupang (tree level)") +
  theme_classic()



###############################################################################
# E) Deal with zeros prior to Beta GLMM (boundary adjustment)
# Why:
#   - Beta distributions are defined on (0,1), so exact 0 or 1 values cannot be used.
#   - With few zeros, a standard approach is the Smithson–Verkuilen transformation:
#       y* = (y*(n-1) + 0.5)/n
#     which pushes 0 slightly above 0 and 1 slightly below 1.
###############################################################################

# Count number of exact zeros on the proportion scale (Likupang, tree level)
sum(dat_likupang$HERBIVORY_PROPORTION == 0) #two absolute zeros for this dataset

# Apply Smithson–Verkuilen boundary adjustment
n <- nrow(dat_likupang)

dat_likupang <- dat_likupang %>%
  mutate(
    herb_beta = (HERBIVORY_PROPORTION * (n-1) + 0.5) / n
  )



################################################################################
# F) Spatial autocorrelation test (Moran’s I at quadrat scale)
# Why:
#   - Quadrats are spatially located; if residuals are spatially autocorrelated,
#     standard errors and inference may be biased.
# How:
#   1) Fit a baseline GLMM (simple structure)
#   2) Extract DHARMa scaled residuals
#   3) Average residuals per quadrat (since coords are quadrat-level)
#   4) Build a neighbour graph (k-nearest neighbours)
#   5) Run Moran’s I test on quadrat-mean residuals
################################################################################

library(glmmTMB)

# Baseline Beta GLMM, Fixed effect: QUADRAT_TYPE, Random effects: PLANT_SPECIES 
#and QUADRAT
m0 <- glmmTMB(
  herb_beta ~ QUADRAT_TYPE +
    (1 | PLANT_SPECIES) +
    (1 | QUADRAT),
  family = beta_family(),
  data = dat_likupang
)

# DHARMa residual simulation (n=2000 gives stable tests)
library(DHARMa)
res0 <- simulateResiduals(m0, n = 2000)
plot(res0)
# Add scaled residuals back to the dataset (one residual per tree)
dat_likupang$resid <- residuals(res0, type = "scaled")

# Collapse to quadrat mean residuals + mean coordinates
quad_resid <- dat_likupang %>%
  group_by(QUADRAT) %>%
  summarise(
    lon = mean(longitude_x, na.rm = TRUE),
    lat = mean(latitude_y, na.rm = TRUE),
    resid = mean(resid),
    .groups = "drop"
  )


print(quad_resid)


# Moran’s I test

library(spdep)
coords <- as.matrix(quad_resid[, c("lon", "lat")])

# Define neighbour structure using 4 nearest neighbours per quadrat
nb <- knn2nb(knearneigh(coords, k = 4))

# Convert neighbour object to weights (row-standardised)
lw <- nb2listw(nb, style = "W")

# Moran’s I: tests for positive spatial autocorrelation in residuals
moran.test(quad_resid$resid, lw)  ### p > 0.05 so no spatial terms needed



################################################################################
# G) Model selection and refinement (Beta GLMM with dispersion submodel)
# Key modelling choices:
#   - Response: herb_beta (boundary-adjusted proportion)
#   - Random effect: QUADRAT (sampling structure; trees nested in quadrats)
#   - Fixed effects (a priori): QUADRAT_TYPE, PLANT_SPECIES
#   - Candidate covariates: leaf area, DBH, CaCo (scaled)
#   - Dispersion (precision) model: allow dispersion to vary by QUADRAT_TYPE and CaCo
# Diagnostics:
#   - DHARMa tests: dispersion, uniformity, outliers
#   - AIC used for comparing candidate fixed-effect structures (same data rows)
################################################################################

library(dplyr)
library(glmmTMB)
library(DHARMa)
library(performance)
library(emmeans)

# Build modelling data frame with scaled predictors
dat <- dat_likupang %>%
  mutate(
    QUADRAT = factor(QUADRAT),
    QUADRAT_TYPE = factor(QUADRAT_TYPE),
    PLANT_SPECIES = factor(PLANT_SPECIES),
    
    # Scale continuous predictors to mean 0, SD 1
    # This improves numerical stability and convergence for GLMMs
    leafarea_z = as.numeric(scale(TOTAL_RESTORED_LEAF_AREA_mm)),
    DBH_z      = as.numeric(scale(as.numeric(DBH_cm))),
    CaCo_z     = as.numeric(scale(CaCo_INDEX / 100))  # percent -> proportion -> scale
  )


### Base model: decide PLANT_SPECIES as fixed vs random effect

# Species as random intercept:
m_base_RE <- glmmTMB(
  herb_beta ~ QUADRAT_TYPE +
    (1 | QUADRAT) +
    (1 | PLANT_SPECIES),
  family = beta_family(link = "logit"),
  data = dat
)

# Species as fixed effect:
m_base_FE <- glmmTMB(
  herb_beta ~ QUADRAT_TYPE + PLANT_SPECIES +
    (1 | QUADRAT),
  family = beta_family(link = "logit"),
  data = dat
)

# Compare AIC (lower is better fit-adjusted-for-complexity)
AIC(m_base_RE, m_base_FE)   ### suggests plant as fixed is better

# Diagnose both to understand which misspecification is occurring
res_FE <- simulateResiduals(m_base_FE, n = 2000)  
plot(res_FE)
testDispersion(res_FE)
testUniformity(res_FE)
testOutliers(res_FE)

res_RE <- simulateResiduals(m_base_RE, n = 2000)
plot(res_RE)
testDispersion(res_RE)
testUniformity(res_RE)
testOutliers(res_RE)

# Decision rationale:
#   - Proceed with plant as fixed effect, then improve dispersion via dispformula.


###Improve dispersion: allow dispersion to differ by QUADRAT_TYPE

m_base_FE_disp <- glmmTMB(
  herb_beta ~ QUADRAT_TYPE + PLANT_SPECIES +
    (1 | QUADRAT),
  dispformula = ~ QUADRAT_TYPE,
  family = beta_family(),
  data = dat
)

# Diagnostic checks after adding dispersion structure
res_FE_disp <- simulateResiduals(m_base_FE_disp, n = 2000)
plot(res_FE_disp)
testDispersion(res_FE_disp)
testUniformity(res_FE_disp)
testOutliers(res_FE_disp)


###Add candidate covariates and compare via AIC
# Important:
#   - AIC comparisons must use models fit to the SAME ROWS.
#   - Therefore, filter to complete cases for the candidate covariates first.

dat_sel <- dat %>%
  filter(
    !is.na(leafarea_z),
    !is.na(DBH_z),
    !is.na(CaCo_z)
  )

# Refit baseline dispersion model on the complete-case dataset
m_base_FE_disp2 <- update(m_base_FE_disp, data = dat_sel)

# Add predictors one-at-a-time, then the full model
m_FE_disp_leaf <- update(m_base_FE_disp, . ~ . + leafarea_z)
m_FE_disp_leaf2 <- update(m_FE_disp_leaf, data = dat_sel)

m_FE_disp_dbh <- update(m_base_FE_disp, . ~ . + DBH_z)
m_FE_disp_dbh2  <- update(m_FE_disp_dbh,  data = dat_sel)

m_FE_disp_caco <- update(m_base_FE_disp, . ~ . + CaCo_z)
m_FE_disp_caco2 <- update(m_FE_disp_caco, data = dat_sel)

m_FE_disp_full <- update(m_base_FE_disp,. ~ . + leafarea_z + DBH_z + CaCo_z)
m_FE_disp_full2 <- update(m_FE_disp_full, data = dat_sel)

# AIC comparison across candidate models (full model lowest AIC in the run)
AIC(
  m_base_FE_disp2,
  m_FE_disp_leaf2,
  m_FE_disp_dbh2,
  m_FE_disp_caco2,
  m_FE_disp_full2
)

# Choose best mean-structure model
m_final <- m_FE_disp_full2

# Diagnose the candidate final model
res_final <- simulateResiduals(m_final, n = 2000)
plot(res_final)
testDispersion(res_final)  #still a small dispersion deviation
testUniformity(res_final)
testOutliers(res_final)


# Final refinement:
#   - Add CaCo_z into dispersion model (dispformula)
#   - This can address heteroskedasticity that depends on canopy cover / quadrat context
m_final3 <- update(
  m_final,
  dispformula = ~ QUADRAT_TYPE + CaCo_z
)

# Diagnose the updated model
res3 <- simulateResiduals(m_final3, n = 2000)
plot(res3)
testDispersion(res3)
testUniformity(res3)
testOutliers(res3)

# At this stage: diagnostics indicate an acceptable fit in the workflow.

################################################################################
#H)Final model refit with explicit reference levels (interpretation)
# Why:
#   - Releveling factors changes which category is the baseline comparison.
#   - It does not change fitted values, but makes coefficient interpretation
#     aligned with your narrative (reference forest baseline, common species baseline).
################################################################################

dat_sel <- dat_sel %>%
  mutate(
    PLANT_SPECIES = relevel(PLANT_SPECIES, ref = "Rhizophora apiculata"),
    QUADRAT_TYPE  = relevel(QUADRAT_TYPE,  ref = "Reference Forest")
  )

# Refit the final model using the relevelled data frame
m_final3 <- glmmTMB(
  herb_beta ~ QUADRAT_TYPE + PLANT_SPECIES +
    leafarea_z + DBH_z + CaCo_z +
    (1 | QUADRAT),
  dispformula = ~ QUADRAT_TYPE + CaCo_z,
  family = beta_family(link = "logit"),
  data = dat_sel
)

# Print model summary (coefficients on link scale + random effect variance)
summary(m_final3)

# Create a formatted summary table
library(sjPlot)
tab_model(m_final3, show.est = TRUE, show.se = TRUE, show.stat = TRUE)

# Save model object for reproducibility (reload with readRDS() later)
saveRDS(m_final3, "Likupang_tree_betaGLMM_final.rds")

################################################################################
# I) Export results in publishable/reusable formats
################################################################################

#Export fixed-effect coefficients (tidy format, with Wald confidence intervals)
library(broom.mixed)

coef_tab <- tidy(
  m_final3,
  effects = "fixed",
  conf.int = TRUE,
  conf.method = "Wald"
)
coef_tab
write.csv(coef_tab,
          "Likupang_betaGLMM_coefficients.csv",
          row.names = FALSE)


###Compute R2 for mixed models (Nakagawa marginal + conditional)
#     - Marginal R2: variance explained by fixed effects
#     - Conditional R2: variance explained by fixed + random effects
library(performance)
r2_vals <- r2_nakagawa(m_final3)
r2_vals


###Estimated marginal means (predicted means adjusted for covariates)
#     - Useful for clear reporting of differences among QUADRAT_TYPE and species
library(emmeans)

emm_qt <- emmeans(m_final3, ~ QUADRAT_TYPE, type = "response")
emm_sp <- emmeans(m_final3, ~ PLANT_SPECIES, type = "response")
emm_qt
emm_sp

write.csv(as.data.frame(emm_qt),
          "Likupang_betaGLMM_emmeans_quadrat_type.csv",
          row.names = FALSE)

write.csv(as.data.frame(emm_sp),
          "Likupang_betaGLMM_emmeans_species.csv",
          row.names = FALSE)


###Export DHARMa diagnostic plots (for supplement / documentation)
library(DHARMa)

res_final <- simulateResiduals(m_final3, n = 2000)
plot(res_final)

pdf("Likupang_betaGLMM_DHARMa_diagnostics.pdf", width = 8, height = 8)
plot(res_final)
dev.off()


###Figures
# Figure types included:
#   1) Model-predicted means with confidence intervals (emmeans)
#   2) Raw data distribution plots (boxplot + jitter)
#   3) Continuous predictor response curves from prediction grids (predict.glmmTMB)

library(ggplot2)


###Quadrat/forest type effect figures

##Model-predicted differences among QUADRAT_TYPE (emmeans, response scale)
emm_qt_df <- as.data.frame(emm_qt)

p1a <- ggplot(emm_qt_df,
              aes(x = QUADRAT_TYPE,
                  y = response)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.15
  ) +
  ylab("Predicted herbivory proportion") +
  xlab("Quadrat type") +
  labs(title = "Likupang") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )
p1a

ggsave("Fig1a_predicted_herbivory_quadrat_type.png",
       p1a, width = 6, height = 4, dpi = 600)

##Raw distribution by QUADRAT_TYPE (boxplot + jitter) with consistent colours
p1b <- ggplot(dat_sel,
              aes(x = QUADRAT_TYPE,
                  y = HERBIVORY_PROPORTION,
                  fill = QUADRAT_TYPE,
                  colour = QUADRAT_TYPE)) +
  geom_boxplot(outlier.shape = NA,
               alpha = 0.6,
               width = 0.65) +
  geom_jitter(width = 0.15,
              alpha = 0.4,
              size = 1.5) +
  ylab("Observed herbivory proportion") +
  xlab("Quadrat type") +
  labs(title = "Likupang") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, family = "Times New Roman", face = "bold"),
    text = element_text(family = "Times New Roman")
  )

p1b

ggsave("Fig1b_herbivory_quadrat_type.png",
       p1a, width = 6, height = 4, dpi = 600)


###Plant species effect figures

##Model-predicted differences among PLANT_SPECIES (emmeans)
emm_sp_df <- as.data.frame(emm_sp)

p2a <- ggplot(emm_sp_df,
              aes(x = PLANT_SPECIES,
                  y = response)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.15
  ) +
  ylab("Predicted herbivory proportion") +
  xlab("Plant species") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               vjust = 1),
    axis.title.x = element_text(margin = margin(t = 12))
  )

p2a

##Raw distribution by PLANT_SPECIES (boxplot + jitter)
p2b <- ggplot(dat_sel,
              aes(x = PLANT_SPECIES,
                  y = HERBIVORY_PROPORTION)) +
  geom_boxplot(outlier.shape = NA,
               fill = "grey90",
               colour = "black",
               width = 0.6) +
  geom_jitter(width = 0.15,
              alpha = 0.5,
              size = 1.8) +
  ylab("Observed herbivory proportion") +
  xlab("Plant species") +
  labs(title = "Likupang") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               vjust = 1),
    axis.title.x = element_text(margin = margin(t = 12)),
    plot.title = element_text(hjust = 0.5, family = "Times New Roman", face = "bold"),
    text = element_text(family = "Times New Roman")
  )

p2b

ggsave("Fig2b_herbivory_plant_species.png",
       p2a, width = 6, height = 4, dpi = 600)


###Continuous predictor effects using prediction grids (CaCo and DBH)
# Why a prediction grid?
#   - One predictor is caried across across its observed range (e.g. CaCo_z),
#     while holding other predictors at typical values (0 for scaled covariates).
#   - We then predict fitted means and approximate 95% CIs using se.fit.
# Important detail:
#   - CaCo_z and DBH_z are scaled. We back-transform to original units for plotting.


###Plot the impact of vegetation density (CaCo)

##Recover scaling parameters explicitly:
    #CaCo_z was created as scale(CaCo_INDEX/100), so we reconstruct mean and sd
caco_raw_prop <- dat_sel$CaCo_INDEX / 100
caco_mean <- mean(caco_raw_prop, na.rm = TRUE)
caco_sd   <- sd(caco_raw_prop,   na.rm = TRUE)


###Build prediction grid:
#  - CaCo_z varies across observed range
#  - leafarea_z and DBH_z set to 0 (mean of scaled variables)
#  - QUADRAT_TYPE varies to produce separate lines by forest type
#  - PLANT_SPECIES fixed to the first factor level (reference for this plot)
#  - QUADRAT = NA so predictions are at population level (no specific quadrat)

newdat <- expand.grid(
  CaCo_z = seq(min(dat_sel$CaCo_z, na.rm = TRUE),
               max(dat_sel$CaCo_z, na.rm = TRUE),
               length.out = 100),
  leafarea_z   = 0,
  DBH_z        = 0,
  QUADRAT_TYPE = levels(dat_sel$QUADRAT_TYPE),
  PLANT_SPECIES = levels(dat_sel$PLANT_SPECIES)[1],
  QUADRAT      = NA
)


##Predict fitted values and SE on the response scale
pred <- predict(m_final3, newdata = newdat, type = "response", se.fit = TRUE)


##Attach predictions + compute approximate 95% CI
#Also back-transform CaCo_z to CaCo_INDEX (%) for interpretability
newdat <- newdat %>%
  mutate(
    fit = pred$fit,
    se  = pred$se.fit,
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se,
    CaCo_INDEX_pct = 100 * (CaCo_z * caco_sd + caco_mean)
  )


##Plot fitted lines + CI ribbons, colour by QUADRAT_TYPE
p3 <- ggplot(newdat,
             aes(x = CaCo_INDEX_pct,
                 y = fit,
                 colour = QUADRAT_TYPE,
                 fill = QUADRAT_TYPE)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              alpha = 0.2,
              colour = NA) +
  xlab("Caco - Canopy cover index (%)") +
  ylab("Predicted herbivory proportion") +
  theme_classic(base_size = 14)

p3

##Save figure in high-resolution PNG and vector PDF
ggsave("Fig3_predicted_herbivory_CaCo.png",
       p3, width = 6.5, height = 4.5, dpi = 600)
ggsave("Fig3_predicted_herbivory_CaCo.pdf",
       p3, width = 6.5, height = 4.5)



##canopy cover by forest type 

library(dplyr)

caco_by_type <- dat_sel %>%
  group_by(QUADRAT_TYPE) %>%
  summarise(
    mean_caco = mean(CaCo_INDEX, na.rm = TRUE),
    sd_caco = sd(CaCo_INDEX, na.rm = TRUE),
    n = n()
  )

caco_by_type
ggplot(dat_sel, aes(x = QUADRAT_TYPE, y = CaCo_INDEX)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  ylab("Canopy cover index (%)") +
  xlab("Forest type") +
  theme_classic()


###Plot the impact of tree size (DBH)

##Recover scaling parameters for DBH_z = scale(DBH_cm)
dbh_mean <- mean(dat_sel$DBH_cm, na.rm = TRUE)
dbh_sd   <- sd(dat_sel$DBH_cm,   na.rm = TRUE)

#Prediction grid:
#  - DBH_z varies; CaCo_z fixed at 0 (mean)
#  - leafarea_z fixed at 0 (mean)
#  - QUADRAT_TYPE varies (separate lines)
#  - species fixed to reference level for this plot
newdat_dbh <- expand.grid(
  DBH_z = seq(min(dat_sel$DBH_z, na.rm = TRUE),
              max(dat_sel$DBH_z, na.rm = TRUE),
              length.out = 100),
  leafarea_z   = 0,
  CaCo_z       = 0,
  QUADRAT_TYPE = levels(dat_sel$QUADRAT_TYPE),
  PLANT_SPECIES = levels(dat_sel$PLANT_SPECIES)[1],
  QUADRAT      = NA
)


##Predictions + SE on response scale
pred_dbh <- predict(m_final3, newdata = newdat_dbh, type = "response", se.fit = TRUE)

##Attach predictions + CI; back-transform DBH_z to DBH in cm
newdat_dbh <- newdat_dbh %>%
  mutate(
    fit = pred_dbh$fit,
    se  = pred_dbh$se.fit,
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se,
    DBH_cm = (DBH_z * dbh_sd + dbh_mean)
  )


##Plot DBH effect with quadrat-type colours
p4 <- ggplot(newdat_dbh,
             aes(x = DBH_cm,
                 y = fit,
                 colour = QUADRAT_TYPE,
                 fill = QUADRAT_TYPE)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              alpha = 0.2,
              colour = NA) +
  xlab("DBH (cm)") +
  ylab("Predicted herbivory proportion") +
  theme_classic(base_size = 14)

p4

##Save figure outputs
ggsave("Fig4_predicted_herbivory_DBH.png",
       p4, width = 6.5, height = 4.5, dpi = 600)
ggsave("Fig4_predicted_herbivory_DBH.pdf",
       p4, width = 6.5, height = 4.5)


##DBH cover by forest type 

library(dplyr)

dbh_by_type <- dat_sel %>%
  group_by(QUADRAT_TYPE) %>%
  summarise(
    mean_dbh = mean(DBH_cm, na.rm = TRUE),
    sd_dbh = sd(DBH_cm, na.rm = TRUE),
    n = n()
  )

dbh_by_type
ggplot(dat_sel, aes(x = QUADRAT_TYPE, y = DBH_cm)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  ylab("DBH (cm)") +
  xlab("Forest type") +
  theme_classic()


###sjPlot figures (alternative output style)

library(sjPlot)
library(ggplot2)
library(dplyr)
library(emmeans)


# Define a re-usable publication theme for sjPlot-based ggplots
pub_theme <- function() {
  theme_classic(base_size = 14) +
    theme(
      legend.key.size = unit(0.9, "cm"),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text  = element_text(size = 11),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x  = element_text(size = 12, colour = "black"),
      axis.text.y  = element_text(size = 12, colour = "black"),
      strip.text   = element_text(size = 12, face = "bold", colour = "black")
    )
}

### CaCo effect plot (x-axis on 0–100%)

# Compute mean/sd for back-transforming z-scale to percent scale for labelling
caco_prop <- dat_sel$CaCo_INDEX / 100
caco_mean <- mean(caco_prop, na.rm = TRUE)
caco_sd   <- sd(caco_prop,   na.rm = TRUE)

# Create z-scale tick positions and map them to percent labels
z_breaks <- pretty(range(dat_sel$CaCo_z, na.rm = TRUE), n = 5)
pct_labels <- round(100 * (z_breaks * caco_sd + caco_mean), 1)

# sjPlot prediction plot (note: dispersion submodel may be ignored in plot_model)
p_caco <- plot_model(
  m_final3,
  type  = "pred",
  terms = c("CaCo_z [all]", "QUADRAT_TYPE")
) +
  labs(x = "Canopy cover index (%)",
       y = "Predicted herbivory proportion",
       colour = "Quadrat type",
       fill = "Quadrat type",
       title = "Likupang") +
  scale_x_continuous(breaks = z_breaks, labels = pct_labels) +
  pub_theme() +
  theme(
    plot.title = element_text(hjust = 0.5, family = "Times New Roman", face = "bold"),
    text = element_text(family = "Times New Roman")
  )

p_caco
ggsave("Fig_effect_CaCo_quadratType.png", p_caco, width = 7, height = 5, dpi = 600)
ggsave("Fig_effect_CaCo_quadratType.pdf", p_caco, width = 7, height = 5)


### DBH

# Compute mean/sd for back-transforming DBH_z ticks to DBH_cm
dbh_mean <- mean(dat_sel$DBH_cm, na.rm = TRUE)
dbh_sd   <- sd(dat_sel$DBH_cm,   na.rm = TRUE)

z_breaks_dbh <- pretty(range(dat_sel$DBH_z, na.rm = TRUE), n = 5)
cm_labels <- round(z_breaks_dbh * dbh_sd + dbh_mean, 1)

# sjPlot prediction plot for DBH
p_dbh <- plot_model(
  m_final3,
  type  = "pred",
  terms = c("DBH_z [all]", "QUADRAT_TYPE")
) +
  labs(x = "DBH (cm)",
       y = "Predicted herbivory proportion",
       colour = "Quadrat type",
       fill = "Quadrat type",
       title = "Likupang") +
  scale_x_continuous(breaks = z_breaks_dbh, labels = cm_labels) +
  pub_theme() +
  theme(
    plot.title = element_text(hjust = 0.5, family = "Times New Roman", face = "bold"),
    text = element_text(family = "Times New Roman")
  )

p_dbh
ggsave("Fig_effect_DBH_quadratType.png", p_dbh, width = 7, height = 5, dpi = 600)
ggsave("Fig_effect_DBH_quadratType.pdf", p_dbh, width = 7, height = 5)


#### quadrat

# sjPlot predicted means for quadrat type (categorical predictor)
p_qtype <- plot_model(
  m_final3,
  type = "pred",
  terms = "QUADRAT_TYPE"
) +
  labs(x = "Quadrat type",
       y = "Predicted herbivory proportion") +
  pub_theme() +
  theme(legend.position = "none") +
  ggtitle("")

p_qtype
ggsave("Fig_effect_QuadratType.png", p_qtype, width = 6, height = 4, dpi = 600)
ggsave("Fig_effect_QuadratType.pdf", p_qtype, width = 6, height = 4)


### plant species

# For many species levels, an emmeans dot + CI plot is usually clearer
emm_sp <- emmeans(m_final3, ~ PLANT_SPECIES, type = "response")
emm_sp_df <- as.data.frame(emm_sp)

p_species <- ggplot(emm_sp_df,
                    aes(x = reorder(PLANT_SPECIES, response),
                        y = response)) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  coord_flip() +
  labs(x = "Plant species",
       y = "Predicted herbivory proportion",
       title = "Likupang") +
  pub_theme() +
  theme(
    plot.title = element_text(hjust = 0.5, family = "Times New Roman", face = "bold"),
    text = element_text(family = "Times New Roman")
  )

p_species
ggsave("Fig_effect_PlantSpecies.png", p_species, width = 7.5, height = 6.5, dpi = 600)
ggsave("Fig_effect_PlantSpecies.pdf", p_species, width = 7.5, height = 6.5)




###############################################################################
# J) Tiwoho analyses (site-specific)
###############################################################################

# Subset to Tiwoho (your analysis focus is site = Tiwoho)
dat_tiwoho <- dat_tree %>%
  filter(SITE == "Tiwoho")
nrow(dat_tiwoho)


###Explore distribution of herbivory (skew check). If distribution is 
#strongly right-skewed and bounded, Beta regression is generally appropriate 
#for the proportion scale.

library(ggplot2)

ggplot(
  dat_tiwoho,
  aes(x = HERBIVORY_PERCENTAGE)
) +
  geom_histogram(bins = 100) +
  labs(x = "Herbivory (%)",
       y = "Frequency",
       title = "Distribution of leaf herbivory — Tiwoho (tree level)") +
  theme_classic()



###############################################################################
# K) Deal with zeros prior to Beta GLMM (boundary adjustment)
# Why:
#   - Beta distributions are defined on (0,1), so exact 0 or 1 values cannot be used.
#   - With few zeros, a standard approach is the Smithson–Verkuilen transformation:
#       y* = (y*(n-1) + 0.5)/n
#     which pushes 0 slightly above 0 and 1 slightly below 1.
###############################################################################

# Count number of exact zeros on the proportion scale (Tiwoho, tree level)
sum(dat_tiwoho$HERBIVORY_PROPORTION == 0) #two absolute zeros for this dataset

# Apply Smithson–Verkuilen boundary adjustment
n <- nrow(dat_tiwoho)

dat_tiwoho <- dat_tiwoho %>%
  mutate(
    herb_beta = (HERBIVORY_PROPORTION * (n-1) + 0.5) / n
  )



################################################################################
# L) Spatial autocorrelation test (Moran’s I at quadrat scale)
# Why:
#   - Quadrats are spatially located; if residuals are spatially autocorrelated,
#     standard errors and inference may be biased.
# How:
#   1) Fit a baseline GLMM (simple structure)
#   2) Extract DHARMa scaled residuals
#   3) Average residuals per quadrat (since coords are quadrat-level)
#   4) Build a neighbour graph (k-nearest neighbours)
#   5) Run Moran’s I test on quadrat-mean residuals
################################################################################

library(glmmTMB)

# Baseline Beta GLMM, Fixed effect: QUADRAT_TYPE, Random effects: PLANT_SPECIES 
#and QUADRAT
m0 <- glmmTMB(
  herb_beta ~ QUADRAT_TYPE +
    (1 | PLANT_SPECIES) +
    (1 | QUADRAT),
  family = beta_family(),
  data = dat_tiwoho
)

# DHARMa residual simulation (n=2000 gives stable tests)
library(DHARMa)
res0 <- simulateResiduals(m0, n = 2000)
plot(res0)
# Add scaled residuals back to the dataset (one residual per tree)
dat_tiwoho$resid <- residuals(res0, type = "scaled")

# Collapse to quadrat mean residuals + mean coordinates
quad_resid <- dat_tiwoho %>%
  group_by(QUADRAT) %>%
  summarise(
    lon = mean(longitude_x, na.rm = TRUE),
    lat = mean(latitude_y, na.rm = TRUE),
    resid = mean(resid),
    .groups = "drop"
  )


print(quad_resid)


# Moran’s I test

library(spdep)
coords <- as.matrix(quad_resid[, c("lon", "lat")])

# Define neighbour structure using 4 nearest neighbours per quadrat
nb <- knn2nb(knearneigh(coords, k = 4))

# Convert neighbour object to weights (row-standardised)
lw <- nb2listw(nb, style = "W")

# Moran’s I: tests for positive spatial autocorrelation in residuals
moran.test(quad_resid$resid, lw)  ### p > 0.05 so no spatial terms needed



################################################################################
# M) Model selection and refinement (Beta GLMM with dispersion submodel)
#
# Key modelling choices:
#   - Response: herb_beta (boundary-adjusted proportion)
#   - Random effect: QUADRAT (sampling structure; trees nested in quadrats)
#   - Fixed effects (a priori): QUADRAT_TYPE, PLANT_SPECIES
#   - Candidate covariates: leaf area, DBH, CaCo (scaled)
#   - Dispersion (precision) model: allow dispersion to vary by QUADRAT_TYPE and CaCo
#
# Diagnostics:
#   - DHARMa tests: dispersion, uniformity, outliers
#   - AIC used for comparing candidate fixed-effect structures (same data rows)
################################################################################

library(dplyr)
library(glmmTMB)
library(DHARMa)
library(performance)
library(emmeans)

# Build modelling data frame with scaled predictors
dat <- dat_tiwoho %>%
  mutate(
    QUADRAT = factor(QUADRAT),
    QUADRAT_TYPE = factor(QUADRAT_TYPE),
    PLANT_SPECIES = factor(PLANT_SPECIES),
    
    # Scale continuous predictors to mean 0, SD 1
    # This improves numerical stability and convergence for GLMMs
    leafarea_z = as.numeric(scale(TOTAL_RESTORED_LEAF_AREA_mm)),
    DBH_z      = as.numeric(scale(as.numeric(DBH_cm))),
    CaCo_z     = as.numeric(scale(CaCo_INDEX / 100))  # percent -> proportion -> scale
  )


### Base model: decide PLANT_SPECIES as fixed vs random effect

# Species as random intercept:
m_base_RE <- glmmTMB(
  herb_beta ~ QUADRAT_TYPE +
    (1 | QUADRAT) +
    (1 | PLANT_SPECIES),
  family = beta_family(link = "logit"),
  data = dat
)

# Species as fixed effect:
m_base_FE <- glmmTMB(
  herb_beta ~ QUADRAT_TYPE + PLANT_SPECIES +
    (1 | QUADRAT),
  family = beta_family(link = "logit"),
  data = dat
)

# Compare AIC (lower is better fit-adjusted-for-complexity)
AIC(m_base_RE, m_base_FE)   ### suggests plant as fixed is better

# Diagnose both to understand which misspecification is occurring
res_FE <- simulateResiduals(m_base_FE, n = 2000)  
plot(res_FE)
testDispersion(res_FE)
testUniformity(res_FE)
testOutliers(res_FE)

res_RE <- simulateResiduals(m_base_RE, n = 2000)
plot(res_RE)
testDispersion(res_RE)
testUniformity(res_RE)
testOutliers(res_RE)

# Decision 
#   - Proceed with plant as fixed effect, then improve dispersion via dispformula.


###Improve dispersion: allow dispersion to differ by QUADRAT_TYPE

m_base_FE_disp <- glmmTMB(
  herb_beta ~ QUADRAT_TYPE + PLANT_SPECIES +
    (1 | QUADRAT),
  dispformula = ~ QUADRAT_TYPE,
  family = beta_family(),
  data = dat
)

# Diagnostic checks after adding dispersion structure
res_FE_disp <- simulateResiduals(m_base_FE_disp, n = 2000)
plot(res_FE_disp)
testDispersion(res_FE_disp)
testUniformity(res_FE_disp)
testOutliers(res_FE_disp)


###Add candidate covariates and compare via AIC
# Important:
#   - AIC comparisons must use models fit to the SAME ROWS.
#   - Therefore, filter to complete cases for the candidate covariates first.

dat_sel <- dat %>%
  filter(
    !is.na(leafarea_z),
    !is.na(DBH_z),
    !is.na(CaCo_z)
  )

# Refit baseline dispersion model on the complete-case dataset
m_base_FE_disp2 <- update(m_base_FE_disp, data = dat_sel)

# Add predictors one-at-a-time, then the full model
m_FE_disp_leaf <- update(m_base_FE_disp, . ~ . + leafarea_z)
m_FE_disp_leaf2 <- update(m_FE_disp_leaf, data = dat_sel)

m_FE_disp_dbh <- update(m_base_FE_disp, . ~ . + DBH_z)
m_FE_disp_dbh2  <- update(m_FE_disp_dbh,  data = dat_sel)

m_FE_disp_caco <- update(m_base_FE_disp, . ~ . + CaCo_z)
m_FE_disp_caco2 <- update(m_FE_disp_caco, data = dat_sel)

m_FE_disp_full <- update(m_base_FE_disp,. ~ . + leafarea_z + DBH_z + CaCo_z)
m_FE_disp_full2 <- update(m_FE_disp_full, data = dat_sel)

# AIC comparison across candidate models (full model lowest AIC in the run)
AIC(
  m_base_FE_disp2,
  m_FE_disp_leaf2,
  m_FE_disp_dbh2,
  m_FE_disp_caco2,
  m_FE_disp_full2
)

# Choose best mean-structure model
m_final <- m_FE_disp_full2

# Diagnose the candidate final model
res_final <- simulateResiduals(m_final, n = 2000)
plot(res_final)
testDispersion(res_final)  #still a small dispersion deviation
testUniformity(res_final)
testOutliers(res_final)


# Final refinement:
#   - Add CaCo_z into dispersion model (dispformula)
#   - This can address heteroskedasticity that depends on canopy cover / quadrat context
m_final3 <- update(
  m_final,
  dispformula = ~ QUADRAT_TYPE + CaCo_z
)

# Diagnose the updated model
res3 <- simulateResiduals(m_final3, n = 2000)
plot(res3)
testDispersion(res3)
testUniformity(res3)
testOutliers(res3)

# At this stage: diagnostics indicate an acceptable fit in the workflow.

################################################################################
# N) Final model refit with explicit reference levels (interpretation)
# Why:
#   - Releveling factors changes which category is the baseline comparison.
#   - It does not change fitted values, but makes coefficient interpretation
#     aligned with your narrative (reference forest baseline, common species baseline).
################################################################################

dat_sel <- dat_sel %>%
  mutate(
    PLANT_SPECIES = relevel(PLANT_SPECIES, ref = "Rhizophora apiculata"),
    QUADRAT_TYPE  = relevel(QUADRAT_TYPE,  ref = "Reference Forest")
  )

# Refit the final model using the relevelled data frame
m_final3 <- glmmTMB(
  herb_beta ~ QUADRAT_TYPE + PLANT_SPECIES +
    leafarea_z + DBH_z + CaCo_z +
    (1 | QUADRAT),
  dispformula = ~ QUADRAT_TYPE + CaCo_z,
  family = beta_family(link = "logit"),
  data = dat_sel
)

# Print model summary (coefficients on link scale + random effect variance)
summary(m_final3)

# Create a formatted summary table
library(sjPlot)
tab_model(m_final3, show.est = TRUE, show.se = TRUE, show.stat = TRUE)

# Save model object for reproducibility (reload with readRDS() later)
saveRDS(m_final3, "Tiwoho_tree_betaGLMM_final.rds")

################################################################################
# O) Export results in publishable/reusable formats
################################################################################

#Export fixed-effect coefficients (tidy format, with Wald confidence intervals)
library(broom.mixed)

coef_tab <- tidy(
  m_final3,
  effects = "fixed",
  conf.int = TRUE,
  conf.method = "Wald"
)
coef_tab
write.csv(coef_tab,
          "Tiwoho_betaGLMM_coefficients.csv",
          row.names = FALSE)


###Compute R2 for mixed models (Nakagawa marginal + conditional)
#     - Marginal R2: variance explained by fixed effects
#     - Conditional R2: variance explained by fixed + random effects
library(performance)
r2_vals <- r2_nakagawa(m_final3)
r2_vals


###Estimated marginal means (predicted means adjusted for covariates)
#     - Useful for clear reporting of differences among QUADRAT_TYPE and species
library(emmeans)

emm_qt <- emmeans(m_final3, ~ QUADRAT_TYPE, type = "response")
emm_sp <- emmeans(m_final3, ~ PLANT_SPECIES, type = "response")
emm_qt
emm_sp

write.csv(as.data.frame(emm_qt),
          "Tiwoho_betaGLMM_emmeans_quadrat_type.csv",
          row.names = FALSE)

write.csv(as.data.frame(emm_sp),
          "Tiwoho_betaGLMM_emmeans_species.csv",
          row.names = FALSE)


###Export DHARMa diagnostic plots (for supplement / documentation)
library(DHARMa)

res_final <- simulateResiduals(m_final3, n = 2000)
plot(res_final)

pdf("Tiwoho_betaGLMM_DHARMa_diagnostics.pdf", width = 8, height = 8)
plot(res_final)
dev.off()


###Figures
# Figure types included:
#   1) Model-predicted means with confidence intervals (emmeans)
#   2) Raw data distribution plots (boxplot + jitter)
#   3) Continuous predictor response curves from prediction grids (predict.glmmTMB)

library(ggplot2)


###Quadrat/forest type effect figures

##Model-predicted differences among QUADRAT_TYPE (emmeans, response scale)
emm_qt_df <- as.data.frame(emm_qt)

p1a <- ggplot(emm_qt_df,
              aes(x = QUADRAT_TYPE,
                  y = response)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.15
  ) +
  ylab("Predicted herbivory proportion") +
  xlab("Quadrat type") +
  labs(title = "Tiwoho") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )
p1a

ggsave("Fig1a_predicted_herbivory_quadrat_type.png",
       p1a, width = 6, height = 4, dpi = 600)

##Raw distribution by QUADRAT_TYPE (boxplot + jitter) with consistent colours
p1b <- ggplot(dat_sel,
              aes(x = QUADRAT_TYPE,
                  y = HERBIVORY_PROPORTION,
                  fill = QUADRAT_TYPE,
                  colour = QUADRAT_TYPE)) +
  geom_boxplot(outlier.shape = NA,
               alpha = 0.6,
               width = 0.65) +
  geom_jitter(width = 0.15,
              alpha = 0.4,
              size = 1.5) +
  ylab("Observed herbivory proportion") +
  xlab("Quadrat type") +
  labs(title = "Tiwoho") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
p1b

ggsave("Fig1b_herbivory_quadrat_type.png",
       p1a, width = 6, height = 4, dpi = 600)


###Plant species effect figures

##Model-predicted differences among PLANT_SPECIES (emmeans)
emm_sp_df <- as.data.frame(emm_sp)

p2a <- ggplot(emm_sp_df,
              aes(x = PLANT_SPECIES,
                  y = response)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.15
  ) +
  ylab("Predicted herbivory proportion") +
  xlab("Plant species") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               vjust = 1),
    axis.title.x = element_text(margin = margin(t = 12))
  )

p2a

##Raw distribution by PLANT_SPECIES (boxplot + jitter)
p2b <- ggplot(dat_sel,
              aes(x = PLANT_SPECIES,
                  y = HERBIVORY_PROPORTION)) +
  geom_boxplot(outlier.shape = NA,
               fill = "grey90",
               colour = "black",
               width = 0.6) +
  geom_jitter(width = 0.15,
              alpha = 0.5,
              size = 1.8) +
  ylab("Observed herbivory proportion") +
  xlab("Plant species") +
  labs(title = "Tiwoho") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               vjust = 1),
    axis.title.x = element_text(margin = margin(t = 12)),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p2b

ggsave("Fig2b_herbivory_plant_species.png",
       p2a, width = 6, height = 4, dpi = 600)


###Continuous predictor effects using prediction grids (CaCo and DBH)
# Why a prediction grid?
#   - One predictor is caried across across its observed range (e.g. CaCo_z),
#     while holding other predictors at typical values (0 for scaled covariates).
#   - We then predict fitted means and approximate 95% CIs using se.fit.
# Important detail:
#   - CaCo_z and DBH_z are scaled. We back-transform to original units for plotting.


###Plot the impact of vegetation density (CaCo)

##Recover scaling parameters explicitly:
#CaCo_z was created as scale(CaCo_INDEX/100), so we reconstruct mean and sd
caco_raw_prop <- dat_sel$CaCo_INDEX / 100
caco_mean <- mean(caco_raw_prop, na.rm = TRUE)
caco_sd   <- sd(caco_raw_prop,   na.rm = TRUE)


###Build prediction grid:
#  - CaCo_z varies across observed range
#  - leafarea_z and DBH_z set to 0 (mean of scaled variables)
#  - QUADRAT_TYPE varies to produce separate lines by forest type
#  - PLANT_SPECIES fixed to the first factor level (reference for this plot)
#  - QUADRAT = NA so predictions are at population level (no specific quadrat)

newdat <- expand.grid(
  CaCo_z = seq(min(dat_sel$CaCo_z, na.rm = TRUE),
               max(dat_sel$CaCo_z, na.rm = TRUE),
               length.out = 100),
  leafarea_z   = 0,
  DBH_z        = 0,
  QUADRAT_TYPE = levels(dat_sel$QUADRAT_TYPE),
  PLANT_SPECIES = levels(dat_sel$PLANT_SPECIES)[1],
  QUADRAT      = NA
)


##Predict fitted values and SE on the response scale
pred <- predict(m_final3, newdata = newdat, type = "response", se.fit = TRUE)


##Attach predictions + compute approximate 95% CI
#Also back-transform CaCo_z to CaCo_INDEX (%) for interpretability
newdat <- newdat %>%
  mutate(
    fit = pred$fit,
    se  = pred$se.fit,
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se,
    CaCo_INDEX_pct = 100 * (CaCo_z * caco_sd + caco_mean)
  )


##Plot fitted lines + CI ribbons, colour by QUADRAT_TYPE
p3 <- ggplot(newdat,
             aes(x = CaCo_INDEX_pct,
                 y = fit,
                 colour = QUADRAT_TYPE,
                 fill = QUADRAT_TYPE)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              alpha = 0.2,
              colour = NA) +
  xlab("Caco - Canopy cover index (%)") +
  ylab("Predicted herbivory proportion") +
  theme_classic(base_size = 14)

p3

##Save figure in high-resolution PNG and vector PDF
ggsave("Fig3_predicted_herbivory_CaCo.png",
       p3, width = 6.5, height = 4.5, dpi = 600)
ggsave("Fig3_predicted_herbivory_CaCo.pdf",
       p3, width = 6.5, height = 4.5)



##canopy cover by forest type 

library(dplyr)

caco_by_type <- dat_sel %>%
  group_by(QUADRAT_TYPE) %>%
  summarise(
    mean_caco = mean(CaCo_INDEX, na.rm = TRUE),
    sd_caco = sd(CaCo_INDEX, na.rm = TRUE),
    n = n()
  )

caco_by_type
ggplot(dat_sel, aes(x = QUADRAT_TYPE, y = CaCo_INDEX)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  ylab("Canopy cover index (%)") +
  xlab("Forest type") +
  theme_classic()


###Plot the impact of tree size (DBH)

##Recover scaling parameters for DBH_z = scale(DBH_cm)
dbh_mean <- mean(dat_sel$DBH_cm, na.rm = TRUE)
dbh_sd   <- sd(dat_sel$DBH_cm,   na.rm = TRUE)

#Prediction grid:
#  - DBH_z varies; CaCo_z fixed at 0 (mean)
#  - leafarea_z fixed at 0 (mean)
#  - QUADRAT_TYPE varies (separate lines)
#  - species fixed to reference level for this plot
newdat_dbh <- expand.grid(
  DBH_z = seq(min(dat_sel$DBH_z, na.rm = TRUE),
              max(dat_sel$DBH_z, na.rm = TRUE),
              length.out = 100),
  leafarea_z   = 0,
  CaCo_z       = 0,
  QUADRAT_TYPE = levels(dat_sel$QUADRAT_TYPE),
  PLANT_SPECIES = levels(dat_sel$PLANT_SPECIES)[1],
  QUADRAT      = NA
)


##Predictions + SE on response scale
pred_dbh <- predict(m_final3, newdata = newdat_dbh, type = "response", se.fit = TRUE)

##Attach predictions + CI; back-transform DBH_z to DBH in cm
newdat_dbh <- newdat_dbh %>%
  mutate(
    fit = pred_dbh$fit,
    se  = pred_dbh$se.fit,
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se,
    DBH_cm = (DBH_z * dbh_sd + dbh_mean)
  )


##Plot DBH effect with quadrat-type colours
p4 <- ggplot(newdat_dbh,
             aes(x = DBH_cm,
                 y = fit,
                 colour = QUADRAT_TYPE,
                 fill = QUADRAT_TYPE)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              alpha = 0.2,
              colour = NA) +
  xlab("DBH (cm)") +
  ylab("Predicted herbivory proportion") +
  theme_classic(base_size = 14)

p4

##Save figure outputs
ggsave("Fig4_predicted_herbivory_DBH.png",
       p4, width = 6.5, height = 4.5, dpi = 600)
ggsave("Fig4_predicted_herbivory_DBH.pdf",
       p4, width = 6.5, height = 4.5)


##DBH cover by forest type 

library(dplyr)

dbh_by_type <- dat_sel %>%
  group_by(QUADRAT_TYPE) %>%
  summarise(
    mean_dbh = mean(DBH_cm, na.rm = TRUE),
    sd_dbh = sd(DBH_cm, na.rm = TRUE),
    n = n()
  )

dbh_by_type
ggplot(dat_sel, aes(x = QUADRAT_TYPE, y = DBH_cm)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  ylab("DBH (cm)") +
  xlab("Forest type") +
  theme_classic()


###sjPlot figures (alternative output style)

library(sjPlot)
library(ggplot2)
library(dplyr)
library(emmeans)


# Define a re-usable publication theme for sjPlot-based ggplots
pub_theme <- function() {
  theme_classic(base_size = 14) +
    theme(
      legend.key.size = unit(0.9, "cm"),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text  = element_text(size = 11),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x  = element_text(size = 12, colour = "black"),
      axis.text.y  = element_text(size = 12, colour = "black"),
      strip.text   = element_text(size = 12, face = "bold", colour = "black")
    )
}


### CaCo effect plot (x-axis on 0–100%)

# Compute mean/sd for back-transforming z-scale to percent scale for labelling
caco_prop <- dat_sel$CaCo_INDEX / 100
caco_mean <- mean(caco_prop, na.rm = TRUE)
caco_sd   <- sd(caco_prop,   na.rm = TRUE)

# Create z-scale tick positions and map them to percent labels
z_breaks <- pretty(range(dat_sel$CaCo_z, na.rm = TRUE), n = 5)
pct_labels <- round(100 * (z_breaks * caco_sd + caco_mean), 1)

# sjPlot prediction plot (note: dispersion submodel may be ignored in plot_model)
p_caco <- plot_model(
  m_final3,
  type  = "pred",
  terms = c("CaCo_z [all]", "QUADRAT_TYPE")
) +
  labs(x = "Canopy cover index (%)",
       y = "Predicted herbivory proportion",
       colour = "Quadrat type",
       fill = "Quadrat type",
       title = "Tiwoho") +
  scale_x_continuous(breaks = z_breaks, labels = pct_labels) +
  pub_theme() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_caco
ggsave("Fig_effect_CaCo_quadratType.png", p_caco, width = 7, height = 5, dpi = 600)
ggsave("Fig_effect_CaCo_quadratType.pdf", p_caco, width = 7, height = 5)


### DBH

# Compute mean/sd for back-transforming DBH_z ticks to DBH_cm
dbh_mean <- mean(dat_sel$DBH_cm, na.rm = TRUE)
dbh_sd   <- sd(dat_sel$DBH_cm,   na.rm = TRUE)

z_breaks_dbh <- pretty(range(dat_sel$DBH_z, na.rm = TRUE), n = 5)
cm_labels <- round(z_breaks_dbh * dbh_sd + dbh_mean, 1)

# sjPlot prediction plot for DBH
p_dbh <- plot_model(
  m_final3,
  type  = "pred",
  terms = c("DBH_z [all]", "QUADRAT_TYPE")
) +
  labs(x = "DBH (cm)",
       y = "Predicted herbivory proportion",
       colour = "Quadrat type",
       fill = "Quadrat type",
       title = "Tiwoho") +
  scale_x_continuous(breaks = z_breaks_dbh, labels = cm_labels) +
  pub_theme() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_dbh
ggsave("Fig_effect_DBH_quadratType.png", p_dbh, width = 7, height = 5, dpi = 600)
ggsave("Fig_effect_DBH_quadratType.pdf", p_dbh, width = 7, height = 5)


#### quadrat

# sjPlot predicted means for quadrat type (categorical predictor)
p_qtype <- plot_model(
  m_final3,
  type = "pred",
  terms = "QUADRAT_TYPE"
) +
  labs(x = "Quadrat type",
       y = "Predicted herbivory proportion") +
  pub_theme() +
  theme(legend.position = "none") +
  ggtitle("")

p_qtype
ggsave("Fig_effect_QuadratType.png", p_qtype, width = 6, height = 4, dpi = 600)
ggsave("Fig_effect_QuadratType.pdf", p_qtype, width = 6, height = 4)


### plant species

# For many species levels, an emmeans dot + CI plot is usually clearer
emm_sp <- emmeans(m_final3, ~ PLANT_SPECIES, type = "response")
emm_sp_df <- as.data.frame(emm_sp)

p_species <- ggplot(emm_sp_df,
                    aes(x = reorder(PLANT_SPECIES, response),
                        y = response)) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  coord_flip() +
  labs(x = "Plant species",
       y = "Predicted herbivory proportion",
       title = "Tiwoho") +
  pub_theme() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p_species
ggsave("Fig_effect_PlantSpecies.png", p_species, width = 7.5, height = 6.5, dpi = 600)
ggsave("Fig_effect_PlantSpecies.pdf", p_species, width = 7.5, height = 6.5)

  