#================================================================================================================#
# This fits linear mixed models using a transformation and a model structure obtained separately
# ("Analysis.notebook.html"). 
#================================================================================================================#

# Load packages

# Data processing
library(magrittr)
library(tidyverse)
library(data.table)

# Plotting
library(lattice)
library(gridExtra)
library(extrafont)

# Simple models
library(car)
library(janitor)

# Mixed-effects modeling
library(afex)
library(emmeans)
library(lme4)
library(MuMIn)
library(buildmer)
library(broom.mixed)
library(marginaleffects)

# Bayesian models
library(brms)
library(bayestestR)
library(BayesFactor)

# Assess models and results
library(effects)
library(ggeffects)
library(performance)
library(knitr)
library(kableExtra)
library(sjPlot)
library(prediction)


# Save defaults
graphical_defaults <- par()
options_defaults <- options()

# Load fonts from extrafonts
loadfonts()

# Prepare data 
source('./Analysis/Prepare_data.R')

# Prepare theme for plots
source("./Analysis/plot_theme.R")
theme_set(mytheme)

# emm options
emm_options(lmer.df = "satterthwaite", lmerTest.limit = 12000)


#==========================================================================================#
#================================= 1. Explore individual data ==============================
#==========================================================================================#


# Plot RT by FP by participant and model using complete pooling 
FPfitAll=lm(meanRT ~ foreperiod,
            data=summaryData2)

fit.params=tidy(FPfitAll)

summary(FPfitAll)


ggplot(data=summaryData2,
       aes(x=foreperiod,
           y=meanRT)) +
  stat_summary(fun="mean", geom="point", size=1.5)+
  geom_abline(intercept=fit.params$estimate[1],
              slope=fit.params$estimate[2],
              color="blue")+
  facet_wrap(~ participant, ncol=6)


# Plot RT by FP by participant and model using individual data (no pooling)
dataGroupedByRT <- summaryData2 %>% 
  group_by(participant,foreperiod) %>% 
  summarise(meanRT=mean(meanRT)) %>%
  ungroup()

data.no_pooling <- dataGroupedByRT %>%
  group_by(participant) %>%
  nest(data = c(foreperiod, meanRT)) %>% # Nest data frames in columns for each row combination of values
  mutate(fit = map(data, ~ lm(meanRT ~ foreperiod, data = .)), # fit linear model for each nested dataset
         params = map(fit, tidy)) %>%
  ungroup() %>%
  unnest(c(params)) %>% # transform nested parameters to columns
  select(participant, term, estimate) %>%
  # creates all combinations of levels for each factor inputted as argument, regardless of whether they appear in the data
  complete(participant, term, fill = list(estimate = 0)) %>%
  pivot_wider(names_from = term,
              values_from = estimate) %>% 
  clean_names() %>% # remove non-ascii symbols from names
  rename(foreperiod=foreperiod2_8)

# Plot slopes for each participant along with mean RT by FP
ggplot(data = dataGroupedByRT,
       aes(x = foreperiod, y = meanRT)) + 
  geom_abline(data = data.no_pooling,
              aes(intercept = intercept,
                  slope = foreperiod),
              color = "blue") +
  geom_point() +
  facet_wrap(~participant, ncol=6) + 
  theme(strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))

fp_no_pooling <- data.no_pooling$foreperiod

# Compare results to see how much they differ
data_grouped_by_fp <- data %>%
  group_by(participant, foreperiod) %>%
  summarise(meanRT=mean(RT)) %>%
  ungroup()

fit_partial_pooling <- lmer(formula = meanRT ~ foreperiod + 
                              (1 + foreperiod|participant),
                            data = data_grouped_by_fp)

data_partial_pooling <- fit_partial_pooling %>%
  augment() %>%
  select(participant, numForeperiod, RT, .fitted) %>%
  rename(fitted=.fitted)

#================ 2.2. FP n-1 ==================
dataGroupedByRT <- summaryData %>% 
  group_by(ID,oneBackFP) %>% 
  summarise(meanRT=mean(meanRT)) %>%
  ungroup() %>%
  mutate(numOneBackFP = as.numeric(as.character(oneBackFP)))

data.no_pooling <- dataGroupedByRT %>%
  select(-oneBackFP) %>%
  group_by(ID) %>%
  nest(data = c(numOneBackFP, meanRT)) %>%
  mutate(fit = map(data, ~ lm(meanRT ~ numOneBackFP, data = .)),
         params = map(fit, tidy)) %>%
  ungroup() %>%
  unnest(c(params)) %>%
  select(ID, term, estimate) %>%
  complete(ID, term, fill = list(estimate = 0)) %>%
  pivot_wider(names_from = term,
              values_from = estimate) %>% 
  clean_names()


data.no_pooling <- data.no_pooling %>%
  rename(ID=id,
         numOneBackFP=num_one_back_fp)


ggplot(data = dataGroupedByRT,
       aes(x = numOneBackFP, y = meanRT)) + 
  geom_abline(data = data.no_pooling,
              aes(intercept = intercept,
                  slope = numOneBackFP),
              color = "blue") +
  geom_point() +
  facet_wrap(~ID, ncol=6) + 
  scale_x_continuous(breaks = 0:4 * 2) +
  theme(strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))


#==========================================================================================#
#====================================== 2. Prepare model ====================================
#==========================================================================================#

fplmm <- mixed(formula =  logRT ~ 1 + foreperiod + condition + foreperiod:condition + oneBackFP + 
                 foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                 (1 + condition + foreperiod + foreperiod:condition | participant),
               data = data2,
               control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
               progress = TRUE,
               expand_re = TRUE,
               method =  'S',
               REML=TRUE,
               return = "merMod")




#==============================================================================================#
#==================================== 3. Model assessment ======================================
#==============================================================================================#

#============= 3.1. Using logRT ===============
anova(fplmm)

saveRDS(fplmm, file = "./Analysis/Sensitivity_analysis/fplmm.rds")

#Visualize random effects
dotplot(ranef(fplmm, condVar = TRUE))

#========== 3.2. Two-way interactions ==========
# 3.2.1. Compare difference between conditions within each level of FPn
fp_cond_emm <- emmeans(fplmm, ~ condition|foreperiod)
contrast(fp_cond_emm, interaction = c("pairwise"), adjust = "holm")

# 3.2.2. Exame effect of FP within each condition
cond_emm <- emmeans(fplmm, ~ foreperiod|condition)
contrast(cond_emm, interaction = c("pairwise"), adjust = "holm")

# Compare FP effect between conditions
cond_emm <- emmeans(fplmm, ~ foreperiod*condition)
contrast(cond_emm, interaction = c("pairwise"), adjust = "holm")

# 3.2.3. Compare difference between consecutive levels of FPn-1 for each level of FPn
fp_onebackfp_emm <- emmeans(fplmm, ~ oneBackFP|foreperiod)
contrast(fp_onebackfp_emm, interaction = c("consec"), adjust = "holm")

fp_onebackfp_emm <- emmeans(fplmm, ~ foreperiod|oneBackFP)
contrast(fp_onebackfp_emm, interaction = c("poly"), adjust = "holm", max.degree = 1)

#========== 3.3 Three-way interaction =========
# Visualize interactions
emmip(fplmm,
      oneBackFP ~ foreperiod|condition, style = "factor") # Using averaged FP

fpemm <- emmeans(fplmm, ~ foreperiod*oneBackFP|condition)
contrast(fpemm, interaction = c("pairwise", "pairwise"), adjust = "holm")

fpemm <- emmeans(fplmm, ~ foreperiod*oneBackFP*condition)
contrast(fpemm, interaction = c("pairwise", "pairwise", "pairwise"), adjust = "holm")



emmip(fplmm,
      condition ~ foreperiod|oneBackFP, style = "factor") # Using averaged FP

# andré
fpemm <- emmeans(fplmm, ~ foreperiod|oneBackFP*condition)
contrast(fpemm, interaction = c("pairwise"), adjust = "holm")


# andré
fpemm <- emmeans(fplmm, ~ foreperiod*condition|oneBackFP)
contrast(fpemm, interaction = c("pairwise", "pairwise"), adjust = "holm")



#==============================================================================================#
#================================== 4. Choose distribution ====================================
#==============================================================================================#

#======================== 4.1. Visualize distributions for each variable =======================
# Pooled
RTHistograms <- ggplot(data=data2,
                       aes(x=RT))+
  geom_histogram()
RTHistograms

invRTHistograms <- ggplot(data=data2,
                          aes(x=invRT)) +
  geom_histogram()
invRTHistograms

logRTHistograms <- ggplot(data=data2,
                          aes(x=logRT)) +
  geom_histogram()
logRTHistograms


# By participant
indRTHistograms <- ggplot(data=data2,
                          aes(x=RT))+
  geom_histogram()+
  facet_wrap(~ID)
indRTHistograms

indinvRTHistograms <- ggplot(data=data2,
                             aes(x=invRT)) +
  geom_histogram() +
  facet_wrap(~ID)
indinvRTHistograms

indlogRTHistograms <- ggplot(data=data2,
                             aes(x=logRT)) +
  geom_histogram() +
  facet_wrap(~ID)
indlogRTHistograms

# Try models
trimfpgauss <- mixed(formula = RT ~ 1 + foreperiod + condition + oneBackFP + foreperiod:condition + 
                       foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                       (1 + foreperiod + condition | participant),
                        data = data2,
                        #family=gaussian(link = "identity"),
                        control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'KR',
                        return = "merMod")

summary(trimfpgauss)

trimfpinvgauss <- mixed(formula = RT ~ 1 + foreperiod + condition + oneBackFP + foreperiod:condition + 
                             foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                             (1 + foreperiod + condition | participant),
                           data = data2,
                           family=inverse.gaussian(link = "identity"),
                           control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                           progress = TRUE,
                           expand_re = TRUE,
                           method =  'KR',
                           return = "merMod")

summary(trimfpinvgauss)

trimfpgamma <- mixed(formula = RT ~ 1 + foreperiod + condition + oneBackFP + foreperiod:condition + 
                          foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                          (1 + foreperiod + condition | participant),
                        data = data2,
                        family=Gamma(link = "identity"),
                        control = glmerControl(optimizer = c("nloptwrap"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'KR',
                        return = "merMod")

summary(trimfpgamma)

# Compare visualizations
ggpredict(model = trimfpgauss,
          terms = "foreperiod",
          type = 'fe') %>%
  plot()

ggpredict(model = trimfpinvgauss,
          terms = "foreperiod",
          type = 'fe') %>%
  plot()

ggpredict(model = trimfpgamma,
          terms = 'foreperiod',
          type = 'fe') %>%
  plot()

ggpredict(model = trimfpgauss,
          terms = "condition",
          type = 'fe') %>%
  plot()

ggpredict(model = trimfpinvgauss,
          terms = "condition",
          type = 'fe') %>%
  plot()

ggpredict(model = trimfpgamma,
          terms = 'condition',
          type = 'fe') %>%
  plot()


# Compare performance across models
trimfpgauss %>%
  check_model()

trimfpinvgauss %>%
  check_model()

trimfpgamma %>%
  check_model()

# Compare with log model
trimlogfplmm %>%
  check_model

# The log model seems to provide a better fit than the alternatives

#===========================================================================================#
#=================================== 5. Accuracy ============================================
#===========================================================================================#

fpaccglm <- mixed(formula = error_result ~ 1 + foreperiod:condition:oneBackFP + foreperiod + 
                    condition + oneBackFP + foreperiod:condition + foreperiod:oneBackFP + 
                    condition:oneBackFP + (1 | participant),
                    data = dataAcc,
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5)),
                    progress = TRUE,
                    expand_re = FALSE,
                    method = "LRT")


isSingular(fpaccglm)

anova(fpaccglm)

# Compare difference between levels of FPn-1 within conditions
cond_onebackfp_emm_acc <- emmeans(fpaccglm, ~ condition|oneBackFP)
contrast(cond_onebackfp_emm_acc, interaction = c("pairwise"), adjust = "holm")

acc3way <- emmeans(fpaccglm, ~ condition|oneBackFP*foreperiod)
contrast(acc3way, interaction = c("pairwise"), adjust = "holm")


