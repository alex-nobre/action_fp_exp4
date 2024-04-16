#================================================================================================================#
# This fits linear mixed models using a transformation and a model structure obtained separately
# ("Analysis.notebook.html"). 
#================================================================================================================#

# Load packages

# Data processing and plotting
library(magrittr)
library(tidyverse)
library(lattice)
library(gridExtra)
library(data.table)

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

# emm options
emm_options(lmer.df = "satterthwaite", lmerTest.limit = 12000)


#================================== 0. Read data ================================
# Create dataset
source('./Analysis/Prepare_data.R')

# Set contrasts
contrasts(data$foreperiod) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(data$logFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(data$condition) <- c(-1/2, 1/2)
contrasts(data$prevOri) <- c(-1/2, 1/2)
contrasts(data$oneBackFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)

contrasts(data2$foreperiod) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(data2$logFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(data2$condition) <- c(-1/2, 1/2)
contrasts(data2$prevOri) <- c(-1/2, 1/2)
contrasts(data2$oneBackFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)


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

trimlogfplmm <- mixed(formula =  logRT ~ 1 + foreperiod*condition*oneBackFP + 
                         (1 + foreperiod + condition + foreperiod:condition + oneBackFP | participant),
                       data = data,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")



# Systematic comparisons betweeen lmm's via BIC


# Random-intercept only model
trimlogfplmm1v2 <- mixed(logRT ~ 1 + foreperiod + condition + oneBackFP + foreperiod:condition + 
                           foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                         (1 | participant),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")


# Compare BICs and AICs
BIC(trimlogfplmm1, trimlogfplmm1v2) %>%
  kable()

AIC(trimlogfplmm1, trimlogfplmm1v2) %>%
  kable()

cor(fitted(trimlogfplmm1), data2$logRT)^2
cor(fitted(trimlogfplmm1v2), data2$logRT)^2

# Better BICs/AICs for intercept-only model, but by a small margin; those are also less plausible according to the data

#==============================================================================================#
#==================================== 3. Model assessment ======================================
#==============================================================================================#

#============= 3.1. Using logRT ===============
anova(trimlogfplmm)

#Visualize random effects
dotplot(ranef(trimlogfplmm, condVar = TRUE))

# Visualize interactions
emmip(trimlogfplmm,
      oneBackFP ~ condition, style = "factor") # Using averaged FP

# Single slopes tests
fp_by_condition <- slopes(trimlogfplmm, by = "condition", variables = "foreperiod",
                          p_adjust = "holm")

test(emtrends(trimlogfplmm, ~ condition, var="foreperiod")) # equivalent to slopes

fp_by_oneback <- slopes(trimlogfplmm, by = "oneBackFP", variables = "scaledNumForeperiod",
                        p_adjust = "holm")

test(emtrends(trimlogfplmm, ~ oneBackFP, var="scaledNumForeperiod")) # equivalent to slopes

threeway_int <- slopes(trimlogfplmm, by = c("oneBackFP", "condition"), variables = "scaledNumForeperiod",
                       p_adjust = "holm")


# Pairwise comparisons
fp_by_condition_comp <- emtrends(trimlogfplmm, "condition", var = "scaledNumForeperiod")
fp_by_condition_comp
update(pairs(fp_by_condition_comp), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimlogfplmm, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

threeway_int_comp = emtrends(trimlogfplmm, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
threeway_int_comp
update(pairs(threeway_int_comp), by = NULL, adjust = "holm")
pairs(threeway_int_comp, simple = "condition")

# Marginal means for FP by FP n-1 and condition
oneback_by_cond = emmeans(trimlogfplmm, c("condition", "oneBackFP"))
oneback_by_cond = emmeans(trimlogfplmm, ~ oneBackFP * condition)
pairs(oneback_by_cond, simple = "condition")


#============ 3.5. Hierarchical entry ===============
h_trimlogfplmm1 <- mixed(logRT ~ 1 + numForeperiod + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method = 'KR',
                         REML=TRUE,
                         return = "merMod")

h_trimlogfplmm2 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm3 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm4 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm5 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + 
                           condition:numForeperiod + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm6 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + 
                           condition:numForeperiod + condition:numOneBackFP + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm7 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + 
                           condition:numForeperiod + condition:numOneBackFP + condition:numOneBackFP:numForeperiod + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)


h_trimlogfplmm8 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + 
                           condition:numForeperiod + condition:numOneBackFP + condition:numOneBackFP:numForeperiod + (1 + condition | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

anova(h_trimlogfplmm1, h_trimlogfplmm2, h_trimlogfplmm3, h_trimlogfplmm4, 
      h_trimlogfplmm5, h_trimlogfplmm6, h_trimlogfplmm7, h_trimlogfplmm8)


#==================== 3.3. Run model separately for action and external conditions ===================

#======== 3.3.1. External ============#
# FP as numeric
trimlogfplmmext1 <- mixed(logRT ~ 1 + numForeperiod + 
                            numOneBackFP + numForeperiod:numOneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext1)


# FP as categorical
trimlogfplmmext2 <- mixed(logRT ~ 1 + foreperiod + 
                            numOneBackFP + foreperiod:numOneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext2)

# logFP as numerical
trimlogfplmmext3 <- mixed(logRT ~ 1 + numLogFP + 
                            numOneBackFP + numLogFP:numOneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext3)

# logFP as categorical
trimlogfplmmext4 <- mixed(logRT ~ 1 + logFP + 
                            numOneBackFP + logFP:numOneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext4)

# FP as numerical Including only FP
trimlogfplmmext5 <- mixed(logRT ~ 1 + numForeperiod + 
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext5)
summary(trimlogfplmmext1)

BIC(trimlogfplmmext1, trimlogfplmmext5)

# Including FP n-1 improves BIC

# FP as categorical Including only FP
trimlogfplmmext6 <- mixed(logRT ~ 1 + foreperiod + 
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext6)
summary(trimlogfplmmext2)

var(lme4::fixef(trimlogfplmmext6))
var(lme4::fixef(trimlogfplmmext2))

r.squaredGLMM(trimlogfplmmext6)
r.squaredGLMM(trimlogfplmmext2)

BIC(trimlogfplmmext2, trimlogfplmmext6)

# logFP as numerical Including only FP
trimlogfplmmext7 <- mixed(logRT ~ 1 + numLogFP + 
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext7)
summary(trimlogfplmmext3)

BIC(trimlogfplmmext3, trimlogfplmmext7)

# Including FP n-1 improves BIC

# logFP as categorical Including only FP
trimlogfplmmext8 <- mixed(logRT ~ 1 + logFP + 
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext8)
summary(trimlogfplmmext4)

BIC(trimlogfplmmext4, trimlogfplmmext8)

# FP as numeric, FP n-1 as categorical
trimlogfplmmext9 <- mixed(logRT ~ 1 + numForeperiod + 
                            oneBackFP + numForeperiod:oneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext9)
anova(trimlogfplmmext9)

#======== 3.3.2. Action ============#
trimlogfplmm1act <- mixed(logRT ~ 1 + numForeperiod + 
                            numOneBackFP + numForeperiod:numOneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='action',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmm1act)


trimlogfplmm1 <- mixed(logRT ~ 1 + foreperiod +
                         (1 | ID),
                       data=data2[data2$condition=='external',],
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)
summary(trimlogfplmm1)

anova(trimlogfplmm1)

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

fpaccglmer <- buildmer(error_result ~ foreperiod * condition * oneBackFP +
                         (1+foreperiod*condition*oneBackFP|participant), 
                       data=dataAcc,
                       family = binomial(link = "logit"),
                       buildmerControl = buildmerControl(calc.anova = TRUE,
                                                         ddf = "Satterthwaite",
                                                         include = 'foreperiod*condition*oneBackFP'))


formula(fpaccglmer)
isSingular(fpaccglmer)

fpaccglmer <- mixed(formula = error_result ~ foreperiod * condition * oneBackFP +
                      (1 + condition | participant),
                    data = dataAcc,
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5)),
                    progress = TRUE,
                    expand_re = FALSE,
                    return = "merMod",
                    method = "LRT")


isSingular(fpaccglmer)

summary(fpaccglmer)
anova(fpaccglmer)

