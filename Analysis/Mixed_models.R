#================================================================================================================#
# Changes:

# To choose the dependent variable, we fit models using only random intercepts instead of the full
# random-effects structure, as suggested by Salet et al. (2022) in their model structure Rmd
# Removed old code in settings contrasts section
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

# Functions
hist_resid <- function(M,ptitle='Residuals') {
  d <- data.frame(resid=residuals(M)) 
  d  %>% ggplot(aes(x=resid)) + 
    geom_histogram(aes(y=..density..), bins=75, color='black', fill='grey') + 
    geom_density(color='darkred') + 
    ggtitle(ptitle) -> pl
  return(pl)
}

fitstats = function(M,mname='M') {
  QQ<-qqnorm(residuals(M), plot.it=FALSE)
  R2qq <- cor(QQ$x,QQ$y)^2
  dfqq = data.frame(stat='R2qq', V1=R2qq)
  r2tab <- r.squaredGLMM(M)  %>% 
    t  %>% as.data.frame  %>% rownames_to_column(var='stat')  %>% 
    rbind(.,dfqq)
  r2tab$stat = c("$R^2_m$","$R^2_c$",'$R^2_{qq}$' )
  colnames(r2tab) <- c('stat',mname)
  return(r2tab)
}

# Plot RT by FP by participant and model using complete pooling 
FPfitAll=lm(meanRT ~ foreperiod,
            data=summaryData)

fit.params=tidy(FPfitAll)

summary(FPfitAll)


ggplot(data=summaryData,
       aes(x=foreperiod,
           y=meanRT)) +
  stat_summary(fun="mean", geom="point", size=1.5)+
  geom_abline(intercept=fit.params$estimate[1],
              slope=fit.params$estimate[2],
              color="blue")+
  facet_wrap(~ participant, ncol=6)


# Plot RT by FP by participant and model using individual data (no pooling)
dataGroupedByRT <- summaryData %>% 
  group_by(participant,foreperiod) %>% 
  summarise(meanRT=mean(meanRT)) %>%
  ungroup() %>%
  mutate(numForeperiod=as.numeric(as.character(foreperiod)))

data.no_pooling <- dataGroupedByRT %>%
  select(-foreperiod) %>%
  group_by(participant) %>%
  nest(data = c(numForeperiod, meanRT)) %>%
  mutate(fit = map(data, ~ lm(meanRT ~ numForeperiod, data = .)),
         params = map(fit, tidy)) %>%
  ungroup() %>%
  unnest(c(params)) %>%
  select(participant, term, estimate) %>%
  complete(participant, term, fill = list(estimate = 0)) %>%
  pivot_wider(names_from = term,
              values_from = estimate) %>% 
  clean_names()


data.no_pooling <- data.no_pooling %>%
  rename(numForeperiod=num_foreperiod)


ggplot(data = dataGroupedByRT,
       aes(x = numForeperiod, y = meanRT)) + 
  geom_abline(data = data.no_pooling,
              aes(intercept = intercept,
                  slope = numForeperiod),
              color = "blue") +
  geom_point() +
  facet_wrap(~participant, ncol=6) + 
  scale_x_continuous(breaks = 0:4 * 2) +
  theme(strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))

fp_no_pooling <- data.no_pooling$numForeperiod

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

#=========================== 2.1. Choose dependent variable =============================
# We use the strategy of keeping it maximal to find a model that converges and progressively
# remove terms, one of the strategies recommended to avoid overfitting:
# https://rdrr.io/cran/lme4/man/isSingular.html

#  Trying to fit the model using foreperiod and FPn-1 as factors results
# in R hanging during execution; for this reason, we use them as
# numerical variables

# To choose the data transformation that leads to the optimal random effects structure, we fit models including only 
# random intercepts and compare R2 and residuals

# Fit models with RT and inverse RT without trimming
options(scipen = 999)

fplmm1 <- mixed(formula = RT ~ foreperiod*condition*oneBackFP + 
                  (1|participant),
                data = data,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

summary(fplmm1)


# Now we run the same model with inverse RT and logRT as outcomes
invfplmm1 <- mixed(formula = invRT ~ foreperiod*condition*oneBackFP + 
                     (1|participant),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

logfplmm1 <- mixed(formula = logRT ~ foreperiod*condition*oneBackFP + 
                     (1|participant),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

# Let's check that the current structure does not provide a singular fit:
isSingular(fplmm1)
isSingular(invfplmm1)
isSingular(logfplmm1)

# None of them return singular fits!

# Amount of variance accounted for by the model
var <- data.frame(dataset = 'no trim',
                  'RT' = cor(fitted(fplmm1), data$RT)^2,
                  'invRT' = cor(fitted(invfplmm1), data$invRT)^2,
                  'logRT' = cor(fitted(logfplmm1), data$logRT)^2)


# Check normality of residuals
par(mfrow=c(1,3))
qqnorm(resid(fplmm1),
       main="Normal q-qplot fplmm1")
qqnorm(resid(invfplmm1),
       main="Normal q-qplot invfplmm1")
qqnorm(resid(logfplmm1),
       main="Normal q-qplot logfplmm1")

par(graphical_defaults)

# All models show departures from normality, although this is slightly less for log RT

# Plot residuals
par(mfrow=c(3,1))
plot(resid(fplmm1), fitted(fplmm1),
     main="Residuals fplmm1")
plot(resid(invfplmm1), fitted(invfplmm1),
     main="Residuals invfplmm1")
plot(resid(logfplmm1), fitted(logfplmm1),
     main="Residuals logfplmm1")

par(graphical_defaults)
# There are no clear correlations

# Residual histograms
grid.arrange(hist_resid(fplmm1, 'RT'),
             hist_resid(invfplmm1, '1/RT'),
             hist_resid(logfplmm1, 'logRT'),
             ncol=1)

# All appear to be relatively normally distributed, with a slight positive skew

# Fit models with RT and inverse RT without trimming
trimfplmm1 <- mixed(formula = RT ~ foreperiod*condition*oneBackFP + 
                      (1|participant),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")


# Now we run the same model with inverse RT and logRT as outcomes
triminvfplmm1 <- mixed(formula = invRT ~ foreperiod*condition*oneBackFP + 
                         (1|participant),
                       data = data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

# Now we run the same model with inverse RT and logRT as outcomes
trimlogfplmm1 <- mixed(formula = logRT ~ foreperiod*condition*oneBackFP + 
                         (1|participant),
                       data = data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

isSingular(trimfplmm1)
isSingular(triminvfplmm1)
isSingular(trimlogfplmm1)

# No singular fits here either

# Amount of variance accounted for by the model
var <- rbind(var,
             data.frame(dataset = 'trim',
                        'RT' = cor(fitted(trimfplmm1), data2$RT)^2,
                        'invRT' = cor(fitted(triminvfplmm1), data2$invRT)^2,
                        'logRT'= cor(fitted(trimlogfplmm1), data2$logRT)^2))


var
# Again, the log model accounts for a larger amount of the variance, although this difference is pretty small

# check normality
par(mfrow=c(2,3))
qqnorm(resid(fplmm1),
       main="Normal q-qplot fplmm1")
qqnorm(resid(invfplmm1),
       main="Normal q-qplot invfplmm1")
qqnorm(resid(logfplmm1),
       main="Normal q-qplot logfplmm1")
qqnorm(resid(trimfplmm1),
       main="Normal q-qplot trimfplmm")
qqnorm(resid(triminvfplmm1),
       main="Normal q-qplot triminvfplmm")
qqnorm(resid(trimlogfplmm1),
       main="Normal q-qplot trimlogfplmm")
par(graphical_defaults)

# The log plots show the least deviations from normality; trimming does seem to improve the fit; if anything, it
# appears to worsen it a bit

# Plot residuals
par(mfrow=c(3,2))
plot(resid(fplmm1), fitted(fplmm1),
     main="Residuals RT")
plot(resid(invfplmm1), fitted(invfplmm1),
     main="Residuals 1/RT")
plot(resid(logfplmm1), fitted(logfplmm1),
     main="Residuals logRT")
plot(resid(trimfplmm1), fitted(trimfplmm1),
     main="Residuals trimmed RT")
plot(resid(triminvfplmm1), fitted(triminvfplmm1),
     main="Residuals trimmed 1/RT")
plot(resid(trimlogfplmm1), fitted(trimlogfplmm1),
     main="Residuals trimmed logRT")
par(graphical_defaults)

# logRT still appears to perform better, and trimming seems to yield less outliers

grid.arrange(hist_resid(fplmm1, 'RT'),
             hist_resid(invfplmm1, '1/RT'),
             hist_resid(logfplmm1, 'logRT'),
             hist_resid(trimfplmm1, 'trimmed RT'),
             hist_resid(triminvfplmm1, 'trimmed 1/RT'),
             hist_resid(trimlogfplmm1, 'trimmed logRT'),
             ncol=2,
             as.table=FALSE)

# logRT with trimming results in the closest to normal residual distribution

R2table <- fitstats(fplmm1, 'RT') %>%
  plyr::join(., fitstats(invfplmm1, '1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(logfplmm1, 'log(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimfplmm1, 'trim RT'), by='stat') %>%
  plyr::join(., fitstats(triminvfplmm1, 'trim 1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimlogfplmm1, 'trim log(RT)'), by='stat') %>%
  kable(digits=4)

R2table

# The model using trimmed logRT performs best according to fit. All R2 are very low 



# Variance explained is higher with trimming, although this is probably not significant
# Q-q plots may be worse trimming
# Residuals are less correlated with fitted values with trimming
# Residuals are more normally distributed with trimming

# The model with log RT performs better than the others

#================================= 2.2. Find random effects structure ==========================

#=========================== 2.2.1. FP and FP n-1 as numerical ============================
trimlogfplmm1 <- buildmer(logRT ~ foreperiod * condition * oneBackFP + 
                            (1+foreperiod*condition*oneBackFP|participant), 
                          data=data2,
                          buildmerControl = list(crit='LRT',#ddf = "Satterthwaite",
                                                 family=gaussian(link = 'identity'),
                                                 calc.anova = TRUE))

isSingular(trimlogfplmm1)
formula(trimlogfplmm1)
summary(trimlogfplmm1)
# Systematic comparisons betweeen lmm's via BIC

# Model obtained with buildmer
trimlogfplmm1 <- mixed(logRT ~ 1 + foreperiod + condition + oneBackFP + foreperiod:condition + 
                         foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                         (1 + foreperiod + condition | participant),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

isSingular(trimlogfplmm1)

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
trimlogfplmm <- mixed(logRT ~ 1 + foreperiod + condition + oneBackFP + foreperiod:condition + 
                        foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                        (1 + foreperiod + condition | participant),
                      data=data2,
                      control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = FALSE,
                      method =  'KR',
                      REML=TRUE,
                      return = "merMod")

anova(trimlogfplmm)

#Visualize random effects
dotplot(ranef(trimlogfplmm, condVar = TRUE))

# Visualize interactions
emmip(trimlogfplmm,
      oneBackFP ~ condition, style = "factor") # Using averaged FP

# Single slopes tests
fp_by_condition <- slopes(trimlogfplmm, by = "condition", variables = "foreperiod",
                          p_adjust = "holm")

test(emtrends(trimlogfplmm, ~ condition, var="scaledNumForeperiod")) # equivalent to slopes

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

#============= 3.2. Using RT instead of logRT ================
trimfplmm3 <- mixed(RT ~ 1 + foreperiod + condition + oneBackFP + foreperiod:condition + 
                      foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                      (1 + foreperiod + condition | participant),
                    data=data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod",
                    check_contrasts = FALSE)

isSingular(trimfplmm3)
summary(trimfplmm3)
anova(trimfplmm3)

# 3.3.2.1. Pairwise comparisons by FP (estimates consecutive differences)
emmip(trimfplmm3, condition ~ oneBackFP|foreperiod, CIs = TRUE, style = "factor",
      xlab = "FP n-1") +
  labs(title = "RT pairwise comparisons") +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("orange", "blue"))
ggsave("./Analysis/Plots/mixed_model_pairwise.png",
       width = 8.5,
       height = 5.7)

trimfplmm3emm <- emmeans(trimfplmm3, ~ oneBackFP * condition|foreperiod)
trimfplmm3emm <- emmeans(trimfplmm3, pairwise ~ oneBackFP * condition|foreperiod)

contrast(trimfplmm3emm[[1]], interaction = c("consec", "consec"), by = "foreperiod", adjust = "mvt")


# 3.3.2.2. Pairwise comparisons by FP n-1 (estimate slopes)
emmip(trimfplmm3, condition ~ foreperiod|oneBackFP, CIs = TRUE)

trimfplmm3emm <- emmeans(trimfplmm3, ~ foreperiod * condition|oneBackFP)
trimfplmm3emm <- emmeans(trimfplmm3, pairwise ~ foreperiod * condition|oneBackFP)

contrast(trimfplmm3emm[[1]], interaction = c("poly", "consec"), by = "oneBackFP", adjust = "mvt")

# 3.3.2.3. Separate models for each level of FP n
seq1000 <- mixed(RT ~ 1 + oneBackFP + condition + oneBackFP:condition +
                   (1 + condition | ID),
                 data=filter(data2, foreperiod == '1000'),
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'KR',
                 REML=TRUE,
                 return = "merMod",
                 check_contrasts = FALSE)

summary(seq1000)

anova(seq1000)

emmeans(seq1000, ~ condition|oneBackFP)
emmip(seq1000, condition ~ oneBackFP)


seq1600 <- mixed(RT ~ 1 + oneBackFP + condition + oneBackFP:condition +
                   (1 + condition | ID),
                 data=filter(data2, foreperiod == '1600'),
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'KR',
                 REML=TRUE,
                 return = "merMod",
                 check_contrasts = FALSE)

summary(seq1600)
anova(seq1600)

emmeans(seq1600, pairwise ~ condition|foreperiod)

seq2200 <- mixed(RT ~ 1 + oneBackFP + condition + oneBackFP:condition +
                   (1 + condition | ID),
                 data=filter(data2, foreperiod == '2200'),
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'KR',
                 REML=TRUE,
                 return = "merMod",
                 check_contrasts = FALSE)

summary(seq2200)
anova(seq2200)

emmeans(seq2200, pairwise ~ condition|foreperiod)

seq2800 <- mixed(RT ~ 1 + oneBackFP + condition + oneBackFP:condition +
                   (1 + condition | ID),
                 data=filter(data2, foreperiod == '2800'),
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'KR',
                 REML=TRUE,
                 return = "merMod",
                 check_contrasts = FALSE)

summary(seq2800)
anova(seq2800)

emmeans(seq2800, pairwise ~ condition|foreperiod)

#================== 3.4. Compare dependent variables using random-effects structure ==================
fplmm1 <- mixed(formula = RT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                  numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                  (1 + condition | ID),
                data = data,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'KR',
                REML=TRUE,
                return = "merMod")


invfplmm1 <- mixed(formula = invRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                     numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                     (1 + condition | ID),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'KR',
                   REML=TRUE,
                   return = "merMod")

logfplmm1 <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                     numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                     (1 + condition | ID),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'KR',
                   REML=TRUE,
                   return = "merMod")

trimfplmm1 <- mixed(formula = RT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                      numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                      (1 + condition | ID),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod")


triminvfplmm1 <- mixed(formula = invRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + condition | ID),
                       data = data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

trimlogfplmm1 <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + condition | ID),
                       data = data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

# Compare model R2 and residuals
# Amount of variance accounted for by the model
var <- data.frame(dataset = 'no trim',
                  'RT' = cor(fitted(fplmm1), data$RT)^2,
                  'invRT' = cor(fitted(invfplmm1), data$invRT)^2,
                  'logRT' = cor(fitted(logfplmm1), data$logRT)^2)

var <- rbind(var,
             data.frame(dataset = 'trim',
                        'RT' = cor(fitted(trimfplmm1), data2$RT)^2,
                        'invRT' = cor(fitted(triminvfplmm1), data2$invRT)^2,
                        'logRT'= cor(fitted(trimlogfplmm1), data2$logRT)^2))

var
# Again, the log model accounts (barely) for a larger amount of the variance

# check normality
par(mfrow=c(2,3))
qqnorm(resid(fplmm1),
       main="Normal q-qplot fplmm1")
qqnorm(resid(invfplmm1),
       main="Normal q-qplot invfplmm1")
qqnorm(resid(logfplmm1),
       main="Normal q-qplot logfplmm1")
qqnorm(resid(trimfplmm1),
       main="Normal q-qplot trimfplmm")
qqnorm(resid(triminvfplmm1),
       main="Normal q-qplot triminvfplmm")
qqnorm(resid(trimlogfplmm1),
       main="Normal q-qplot trimlogfplmm")
par(graphical_defaults)

# The log plots show the least deviations from normality, and trimming does seem to help a bit

# Plot residuals
par(mfrow=c(3,2))
plot(resid(fplmm1), fitted(fplmm1),
     main="Residuals RT")
plot(resid(invfplmm1), fitted(invfplmm1),
     main="Residuals 1/RT")
plot(resid(logfplmm1), fitted(logfplmm1),
     main="Residuals logRT")
plot(resid(trimfplmm1), fitted(trimfplmm1),
     main="Residuals trimmed RT")
plot(resid(triminvfplmm1), fitted(triminvfplmm1),
     main="Residuals trimmed 1/RT")
plot(resid(trimlogfplmm1), fitted(trimlogfplmm1),
     main="Residuals trimmed logRT")
par(graphical_defaults)

grid.arrange(hist_resid(fplmm1, 'RT'),
             hist_resid(invfplmm1, '1/RT'),
             hist_resid(logfplmm1, 'logRT'),
             hist_resid(trimfplmm1, 'trimmed RT'),
             hist_resid(triminvfplmm1, 'trimmed 1/RT'),
             hist_resid(trimlogfplmm1, 'trimmed logRT'),
             ncol=2,
             as.table=FALSE)


R2table <- fitstats(fplmm1, 'RT') %>%
  plyr::join(., fitstats(invfplmm1, '1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(logfplmm1, 'log(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimfplmm1, 'trim RT'), by='stat') %>%
  plyr::join(., fitstats(triminvfplmm1, 'trim 1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimlogfplmm1, 'trim log(RT)'), by='stat') %>%
  kable(digits=4)


# By most measures, logRT with trimming still yields the best fit

# Compare lm and lmm models
trimlogfplm <- lm(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod +
                    numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP,
                  data = data2)

trimlogfplmmML <- mixed(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                          (1 + condition | ID),
                        data=data2,
                        control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'KR',
                        REML=FALSE,
                        return = "merMod",
                        check_contrasts = FALSE)

trimlogfplmmML <- lmer(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + condition | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       method =  'KR',
                       REML=FALSE)

anova(trimlogfplmmML, trimlogfplm)

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

