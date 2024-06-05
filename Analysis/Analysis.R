
# Load packages

# Read and process data
library(tidyverse)
library(broom)
library(magrittr)
library(data.table)

# Plotting
library(lattice)
library(gtable)
library(gridExtra)
library(gridGraphics)
library(ggdist)
library(ggpubr)
library(extrafont)
library(egg)

# Linear models
library(car)
library(codingMatrices)
library(modelr)
library(afex)
library(emmeans)
library(rtdists)

# Mixed modeling
library(lme4)
library(performance)

# Bayesian analysis
library(BayesFactor)
library(bayestestR)


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

#==========================================================================================#
#======================================= 1. Data quality ===================================
#==========================================================================================#
# 1.1.1. RTs across blocks (conditions aggregated)
ggplot(data=summaryData,
       aes(x=block,
           y=meanRT))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',linewidth=1,aes(group=1))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  labs(title='RT by block')

# 1.1.2. RTs across blocks (separated by condition)
ggplot(data=summaryData,
       aes(x=block,
           y=meanRT,
           color=condition))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',linewidth=1,aes(group=condition))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  cond_cols +
  labs(title='RT by block and condition')

# 1.1.3. Plot RT by foreperiod by participant to check for strange patterns
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT)) +
  stat_summary(fun = 'mean', geom = 'point') +
  stat_summary(fun = 'mean', geom = 'line', aes(group = 1)) +
  stat_summary(fun.data = 'mean_cl_boot', width = 0.2, geom = 'errorbar') +
  facet_wrap(~participant)

# Plot data by foreperiod only to ascertain that there is an FP effect
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           group = 1)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar")

# Check for influence of external fixation duration
ggplot(data=filter(data,condition=='external'),
       aes(x=extFixDur,
           y=RT,
           color=FP))+
  geom_jitter() +
  labs(x = "External fixation duration",
       y = "RT") +
  facet_wrap(~foreperiod)
ggsave("./Analysis/Plots/extfixduration.jpg",
       width = 13.4,
       height = 10)

# Check for influence of latency of action key press on RT
ggplot(data=filter(data,condition=='action'),
       aes(x=actionTrigLatency,
           y=RT,
           color=FP))+
  geom_point() +
  labs(x = "Action trigger delay",
       y = "RT") +
  facet_wrap(~foreperiod)
ggsave("./Analysis/Plots/actiontrigpress.png",
       width = 13.4,
       height = 10)

# Histogram by foreperiod
ggplot(data=filter(data,condition=='action'),
       aes(x=actionTrigLatency,
           fill = FP))+
  geom_histogram(bins = 60) +
  labs(x = "Action trigger delay") +
  facet_wrap(~foreperiod)


# Assess delay between action and WS onset
ggplot(data = data,
                    aes(x = delay)) +
  geom_histogram(bins = 50)

# Correlation between delay and RT
ggplot(data = data2,
       aes(x = delay,
           y = RT)) +
  geom_point()

#==========================================================================================#
#================================= 1.2. Stopping-rule ======================================
#==========================================================================================#

# Sequential plots adding participants to check when curves start to stabilize

plotsList <- list()

for(part in 2:length(unique(data2$participant))) {
  parts <- unique(data2$participant)[1:part]
  plotData <- summaryData2 %>%
    filter(participant %in% parts)
  
  if(part == length(unique(data2$participant))){
    thisPlot <- ggplot(data = plotData,
                       aes(x = foreperiod,
                           y = meanRT,
                           color = condition)) +
      stat_summary(fun = "mean", geom = "point") +
      stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
      stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
      labs(title = paste(c(part, "parts"), collapse = " "),
           x = "Foreperiod",
           y = "Mean RT") +
      cond_cols +
      theme(legend.position = "none")
  } else {
    thisPlot <- ggplot(data = plotData,
                       aes(x = foreperiod,
                           y = meanRT,
                           color = condition)) +
      stat_summary(fun = "mean", geom = "point") +
      stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
      stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
      labs(title = paste(c(part, "parts"), collapse = " "),
           x = "Foreperiod",
           y = "Mean RT") +
      cond_cols +
      theme(legend.position = "none") 
  }
  plotsList[[part-1]] <- thisPlot
}

stopRulePlot <- cowplot::plot_grid(plotlist = plotsList)

#====================== 1.2.1. Individual linear models comparison =========================

# Variables used as predictors: numForeperiod and numOneBackFP
# Dependent variable: logRT
# Variables nested by condition and ID

buildmodel <- function(data) {
  lm(logRT ~ foreperiod*oneBackFP,
     data = data)
}

nested_data <- data2 %>%
  select(participant, condition, foreperiod, oneBackFP, logRT) %>%
  group_by(participant, condition) %>%
  nest()

fitted_data <- nested_data %>%
  mutate(fit = map(data, buildmodel),
         params = map(fit, tidy)) %>%
  ungroup() %>%
  unnest(c(params)) %>%
  select(participant, condition, term, estimate) %>%
  pivot_wider(names_from = term,
              values_from = estimate) %>%
  rename(foreperiod = foreperiod2.8,
         oneBackFP = oneBackFP2.8,
         `foreperiod:oneBackFP` = `foreperiod2.8:oneBackFP2.8`)
  

# Foreperiod
fp_bfs <- ttestBF(x = fitted_data$foreperiod[fitted_data$condition=='external'],
        y = fitted_data$foreperiod[fitted_data$condition=='action'],
        paired=TRUE)

onebackfp_bfs <- ttestBF(x = fitted_data$oneBackFP[fitted_data$condition=='external'],
        y = fitted_data$oneBackFP[fitted_data$condition=='action'],
        paired=TRUE)

interact_bfs <- ttestBF(x = fitted_data$`foreperiod:oneBackFP`[fitted_data$condition=='external'],
        y = fitted_data$`foreperiod:oneBackFP`[fitted_data$condition=='action'],
        paired=TRUE)


#============================ 1.2.2. Mixed models BF comparison ============================ 
library(afex)
library(lme4)
library(buildmer)

with_onebackfp <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                        numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                        (1 + condition + numForeperiod | ID),
                      data = data2,
                      control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = TRUE,
                      method = 'S',
                      REML = TRUE,
                      return = 'merMod')

isSingular(no_onebackfp)

no_onebackfp <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod +
                          (1 + condition + numForeperiod | ID),
                        data = data2,
                        control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method = 'S',
                        REML = TRUE,
                        return = 'merMod')


isSingular(with_onebackfp)

BIC(no_onebackfp, with_onebackfp)

bic_to_bf(c(BIC(no_onebackfp),
            BIC(with_onebackfp)),
          denominator = c(BIC(no_onebackfp)))



with_condition <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                          (1 + condition + numForeperiod | ID),
                        data = data2,
                        control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method = 'S',
                        REML = TRUE,
                        return = 'merMod')

isSingular(with_condition)

no_condition <- mixed(formula = logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP +
                      (1 + numForeperiod | ID),
                      data = data2,
                      control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = TRUE,
                      method = 'S',
                      REML = TRUE,
                      return = 'merMod')

isSingular(no_condition)

bic_to_bf(c(BIC(no_condition),
            BIC(with_condition)),
          denominator = c(BIC(no_condition)))

#============================== 1.2.3. Sequential bayes factors ===========================

# Compute BFs for difference in effects between conditions as participants are added

external_fits <- fitted_data[fitted_data$condition=='external',]
action_fits <- fitted_data[fitted_data$condition=='action',]

srange <- 10:nrow(external_fits)

# FP n main effect
fp_bfs <- sapply(srange, function(range) {
  extractBF(ttestBF(x = external_fits$foreperiod[1:range],
                    y = action_fits$foreperiod[1:range],
                    paired=TRUE),
            onlybf = TRUE)
})

plot(srange, fp_bfs)
lines(srange, fp_bfs)

# FP n-1 main effect
onebackfp_bfs <- sapply(srange, function(range) {
  extractBF(ttestBF(x = external_fits$oneBackFP[1:range],
                    y = action_fits$oneBackFP[1:range],
                    paired=TRUE),
            onlybf = TRUE)
})

plot(srange, onebackfp_bfs)
lines(srange, onebackfp_bfs)

# FP x FP n-1 interaction
interact_bfs <- sapply(srange, function(range) {
  extractBF(ttestBF(x = external_fits$`foreperiod:oneBackFP`[1:range],
                    y = action_fits$`foreperiod:oneBackFP`[1:range],
                    paired=TRUE),
            onlybf = TRUE)
})

plot(srange, interact_bfs)
lines(srange, interact_bfs)

#============================== 1.2.4. Mixed models using brms ============================
b_one_back_fp <- brm(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                       numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                       (1 + condition + numForeperiod | ID),
                     data = data2)

b_one_back_fp_full <- brm(formula = logRT ~ numForeperiod * condition * numOneBackFP + 
                            (1+numForeperiod*condition*numOneBackFP|ID),
                     data = data2)


options(digits = 5)
summary(b_one_back_fp)
options(options_defaults)

#==========================================================================================#
#================================ 2. Descriptive analysis ==================================
#==========================================================================================#

#====================== 2.1. Descriptive statistics ===========================
# RT By foreperiod and condition
meandata <- summaryData2 %>%
  group_by(condition, foreperiod) %>%
  summarise(condRT = mean(meanRT),
            varRT = var(meanRT),
            sdRT = sd(meanRT))
meandata  

# RT By foreperiod
meandataFP <- summaryData2 %>%
  group_by(foreperiod) %>%
  summarise(condRT = mean(meanRT),
            varRT = var(meanRT),
            sdRT = sd(meanRT))
meandataFP

# marginal means for errors
margError <- summaryDataAcc %>%
  group_by(foreperiod, condition) %>%
  summarise(condError = mean(errorRate),
            varError = var(errorRate),
            sdError = sd(errorRate))
margError

#====================== 2.2. Plots =======================

#================== 2.2.1. RT =====================
#=========== 2.2.2.1 By FPn and condition =============
# Boxplots
boxplots <- ggplot(data = summaryData2 %>% 
                     group_by(participant, foreperiod, condition) %>% 
                     summarise(meanRT = mean(meanRT)),
                   aes(x=foreperiod,
                       y=meanRT,
                       fill=condition))+
  geom_jitter(height = 0, width = 0.15, size = 3.1, alpha = 0.5, aes(color=condition)) +
  geom_boxplot()+
  cond_fill +
  cond_cols

boxplots


# Distribution of data
dataHists <- ggplot(data=data2,
                    aes(x=RT))+
  geom_histogram()
dataHists


# Individual histograms
indHistograms <- ggplot(data=data2,
                        aes(x=RT))+
  geom_histogram()+
  facet_wrap(~participant)
indHistograms

# QQ plots by participant
qqmath(~RT|participant, data=data2)
qqmath(~invRT|participant, data=data2)


# Do the same by participant, condition and foreperiod
RTstats <- data2 %>%
  group_by(participant, condition, foreperiod) %>%
  summarise(meanRT = mean(RT), sdRT = sd(RT)) %>%
  ungroup()

ggplot(data=RTstats,
       aes(x=meanRT,
           y=sdRT,
           color = condition)) +
  geom_point() +
  geom_smooth(method='lm', aes(fill = condition)) +
  cond_cols +
  cond_fill+
  facet_wrap(~foreperiod) 

View(RTstats[which.max(RTstats$sdRT),])

cor(RTstats_full$meanRT, RTstats_full$sdRT)


# FP and condition
rt_by_condition <- ggplot(data = summaryData2 %>% 
                            group_by(participant, foreperiod, condition) %>% 
                            summarise(meanRT = mean(meanRT)),
                          aes(x = foreperiod,
                              y = meanRT,
                              color = condition)) +

  geom_jitter(height = jthgt, width = jtwdt, size = psz, alpha = alp) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = lsz, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = erlsz, width = erlwth, geom = "errorbar") + 
  labs(title = "RT",
       x = "FP (s)", 
       y = "Mean RT (s)",
       color = "Condition") +
  cond_cols
ggplot2::ggsave("./Analysis/Plots/RT_by_condition.jpg",
                rt_by_condition,
                width = 15,
                height = 5,
                unit = "cm",
                dpi = 300)



#=========== 2.1.2 Sequential effects =============
# Separated by FP n-1
seqEff_by_oneback <- ggplot(data = summaryData2 %>%
                              group_by(participant, foreperiod, condition, oneBackFP) %>%
                              summarise(meanRT = mean(meanRT)),
                            aes(x = foreperiod,
                                y = meanRT,
                                color = condition)) +
  geom_jitter(height = jthgt, width = jtwdt, size = psz, alpha = alp) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = lsz, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = erlsz, width = erlwth, geom = "errorbar") + 
  labs(x = "FP (s)",
       y = "Mean RT (s)",
       color = "Condition") +
  facet_wrap(~oneBackFP,
             labeller = as_labeller(c(`1` = "FP[n-1] == 1.0",
                                      `2.8` = "FP[n-1] == 2.8"),
                                    default = label_parsed)) +
  cond_cols
seqEff_by_oneback <- set_panel_size(seqEff_by_oneback, width = unit(4, "cm"),
                                    height = unit(2.6, "cm"))
ggsave("./Analysis/Plots/SeqEff.jpg",
       seqEff_by_oneback,
       width = 10.5,
       height = 5,
       units = "cm",
       dpi = 300)


#============== 2.2. Accuracy ================
#======= 2.2.1. Mean acc by FP and condition =======

ggplot(data = summaryDataAcc) +
  geom_histogram(aes(x = meanAcc)) +
  facet_grid(foreperiod ~ condition)

error_by_condition <- ggplot(data = summaryDataAcc %>%
                               group_by(participant, foreperiod, condition) %>%
                               summarise(errorRate = mean(errorRate)),
                             aes(x = foreperiod,
                                 y = errorRate,
                                 color = condition)) +
  geom_jitter(height = jthgt, width = jtwdt, size = psz, alpha = alp) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = lsz, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = erlsz, width = erlwth, geom = "errorbar") + 
  labs(x = "FP (s)",
       y = "Mean Error Rate",
       color = "Condition",
       title = "Error Rate") +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5, 1), "pt")) +
  cond_cols
ggsave("./Analysis/Plots/errors_by_fp_condition.jpg",
       error_by_condition,
       width = 16,
       height = 5,
       unit = "cm",
       dpi = 300)

# RT and error rate in single panel
cond_legend <- ggpubr::get_legend(error_by_condition)

xaxis_title <- text_grob(error_by_condition$labels$x, # Extract x axis label to plot in the middle of both plots
                         just = "top",
                         size = mytheme$axis.title$size) # size of grob; 11 is the base size in theme_classic
                         #size = mytheme$axis.title$size * 11) # size of grob; 11 is the base size in theme_classic


xaxis_title_margin <- unit(2/5, "line") # margin to separate from other plot elements


# Visualize
grid.arrange(arrangeGrob(rt_by_condition + theme(legend.position = "none",
                                                 axis.title.x = element_blank()),
                         error_by_condition + theme(legend.position = "none",
                                                    axis.title.x = element_blank()),
                         nrow = 1,
                         widths = c(1/2, 1/2)),
             xaxis_title,
             cond_legend,
             heights = unit.c(unit(1, "null"),
                              grobHeight(xaxis_title) + xaxis_title_margin,
                              grobHeight(cond_legend)),
             nrow = 3)

# Save plots
RT_panel <- rt_by_condition + theme(legend.position = "none",
                                    axis.title.x = element_blank())
RT_panel <- set_panel_size(RT_panel,
                           width = unit(4, "cm"),
                           height = unit(2.6, "cm"))


error_panel <- error_by_condition + theme(legend.position = "none",
                                          axis.title.x = element_blank())
error_panel <- set_panel_size(error_panel, 
                              width = unit(4, "cm"),
                              height = unit(2.6, "cm"))

rt_error_plots <- arrangeGrob(arrangeGrob(RT_panel,
                                          error_panel,
                                          nrow = 1,
                                          widths = c(1/2, 1/2)),
                              xaxis_title,
                              cond_legend,
                              heights = unit.c(unit(1, "null"),
                                               grobHeight(xaxis_title) + xaxis_title_margin,
                                               grobHeight(cond_legend)),
                              nrow = 3)

ggsave("./Analysis/Plots/rt_error_plots.jpg",
       rt_error_plots,
       width = 10.5,
       height = 5,
       unit = "cm",
       dpi = 300)


#==========================================================================================#
#======================================= 3. ANOVAs =========================================
#==========================================================================================#

#==================== 3.1. FP x RT by condition ======================
#================ 3.1.1. RT =================


# Run repeated-measures anova
fpAnova <- aov_ez(id = "participant",
       dv = "meanRT",
       data = summaryData2,
       within = c("foreperiod", "condition"))

### Check assumptions

# Sphericity
check_sphericity(fpAnova)

fpAnova_plot <- afex_plot(fpAnova, x = 'foreperiod', trace = 'condition', error = 'within')

# Normality of residuals
is_norm <- check_normality(fpAnova)

plot(is_norm)

plot(is_norm, type = "qq")
plot(is_norm, type = "qq", detrend = TRUE)

testnormality = function(dfr) return(shapiro.test(dfr$RT)$p.value)
p = as.vector(by(data, data$participant, testnormality))
names(p) = levels(data$participant)
names(p[p < 0.05])


# Try transformations
invfpAnova <- aov_ez(id = "participant",
                  dv = "meanInvRT",
                  data = summaryData2,
                  within = c("foreperiod", "condition"))

### Check assumptions

# Sphericity
check_sphericity(invfpAnova)

# Normality of residuals
is_norm <- check_normality(invfpAnova)

plot(is_norm)

plot(is_norm, type = 'qq')

plot(is_norm, type = 'qq', detrend = TRUE)

# using 1/RT does not solve the problem of non-normality

logfpAnova <- aov_ez(id = "participant",
                     dv = "meanLogRT",
                     data = summaryData2,
                     within = c("foreperiod", "condition"))


### Check assumptions

# Sphericity
check_sphericity(logfpAnova)

# Normality of residuals
is_norm <- check_normality(logfpAnova)

plot(is_norm)

plot(is_norm, type = 'qq')

plot(is_norm, type = 'qq', detrend = TRUE)

# The log-transform does not solve the problem either

#================ 3.1.2. Accuracy =================

AccAnova <- aov_ez(id = "ID",
                      dv = "meanAcc",
                      data = summaryDataAcc,
                      within = c("foreperiod", "condition"))


# ============ 3.1.2. Fit models separately for external condition (safety check) ====================
logfpAnovaExt <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("foreperiod"))

# Sphericity
check_sphericity(logfpAnovaExt)

# Normality of residuals
is_norm <- check_normality(logfpAnovaExt)

plot(is_norm)

plot(is_norm, type = 'qq')

plot(is_norm, type = 'qq', detrend = TRUE)


# Effect is significant

logfpAnovaExt2 <- aov_ez(id = "ID",
                      dv = "meanLogRT",
                      data = summaryData2[summaryData2$condition=='external',],
                      within = c("logFP"))

# Effect is significant

#============================== 3.2. Sequential effects ================================================

#============= 3.2.1. RT ==================
# 2.2.1. Anova for FP n-1
seqEffAnova <- aov_ez(id = "participant",
                  dv = "meanRT",
                  data = summaryData2,
                  within = c("foreperiod", "condition", "oneBackFP"))

logSeqEffAnova <- aov_ez(id = "participant",
                      dv = "meanLogRT",
                      data = summaryData2,
                      within = c("foreperiod", "condition", "oneBackFP"))

nice(seqEffAnova,
     correction='none')
                  
nice(logSeqEffAnova,
     correction = "none")

# 2.2.2. FP n-2
twoBackAnova <- aov_ez(id = "ID",
                       dv = "meanRT",
                       data = summaryData2,
                       within = c("foreperiod", "condition", "twoBackFP"))

nice(twoBackAnova,
     correction = "none")

ntworegression <- lm(meanRT ~ foreperiod * oneBackFP * twoBackFP, 
                     data = summaryData)
summary(ntworegression)
anova(ntworegression)
Anova(ntworegression, type = "II")


#============= 3.2.2. Accuracy ================
AccSeqAnova <- aov_ez(id = "ID",
                   dv = "meanAcc",
                   data = summaryDataAcc,
                   within = c("foreperiod", "condition", "oneBackFP"))

# ============ 3.2.2. Fit models separately for external condition (safety check) ====================
seqEffAnovaExt <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("foreperiod", "oneBackFP"))

# Effect of FP is significant

seqEffAnovaExt2 <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("logFP", "oneBackFP"))

# Effect of FP is significant

# Compare full and reduced models
summary(logfpAnovaExt, correction = 'none')
summary(seqEffAnovaExt, correction = 'none')

summary(logfpAnovaExt2, correction = 'none')
summary(seqEffAnovaExt2, correction = 'none')

#=================== 3.3 Sequential effects with difference between current and previous FP ============

# lm with difference between the durations of FPn and FPn-1 as regressor
fpDiffRegression <- lm(meanRT ~ foreperiod * condition * oneBackFPDiff,
                       data = summaryData)
summary(fpDiffRegression)
anova(fpDiffRegression)


ggplot(data = summaryData2,
       aes(x = oneBackFPDiff,
           y = meanSeqEff)) +
  stat_summary(fun = 'mean', geom = 'point')

ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanSeqEff)) +
  stat_summary(fun = 'mean', geom = 'point')

# ============ 3.3.2. Fit models separately for external condition (safety check) ====================
fpDiffAnovaExt <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("foreperiod", "oneBackFPDiff"))

# Error because there are empty cells

logfpAnova2 <- aov_ez(id = "ID",
                      dv = "meanLogRT",
                      data = summaryData2[summaryData2$condition=='external',],
                      within = c("logFP", "oneBackFPDiff"))

# Error because there are empty cells

#=============================== 3.4 Orientation ===================================

# No apparent difference
oriAnova <- aov_ez(id = 'ID',
                   dv = 'meanRT',
                   data = summaryData2,
                   within = c('orientation', 'condition', 'foreperiod'))

#============================= 4. Quadratic effects  =====================================
# fpEmmeans <- emmeans(fpAnova,
#                      pairwise ~ condition|foreperiod,
#                      adjust = 'none')


fpEmmeans <- emmeans(fpAnova,
                     pairwise ~ foreperiod|condition,
                     adjust = 'bonferroni')

fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
                               interaction=c('poly'),
                               adjust='bonferroni')

fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
                               interaction=c('consec'),
                               adjust='bonferroni')

contrasts(summaryData2$foreperiod) <- contr.poly(4)

# using emmeans
logfpAnova <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("foreperiod"))

fp_posthoc <- emmeans(logfpAnova,
                      specs = "foreperiod")

fp_contrasts <- contrast(fp_posthoc,
                         interaction = 'poly')

# Using lists and anova as aov
logfpAnova <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("foreperiod"),
                     return = 'aov')


summary(logfpAnova, split=list(foreperiod=list(linear=1,quadratic=2,cubic=3)))

# Using regression
logfpregression <- lm(meanLogRT ~ numForeperiod, data = summaryData2[summaryData2$condition=='external',])
logfpregression <- lm(meanLogRT ~ numForeperiod * numOneBackFPDiff, 
                      data = summaryData2[summaryData2$condition=='external',])
summary(logfpregression)
anova(logfpregression)

fp_posthoc <- emmeans(logfpregression,
                      specs = "foreperiod", 
                      adjust="tukey")

fp_contrasts <- contrast(fp_posthoc,
                         interaction = "poly",
                         adjust = "tukey")

#=================================== 5. Learning =============================================
#============================ 5.1. Effects across blocks ==============================
# Foreperiod, condition and block
blocklm <- lm(meanRT ~ foreperiod * counterbalance * block,
              data = summaryData)

anova(blocklm)
Anova(blocklm)



#========================== 5.2. Split anovas by counterbalancing order ===================

fpAnova_ae <- aov_ez(id = "ID",
                     dv = "meanRT",
                     data = summaryData2[summaryData2$counterbalance=='action-external',],
                     within = c("foreperiod", "condition"))

fpAnova_ea <- aov_ez(id = "ID",
                     dv = "meanRT",
                     data = summaryData2[summaryData2$counterbalance=='external-action',],
                     within = c("foreperiod", "condition"))


#===================================================================================================#
#=================================== 6. Analysis with scaled predictors =============================
#===================================================================================================#

fpregression <- lm(meanRT ~ condition * foreperiod, data = summaryData)
summary(fpregression)
anova(fpregression)

logfpregression <- lm(meanRT ~ condition * logFP, data = summaryData)
anova(logfpregression)



# ============ 6.1.2. Fit models separately for external condition (safety check) ====================
logfpAnovaExt <- aov_ez(id = "ID",
                        dv = "meanLogRT",
                        data = summaryData2[summaryData2$condition=='external',],
                        within = c("foreperiod"))

# Effect is significant

logfpAnovaExt2 <- aov_ez(id = "ID",
                         dv = "meanLogRT",
                         data = summaryData2[summaryData2$condition=='external',],
                         within = c("logFP"))

# Effect is significant

#============================== 6.2. Sequential effects ================================================
# 2.2.1. Anova for FP n-1
seqEffAnova <- aov_ez(id = "ID",
                      dv = "meanRT",
                      data = summaryData2,
                      within = c("foreperiod", "condition", "oneBackFP"))

seqEffAnovaInv <- aov_ez(id = "ID",
                      dv = "meanInvRT",
                      data = summaryData,
                      within = c("foreperiod", "condition", "oneBackFP"))

nice(seqEffAnova,
     correction='none')

seqEffAnovaLog <- aov_ez(id = "ID",
                      dv = "meanLogRT",
                      data = summaryData,
                      within = c("foreperiod", "condition", "oneBackFP"))

nice(seqEffAnovaLog,
     correction='none')

seqEFffregression <- lm(meanRT ~ foreperiod * oneBackFP * condition, 
                        data = summaryData)
summary(seqEFffregression)
anova(seqEFffregression)

logseqEffregression <- lm(meanRT~logFP*logoneBackFP*condition,
                          data=summaryData) 
anova(logseqEffregression)

# 2.2.3. Anova for FP n-2
ntworegression <- lm(meanRT ~ foreperiod * oneBackFP * twoBackFP, 
                     data = summaryData)
summary(ntworegression)
anova(ntworegression)
Anova(ntworegression, type = "II")

# ============ 6.2.2. Fit models separately for external condition (safety check) ====================
seqEffAnovaExt <- aov_ez(id = "ID",
                         dv = "meanLogRT",
                         data = summaryData2[summaryData2$condition=='external',],
                         within = c("foreperiod", "oneBackFP"))

# Effect of FP is significant

seqEffAnovaExt2 <- aov_ez(id = "ID",
                          dv = "meanLogRT",
                          data = summaryData2[summaryData2$condition=='external',],
                          within = c("logFP", "oneBackFP"))

# Effect of FP is significant

# Compare full and reduced models
summary(logfpAnovaExt, correction = 'none')
summary(seqEffAnovaExt, correction = 'none')

summary(logfpAnovaExt2, correction = 'none')
summary(seqEffAnovaExt2, correction = 'none')

#=================== 6.3 Sequential effects with difference between current and previous FP ============

# Lm with difference between the durations of FPn and FPn-1 as regressor
fpDiffRegression <- lm(meanRT ~ foreperiod * condition * oneBackFPDiff,
                       data = summaryData)
summary(fpDiffRegression)
anova(fpDiffRegression)


ggplot(data = summaryData2,
       aes(x = oneBackFPDiff,
           y = meanSeqEff)) +
  stat_summary(fun = 'mean', geom = 'point')

ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanSeqEff)) +
  stat_summary(fun = 'mean', geom = 'point')

# ============ 6.3.2. Fit models separately for external condition (safety check) ====================
fpDiffAnovaExt <- aov_ez(id = "ID",
                         dv = "meanLogRT",
                         data = summaryData2[summaryData2$condition=='external',],
                         within = c("foreperiod", "oneBackFPDiff"))

# Error because there are empty cells

logfpAnova2 <- aov_ez(id = "ID",
                      dv = "meanLogRT",
                      data = summaryData2[summaryData2$condition=='external',],
                      within = c("logFP", "oneBackFPDiff"))

# Error because there are empty cells



#============================= 6.5. Quadratic effects  =====================================
# fpEmmeans <- emmeans(fpAnova,
#                      pairwise ~ condition|foreperiod,
#                      adjust = 'none')


fpEmmeans <- emmeans(fpAnova,
                     pairwise ~ foreperiod|condition,
                     adjust = 'bonferroni')

fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
                               interaction=c('poly'),
                               adjust='bonferroni')

fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
                               interaction=c('consec'),
                               adjust='bonferroni')

contrasts(summaryData2$foreperiod) <- contr.poly(4)

# using emmeans
logfpAnova <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("foreperiod"))

fp_posthoc <- emmeans(logfpAnova,
                      specs = "foreperiod")

fp_contrasts <- contrast(fp_posthoc,
                         interaction = 'poly')

# Using lists and anova as aov
logfpAnova <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("foreperiod"),
                     return = 'aov')


summary(logfpAnova, split=list(foreperiod=list(linear=1,quadratic=2,cubic=3)))

# Using regression
logfpregression <- lm(meanLogRT ~ numForeperiod, data = summaryData2[summaryData2$condition=='external',])
logfpregression <- lm(meanLogRT ~ numForeperiod * numOneBackFPDiff, 
                      data = summaryData2[summaryData2$condition=='external',])
summary(logfpregression)
anova(logfpregression)

fp_posthoc <- emmeans(logfpregression,
                      specs = "foreperiod", 
                      adjust="tukey")

fp_contrasts <- contrast(fp_posthoc,
                         interaction = "poly",
                         adjust = "tukey")

#============================= 7. Residual analysis for sequential effects ==========================

seqEffAnovares1 <- aov_ez(id = "ID",
                      dv = "meanLogRT",
                      data = filter(summaryData2, condition == 'external'),
                      within = c("oneBackFP"),
                      fun_aggregate = mean)

seqEffres <- seqEffAnovares1$lm$residuals

seqEffres <- seqEffAnovares1$lm$residuals %>%
  as_tibble() %>%
  pivot_longer(cols = 1:4, names_to = c("oneBackFP"), values_to = "residuals")
seqEffres$oneBackFP <- fct_recode(seqEffres$oneBackFP, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')

aovez_lm <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("foreperiod", "oneBackFP"),
                  fun_aggregate = mean)
  modellm <- model$lm
}

aovez_lm2 <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("oneBackFP"),
                  fun_aggregate = mean)
  modellm <- model$lm
}

seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm(data)))

seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~lm(meanLogRT ~ oneBackFP, data = .x)),
         predictions = map2(data, model, add_predictions))

# Add residuals
seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm2(data)),
         residuals = map2(data, model, add_residuals))

# Add predictions
seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm2(data)),
         predictions = map2(data, model, add_predictions))


seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm(data))) %>%
  unnest(c(data))


# Add predictions using broom
seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm(data)),
         model_data = map(model, broom::augment)) %>%
  unnest(model_data)

aovez_resid <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("foreperiod", "oneBackFP"),
                  fun_aggregate = mean)
  model_res <- model$lm$residuals %>%
    as_tibble() %>%
    pivot_longer(cols = 1:16, names_to = c("foreperiod", "oneBackFP"), names_pattern = "(.*)_(.*)", values_to = "residuals")
  model_res$foreperiod <- fct_recode(model_res$foreperiod, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res$oneBackFP <- fct_recode(model_res$oneBackFP, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res
}

aovez_resid2 <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("oneBackFP"),
                  fun_aggregate = mean)
  model_res <- model$lm$residuals %>%
    as_tibble() %>%
    pivot_longer(cols = 1:4, names_to = c("oneBackFP"), values_to = "residuals")
  model_res$oneBackFP <- fct_recode(model_res$oneBackFP, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res
}

seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm(data)),
         residuals = map(data, ~aovez_resid(.x)))

aggSumData2 <- summaryData2 %>%
  group_by(ID, condition, foreperiod, oneBackFP) %>%
  summarise(meanLogRT = mean(meanLogRT)) %>%
  ungroup()

seqEffResInt <- aggSumData2 %>%
  group_by(condition, foreperiod) %>%
  nest() %>%
  mutate(residuals = map(data, ~aovez_resid2(.x))) %>%
  unnest(residuals) %>%
  ungroup()

seqEffResData <- aggSumData2 %>%
  mutate(residuals = seqEffResInt$residuals)

seqEffAnovares2 <- aov_ez(id = 'ID',
                          dv = 'residuals',
                          data = seqEffResData,
                          within = c('condition'))

summary(seqEffAnovares2)

condAnova <- aov_ez(id = "ID",
                    dv = "meanLogRT",
                    data = summaryData2,
                    within = c("condition"))

#==================================================================================#
# Functions to extract lm and residuals from aov_ez
aovez_lm <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("foreperiod", "oneBackFP"),
                  fun_aggregate = mean)
  modellm <- model$lm
}

aovez_resid <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("foreperiod", "oneBackFP"),
                  fun_aggregate = mean)
  model_res <- model$lm$residuals %>%
    as_tibble() %>%
    pivot_longer(cols = 1:16, names_to = c("foreperiod", "oneBackFP"), names_pattern = "(.*)_(.*)", values_to = "residuals")
  model_res$foreperiod <- fct_recode(model_res$foreperiod, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res$oneBackFP <- fct_recode(model_res$oneBackFP, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res
}

# Aggregate over other cells
aggSumData <- summaryData2 %>%
  group_by(ID, condition, foreperiod, oneBackFP) %>%
  summarise(meanLogRT = mean(meanLogRT)) %>%
  ungroup()

# Save residuals to aggregated dataset
seqEffResInt <- aggSumData %>%
  group_by(condition) %>%
  nest() %>%
  mutate(residuals = map(data, ~aovez_resid(.x))) %>%
  unnest(residuals) %>%
  ungroup()

seqEffResData <- aggSumData %>%
  mutate(residuals = seqEffResInt$residuals)

# Run anova with condition as IV
seqEffAnovares2 <- aov_ez(id = 'ID',
                          dv = 'residuals',
                          data = seqEffResData,
                          within = c('condition'))


condAnova <- aov_ez(id = "ID",
                    dv = "meanLogRT",
                    data = summaryData2,
                    within = c("condition"))


#=====================#
aovez_resid2 <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("oneBackFP"),
                  fun_aggregate = mean)
  model_res <- model$lm$residuals %>%
    as_tibble() %>%
    pivot_longer(cols = 1:4, names_to = c("oneBackFP"), values_to = "residuals")
  model_res$oneBackFP <- fct_recode(model_res$oneBackFP, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res
}

seqEffResData2 <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm(data)),
         residuals = map(data, ~aovez_resid2(.x)))

aggSumData2 <- summaryData2 %>%
  group_by(ID, condition, foreperiod, oneBackFP) %>%
  summarise(meanLogRT = mean(meanLogRT)) %>%
  ungroup()

seqEffResInt2 <- aggSumData2 %>%
  group_by(condition, foreperiod) %>%
  nest() %>%
  mutate(residuals = map(data, ~aovez_resid2(.x))) %>%
  unnest(residuals) %>%
  ungroup()

seqEffResData2 <- aggSumData2 %>%
  mutate(residuals = seqEffResInt2$residuals)

seqEffAnovares2 <- aov_ez(id = 'ID',
                          dv = 'residuals',
                          data = seqEffResData2,
                          within = c('condition'))


condAnova <- aov_ez(id = "ID",
                    dv = "meanLogRT",
                    data = summaryData2,
                    within = c("condition"))

#============================#

aovez_resid3 <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("foreperiod"),
                  fun_aggregate = mean)
  model_res <- model$lm$residuals %>%
    as_tibble() %>%
    pivot_longer(cols = 1:4, names_to = c("foreperiod"), values_to = "residuals")
  model_res$foreperiod <- fct_recode(model_res$foreperiod, 
                                     '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res
}

aggSumData3 <- summaryData2 %>%
  group_by(ID, condition, foreperiod, oneBackFP) %>%
  summarise(meanLogRT = mean(meanLogRT)) %>%
  ungroup()

seqEffResInt3 <- aggSumData3 %>%
  group_by(condition, oneBackFP) %>%
  nest() %>%
  mutate(residuals = map(data, ~aovez_resid3(.x))) %>%
  unnest(residuals) %>%
  ungroup()

seqEffResData3 <- aggSumData3 %>%
  mutate(residuals = seqEffResInt3$residuals)

seqEffAnovares3 <- aov_ez(id = 'ID',
                          dv = 'residuals',
                          data = seqEffResData3,
                          within = c('condition'))

summary(seqEffAnovares4)

