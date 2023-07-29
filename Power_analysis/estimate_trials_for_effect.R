
# Load necessary packages
library(tidyverse)
library(magrittr)
library(lattice)
library(afex)
library(emmeans)
library(car)
library(data.table)
library(BayesFactor)
library(bayestestR)
library(plotrix)
library(pwr)
library(Superpower)
library(InteractionPoweR)

# Save defaults
graphical_defaults <- par()
options_defaults <- options() 


#=================================== 1. Functions ==================================

# convert partial eta squared to f2
pEtaToF2 <- function(pEtaSquared) {
  f2 = sqrt((pEtaSquared)/(1 - pEtaSquared))
}


# Function to get ES for main effects from anova using sampled trials 
trSampleES <- function(dataset, ntrials, variable) {
  sampDS <- dataset %>%
    group_by(ID, foreperiod, condition) %>% # Group by ID, FP, condition
    slice_sample(n = ntrials, replace = TRUE) %>% # Sample trials
    ungroup()
  
  summarySampDS <- sampDS %>%
    group_by(ID, foreperiod, condition) %>%
    summarise(meanRT = mean(RT)) %>%
    ungroup()
  
  sampDSAnova <- aov_ez(data = summarySampDS,
                        id = "ID",
                        dv = "meanRT", within = variable,
                        anova_table = list(es = "pes"))
  DSAnovaTable <- sampDSAnova$anova_table
  DSes <- DSAnovaTable[rownames(DSAnovaTable) == variable,"pes"]#[2,5]
}


# Function to compute power from partial eta squared
trPower <- function(dataset, ntrials, variable, dfVar, dfRes) {
  varEta <- trSampleES(dataset, ntrials, variable)
  varf2 <- pEtaToF2(varEta)
  varPower <- pwr.f2.test(u = dfVar, v = dfRes, f2 = varf2)$power
}

# Function to get p-value to compute rate of false positives
trSampleP <- function(dataset, ntrials, variable) {
  sampDS <- dataset %>%
    group_by(ID, foreperiod, condition) %>% # Group by ID, FP, condition
    slice_sample(n = ntrials, replace = TRUE) %>% # Sample trials
    ungroup()
  
  summarySampDS <- sampDS %>%
    group_by(ID, foreperiod, condition) %>%
    summarise(meanRT = mean(RT)) %>%
    ungroup()
  
  sampDSAnova <- aov_ez(data = summarySampDS,
                        id = "ID",
                        dv = "meanRT", within = variable,
                        anova_table = list(es = "pes"))
  DSAnovaTable <- sampDSAnova$anova_table
  DSes <- DSAnovaTable[rownames(DSAnovaTable) == variable,"Pr(>F)"]#[2,5]
}

# Function to extract effect size for interaction from anova using sampled trials
# needed because the way interactions are passed as arguments to aov_ez is peculiar)
trSampleEsInteraction <- function(dataset, ntrials) {
  sampDS <- dataset %>%
    group_by(ID, foreperiod, condition) %>% # Group by ID, FP, condition
    slice_sample(n = ntrials, replace = TRUE) %>% # Sample trials
    ungroup()
  
  summarySampDS <- sampDS %>%
    group_by(ID, foreperiod, condition) %>%
    summarise(meanRT = mean(RT)) %>%
    ungroup()
  
  sampDSAnova <- aov_ez(data = summarySampDS,
                        id = "ID",
                        dv = "meanRT", within = c("foreperiod", "condition"),
                        anova_table = list(es = "pes"))
  DSAnovaTable <- sampDSAnova$anova_table
  DSes <- DSAnovaTable[3,"pes"] #[2,5]
}


# Same for p-value
trSamplePInteraction <- function(dataset, ntrials) {
  sampDS <- dataset %>%
    group_by(ID, foreperiod, condition) %>% # Group by ID, FP, condition
    slice_sample(n = ntrials, replace = TRUE) %>% # Sample trials
    ungroup()
  
  summarySampDS <- sampDS %>%
    group_by(ID, foreperiod, condition) %>%
    summarise(meanRT = mean(RT)) %>%
    ungroup()
  
  sampDSAnova <- aov_ez(data = summarySampDS,
                        id = "ID",
                        dv = "meanRT", within = c("foreperiod", "condition"),
                        anova_table = list(es = "pes"))
  DSAnovaTable <- sampDSAnova$anova_table
  DSes <- DSAnovaTable[3,"Pr(>F)"] #[2,5]
}

#===============================================================================#
#============================== 2. Constant FP effects ==========================
#===============================================================================#

# Load data
data2 <- read_csv('E:/Post-doc_data/Action_foreperiod/Action_fp_con/Analysis/data2.csv')

#========================== 2. Simulate ES with varying number of trials ==========================
ntrList <- 8:16
nSim <- 1000
esByTr <- vector(mode = "list", length = length(ntrList))

options(dplyr.summarise.inform = FALSE)
for(thisTr in 1:length(esByTr)) {
  ntrials <- ntrList[thisTr]
  esList <- vector(mode = "numeric", length = nSim)
  for(iSim in 1:nSim) {
    esList[iSim] <- trSampleES(data2, ntrials, c("foreperiod", "condition"))
  }
  esByTr[[thisTr]] <- esList
}
options(options_defaults)

# Scatterplot
jpeg("./Analysis/Plots/trials_to_effect.jpeg", width = 1200, height = 1000)
par(mfrow = c(3,3))
for(thisTr in 1:length(ntrList)) {
  esToPlot <- esByTr[[thisTr]]
  plot(1:nSim, esToPlot)
}
par(graphical_defaults)
dev.off()

# Histograms
jpeg("./Analysis/Plots/trials_to_effect_2.jpeg", width = 1200, height = 1000)
par(mfrow = c(3,3))
for(thisTr in 1:length(ntrList)) {
  esToPlot <- esByTr[[thisTr]]
  hist(esToPlot, xlim = c(0,0.15), breaks = 20)
}
par(graphical_defaults)
dev.off()

# Lineplots of means and CIs
jpeg("./Analysis/Plots/trials_to_effect_3.jpeg", width = 1200, height = 1000)
esMeans <- sapply(esByTr, mean)
esUpCIs <- sapply(esByTr, function(x) {
  mean(x) + 1.96 * (sd(x)/sqrt(length(x)))
})
esLowCIs <- sapply(esByTr, function(x) {
  mean(x) - 1.96 * (sd(x)/sqrt(length(x)))
})
plotCI(ntrList, esMeans, ui = esUpCIs, li = esLowCIs)
dev.off()

#====================== 3. Simulate power =============================
ntrList <- 8:16
nSim <- 1000

# Values for pwr function
totalN = 36
dfCondition = 1
dfFP = 1
dfRes = totalN - 1 - dfCondition - dfFP

# Foreperiod
FPEsByTr <- vector(mode = "list", length = length(ntrList))

options(dplyr.summarise.inform = FALSE)
for(thisTr in 1:length(FPEsByTr)) {
  nTrials <- ntrList[thisTr]
  EsList <- vector(mode = "numeric", length = nSim)
  for(iSim in 1:nSim) {
    EsList[iSim] <- trPower(data2, nTrials, "foreperiod", dfFP, dfRes)
  }
  FPEsByTr[[thisTr]] <- EsList
}
options(options_defaults)

# Scatterplot
par(mfrow = c(3,3))
for(thisTr in 1:length(ntrList)) {
  EsToPlot <- FPEsByTr[[thisTr]]
  plot(1:nSim, EsToPlot,
       ylim = c(0, 1),
       main = paste(ntrList[thisTr], "trials"))
  abline(h = 0.05, lty = 3)
}
par(graphical_defaults)

# Plot mean power
powerFPMeans <- sapply(FPEsByTr, mean)
powerFPUpCIs <- sapply(FPEsByTr, function(x) {
  mean(x) + 1.96 * (sd(x)/sqrt(length(x)))
})
powerFPLowCIs <- sapply(FPEsByTr, function(x) {
  mean(x) - 1.96 * (sd(x)/sqrt(length(x)))
})
jpeg("./Power_analysis/mean_power_foreperiod_con.jpeg", width = 1200, height = 1000)
plotCI(ntrList, powerFPMeans, ui = powerFPUpCIs, li = powerFPLowCIs,
       xlab = "Number of trials",
       ylab = "Mean power",
       main = "Statistical power by n of trials in cell")
dev.off()

# Plot of high betas
FPbetas <- lapply(FPEsByTr, function(x) {1-x})

FPhighbetas <- sapply(FPbetas, function(x) {length(x[x > 0.1])/length(x)})

jpeg("./Power_analysis/high_betas_foreperiod_con.jpeg", width = 1200, height = 1000)
plot(ntrList, FPhighbetas, pch = 16,
     ylim = c(0, 0.01),
     main = "Proportion of beta > 0.1 for foreperiod effect",
     xlab = "Number of trials per cell",
     ylab = "Proportion of high betas")
lines(ntrList, FPhighbetas)
dev.off()

#======================== 4. Simulate false positive rate ===============

#========== 4.1. Condition ============
condPByTr <- vector(mode = "list", length = length(ntrList))

options(dplyr.summarise.inform = FALSE)
for(thisTr in 1:length(condPByTr)) {
  nTrials <- ntrList[thisTr]
  PList <- vector(mode = "numeric", length = nSim)
  for(iSim in 1:nSim) {
    PList[iSim] <- trSampleP(data2, nTrials, "condition")
  }
  condPByTr[[thisTr]] <- PList
}
options(options_defaults)

# Scatterplot
par(mfrow = c(3,3))
for(thisTr in 1:length(ntrList)) {
  pToPlot <- condPByTr[[thisTr]]
  plot(1:nSim, pToPlot,
       ylim = c(0, 1),
       main = paste(ntrList[thisTr], "trials"))
  abline(h = 0.05, lty = 3)
}
par(graphical_defaults)

# Plot of false positive rates for condition effects
condFalsePosRate <- sapply(condPByTr, function(x) {length(x[x < 0.05])/length(x)})

options(scipen = 999)
jpeg("./Power_analysis/false_positive_rates_condition_con.jpeg", width = 1200, height = 1000)
plot(ntrList, condFalsePosRate, pch = 16,
     ylim = c(0, 0.0015),
     main = "False positive proportions for condition differences",
     xlab = "Number of trials per cell",
     ylab = "Proportion of positive results")
lines(ntrList, condFalsePosRate)
dev.off()
options(options_defaults)

#========== 4.2. Interaction ============
interactPByTr <- vector(mode = "list", length = length(ntrList))

options(dplyr.summarise.inform = FALSE)
for(thisTr in 1:length(interactPByTr)) {
  ntrials <- ntrList[thisTr]
  PList <- vector(mode = "numeric", length = nSim)
  for(iSim in 1:nSim) {
    PList[iSim] <- trSamplePInteraction(data2, ntrials)
  }
  interactPByTr[[thisTr]] <- PList
}
options(options_defaults)

# Scatterplot
par(mfrow = c(3,3))
for(thisTr in 1:length(ntrList)) {
  pToPlot <- interactPByTr[[thisTr]]
  plot(1:nSim, pToPlot,
       ylim = c(0, 1),
       main = paste(ntrList[thisTr], "trials"))
  abline(h = 0.05, lty = 3)
}
par(graphical_defaults)

# Plot of false positive rates for interaction
interactFalsePosRate <- sapply(interactPByTr, function(x) {length(x[x < 0.05])/length(x)})

jpeg("./Power_analysis/false_positive_proportions_interaction_con.jpeg", width = 1200, height = 1000)
plot(ntrList, interactFalsePosRate, pch = 16,
     ylim = c(0, 0.01),
     main = "False positive proportion for interaction",
     xlab = "Number of trials per cell",
     ylab = "Proportion of positive results")
lines(ntrList, interactFalsePosRate)
dev.off()


#=============================================================================================#
#=================================== 3. Variable FP effect ====================================
#=============================================================================================#

# Load data and keep only fps = 1000 or 2800 ms
data2var <- read_csv('E:/Post-doc_data/Action_foreperiod/Action_foreperiod_exp0_v2/Analysis/data2.csv') %>%
  filter(foreperiod %in% c("1000", "2800"))

# Check remaining n of trials
data2var %>%
  group_by(ID, foreperiod, condition) %>%
  summarise(n())

#========================== 2. Simulate ES with varying number of trials ==========================
ntrList <- 14:30#8:16
nSim <- 1000
esByTr <- vector(mode = "list", length = length(ntrList))

options(dplyr.summarise.inform = FALSE)
for(thisTr in 1:length(esByTr)) {
  ntrials <- ntrList[thisTr]
  esList <- vector(mode = "numeric", length = nSim)
  for(iSim in 1:nSim) {
    esList[iSim] <- trSampleES(data2var, ntrials, c("foreperiod", "condition"))
  }
  esByTr[[thisTr]] <- esList
}
options(options_defaults)

# Scatterplot
jpeg("./Analysis/Plots/trials_to_var_effect.jpeg", width = 1200, height = 1000)
par(mfrow = c(4,5))
for(thisTr in 1:length(ntrList)) {
  esToPlot <- esByTr[[thisTr]]
  plot(1:nSim, esToPlot,
       main = paste(ntrList[thisTr], "trials", sep = " "),
       xlab = "Sim",
       ylab = "ES")
}
par(graphical_defaults)
dev.off()

# Histograms
jpeg("./Analysis/Plots/trials_to_var_effect_2.jpeg", width = 1200, height = 1000)
par(mfrow = c(4,5))
for(thisTr in 1:length(ntrList)) {
  esToPlot <- esByTr[[thisTr]]
  hist(esToPlot,
       xlim = c(0.1,0.8),
       breaks = 20,
       main = paste(ntrList[thisTr], "trials", sep = " "),
       xlab = "ES")
}
par(graphical_defaults)
dev.off()

# Lineplots of means and CIs
jpeg("./Analysis/Plots/trials_to_var_effect_3.jpeg", width = 1200, height = 1000)
esMeans <- sapply(esByTr, mean)
esUpCIs <- sapply(esByTr, function(x) {
  mean(x) + 1.96 * (sd(x)/sqrt(length(x)))
})
esLowCIs <- sapply(esByTr, function(x) {
  mean(x) - 1.96 * (sd(x)/sqrt(length(x)))
})
plotCI(ntrList, esMeans, ui = esUpCIs, li = esLowCIs,
       xlab = "Number of trials",
       ylab = "Mean effect size",
       main = "Effect size by n of trials in cell")
dev.off()

#====================== 3.3. Simulate power =============================
ntrList <- 14:30
nSim <- 1000

# Values for pwr function
totalN = 36
dfCondition = 1
dfFP = 1
dfRes = totalN - 1 - dfCondition - dfFP

#============= 3.3.1. Foreperiod ==============
FPEsByTr <- vector(mode = "list", length = length(ntrList))

options(dplyr.summarise.inform = FALSE)
for(thisTr in 1:length(FPEsByTr)) {
  nTrials <- ntrList[thisTr]
  EsList <- vector(mode = "numeric", length = nSim)
  for(iSim in 1:nSim) {
    EsList[iSim] <- trPower(data2var, nTrials, "foreperiod", dfFP, dfRes)
  }
  FPEsByTr[[thisTr]] <- EsList
}
options(options_defaults)

# Scatterplot
par(mfrow = c(4,5))
for(thisTr in 1:length(ntrList)) {
  EsToPlot <- FPEsByTr[[thisTr]]
  plot(1:nSim, EsToPlot,
       ylim = c(0, 1),
       main = paste(ntrList[thisTr], "trials"))
  abline(h = 0.05, lty = 3)
}
par(graphical_defaults)

# Plot mean power
powerFPMeans <- sapply(FPEsByTr, mean)
powerFPUpCIs <- sapply(FPEsByTr, function(x) {
  mean(x) + 1.96 * (sd(x)/sqrt(length(x)))
})
powerFPLowCIs <- sapply(FPEsByTr, function(x) {
  mean(x) - 1.96 * (sd(x)/sqrt(length(x)))
})
jpeg("./Power_analysis/mean_power_foreperiod_con.jpeg", width = 1200, height = 1000)
plotCI(ntrList, powerFPMeans, ui = powerFPUpCIs, li = powerFPLowCIs,
       xlab = "Number of trials",
       ylab = "Mean power",
       main = "Statistical power for foreperiod effect by n of trials ")
dev.off()

# Plot of high betas
FPbetas <- lapply(FPEsByTr, function(x) {1-x})

FPhighbetas <- sapply(FPbetas, function(x) {length(x[x > 0.1])/length(x)})

jpeg("./Power_analysis/high_betas_foreperiod_con.jpeg", width = 1200, height = 1000)
plot(ntrList, FPhighbetas, pch = 16,
     ylim = c(0, 0.01),
     main = "Proportion of beta > 0.1 for foreperiod effect",
     xlab = "Number of trials per cell",
     ylab = "Proportion of high betas")
lines(ntrList, FPhighbetas)
dev.off()


#========== 3.3.2. Condition ============
condEsByTr <- vector(mode = "list", length = length(ntrList))

options(dplyr.summarise.inform = FALSE)
for(thisTr in 1:length(condPByTr)) {
  nTrials <- ntrList[thisTr]
  PList <- vector(mode = "numeric", length = nSim)
  for(iSim in 1:nSim) {
    EsList[iSim] <- trPower(data2var, nTrials, "condition", dfFP, dfRes)
  }
  condEsByTr[[thisTr]] <- EsList
}
options(options_defaults)

# Scatterplot
par(mfrow = c(4,5))
for(thisTr in 1:length(ntrList)) {
  esToPlot <- condEsByTr[[thisTr]]
  plot(1:nSim, esToPlot,
       ylim = c(0, 1),
       main = paste(ntrList[thisTr], "trials"))
  abline(h = 0.05, lty = 3)
}
par(graphical_defaults)

# Plot mean power
powerCondMeans <- sapply(condEsByTr, mean)
powerCondUpCIs <- sapply(condEsByTr, function(x) {
  mean(x) + 1.96 * (sd(x)/sqrt(length(x)))
})
powerCondLowCIs <- sapply(condEsByTr, function(x) {
  mean(x) - 1.96 * (sd(x)/sqrt(length(x)))
})
jpeg("./Power_analysis/mean_power_condition_var.jpeg", width = 1200, height = 1000)
plotCI(ntrList, powerCondMeans, ui = powerCondUpCIs, li = powerCondLowCIs,
       xlab = "Number of trials",
       ylab = "Mean power",
       main = "Statistical power for condition effect by n of trials")
dev.off()

# Plot of high betas
condbetas <- lapply(condEsByTr, function(x) {1-x})

condhighbetas <- sapply(condbetas, function(x) {length(x[x > 0.1])/length(x)})

jpeg("./Power_analysis/high_betas_condition_var.jpeg", width = 1200, height = 1000)
plot(ntrList, condhighbetas, pch = 16,
     ylim = c(0, 0.10),
     main = "Proportion of beta > 0.1 for condition effect",
     xlab = "Number of trials per cell",
     ylab = "Proportion of high betas")
lines(ntrList, condhighbetas)
dev.off()

#========== 3.3.3 Interaction ============
interactEsByTr <- vector(mode = "list", length = length(ntrList))

options(dplyr.summarise.inform = FALSE)
for(thisTr in 1:length(interactEsByTr)) {
  ntrials <- ntrList[thisTr]
  esList <- vector(mode = "numeric", length = nSim)
  for(iSim in 1:nSim) {
    EsList[iSim] <- trSampleEsInteraction(data2var, ntrials)
  }
  interactEsByTr[[thisTr]] <- EsList
}
options(options_defaults)

# Scatterplot
par(mfrow = c(4,5))
for(thisTr in 1:length(ntrList)) {
  esToPlot <- interactEsByTr[[thisTr]]
  plot(1:nSim, esToPlot,
       ylim = c(0, 1),
       main = paste(ntrList[thisTr], "trials"))
  abline(h = 0.05, lty = 3)
}
par(graphical_defaults)

# Plot mean power
powerInteractMeans <- sapply(interactEsByTr, mean)
powerInteractUpCIs <- sapply(interactEsByTr, function(x) {
  mean(x) + 1.96 * (sd(x)/sqrt(length(x)))
})
powerInteractLowCIs <- sapply(interactEsByTr, function(x) {
  mean(x) - 1.96 * (sd(x)/sqrt(length(x)))
})
jpeg("./Power_analysis/mean_power_interact_var.jpeg", width = 1200, height = 1000)
plotCI(ntrList, powerInteractMeans, ui = powerInteractUpCIs, li = powerInteractLowCIs,
       xlab = "Number of trials",
       ylab = "Mean power",
       main = "Statistical power for interaction by n of trials")
dev.off()

# Plot of high betas
interactbetas <- lapply(interactEsByTr, function(x) {1-x})

interacthighbetas <- sapply(interactbetas, function(x) {length(x[x > 0.1])/length(x)})

jpeg("./Power_analysis/high_betas_condition_var.jpeg", width = 1200, height = 1000)
plot(ntrList, interacthighbetas, pch = 16,
     #ylim = c(0, 0.01),
     main = "Proportion of beta > 0.1 for interaction",
     xlab = "Number of trials per cell",
     ylab = "Proportion of high betas")
lines(ntrList, interacthighbetas)
dev.off()
