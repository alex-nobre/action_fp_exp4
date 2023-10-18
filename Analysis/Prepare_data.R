

# Load necessary packages
library(readr)
library(ggplot2)
library(magrittr)
library(dplyr)
library(forcats)

source("./Analysis/helper_functions.R")

# Read data
data <- read_csv("./Analysis/dataActionFPAll.csv")

# Coerce to factors
data <- data %>%
  mutate(across(c(participant, block, condition, FPType, conFPDur, orientation, FP, oneBackFP, countBalance, handedness), as_factor))

data$condition <- data$condition %>%
  fct_relevel(c("external", "action"))

# Column for testing position based on counterbalancing
data <- data %>%
  mutate(condTestPos = case_when((condition == 'action' & countBalance %in% c('1', '3')) ~ '1',
                             (condition == 'external' & countBalance %in% c('1', '3')) ~ '2',
                             (condition == 'action' & countBalance %in% c('2', '4')) ~ '2',
                             (condition == 'external' & countBalance %in% c('2', '4')) ~ '1')) %>%
  mutate(condTestPos = as.factor(condTestPos))

# Create numeric versions of foreperiod and FP n-1
data$numFP <- as.numeric(as.character(data$FP))
data$numOneBackFP <- as.numeric(as.character(data$oneBackFP))

# Quadratic term for numForeperiod
data$squaredNumFP <- data$numFP^2

# Create column for previous orientation and for comparison of current and previous orientations
data <- data %>%
  mutate(seqOri = ifelse(lag(orientation)==orientation, 'same', 'different'),
         prevOri = lag(orientation)) %>%
  mutate(seqOri = as.factor(seqOri),
         prevOri = as.factor(prevOri))
  

# Remove trials without n-1 FP values (i.e., first of each block)
data <- data %>%
  filter((FPType == 'variable' & !is.na(oneBackFP)) | (FPType == 'constant'))

# Save data with error trials to assess accuracy
dataAll <- data

# Keep only trials with correct responses to analyze RT
data <- data %>%
  filter(!is.na(RT), Acc == 1)

# Create log10 of continuous indenpendent variables
data$numLogFP <- log10(data$numFP)
data$logFP <- as.factor(data$numLogFP)
data$logOneBackFP <- log10(data$numOneBackFP)

dataAll$numLogFP <- log10(dataAll$numFP)
dataAll$logFP <- as.factor(dataAll$numLogFP)
dataAll$logOneBackFP <- log10(dataAll$numOneBackFP)

# Remove extreme values
data <- data %>%
  filter(RT < 1.0) %>%
  filter(RT > 0.15)

# Transform RT to reduce skew
data$logRT <- ifelse(!is.na(data$RT), log10(data$RT), NA) # log-transform
data$invRT <- ifelse(!is.na(data$RT), 1/data$RT, NA)


# Trimming
data2 <- data %>%
  group_by(participant) %>%
  mutate(RTzscore=ifelse(!is.na(RT), compute_zscore(RT), NA),
         logRTzscore=ifelse(!is.na(RT), compute_zscore(logRT), NA)) %>%
  filter(abs(logRTzscore) < 3) %>%
  ungroup()

# No trimming
data <- data %>%
  group_by(participant) %>%
  mutate(RTzscore=ifelse(!is.na(RT), compute_zscore(RT), NA),
         logRTzscore=ifelse(!is.na(RT), compute_zscore(logRT), NA)) %>%
  #filter(abs(logRTzscore) < 3) %>%
  ungroup()

# Create scaled predictors
data$scaledNumFP <- scale(data$numFP, scale = FALSE)[,1]
data$squaredScaledNumFP <- data$scaledNumFP^2
data$scaledNumOneBackFP <- scale(data$numOneBackFP, scale = FALSE)[,1]

data2$scaledNumFP <- scale(data2$numFP, scale = FALSE)[,1]
data2$squaredScaledNumFP <- data2$scaledNumFP^2
data2$scaledNumOneBackFP <- scale(data2$numOneBackFP, scale = FALSE)[,1]

dataAll$scaledNumFP <- scale(dataAll$numFP, scale = FALSE)[,1]
dataAll$squaredScaledNumFP <- dataAll$scaledNumFP^2
dataAll$scaledNumOneBackFP <- scale(dataAll$numOneBackFP, scale = FALSE)[,1]


# Average data
summaryData <- data %>%
  group_by(participant, FP, logFP, condition, FPType, conFPDur,
           orientation, prevOri, seqOri,
           oneBackFP, block, countBalance, condTestPos) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT)) %>%
  ungroup() %>%
  mutate(numFP = as.numeric(as.character(FP)),
         numOneBackFP = as.numeric(as.character(oneBackFP)),
         numLogFP = as.numeric(as.character(logFP))) %>%
  mutate(squaredNumFP = numFP^2,
         squaredNumLogFP = numLogFP^2,
         scaledNumFP = scale(numFP)[,1],
         squaredScaledNumFP = scaledNumFP^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1])

summaryData2 <- data2 %>%
  group_by(participant, FP, logFP, condition, FPType, conFPDur,
           orientation, prevOri, seqOri,
           oneBackFP, block, countBalance, condTestPos) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT)) %>%
  ungroup() %>%
  mutate(numFP = as.numeric(as.character(FP)),
         numOneBackFP = as.numeric(as.character(oneBackFP)),
         numLogFP = as.numeric(as.character(logFP))) %>%
  mutate(squaredNumFP = numFP^2,
         squaredNumLogFP = numLogFP^2,
         scaledNumFP = scale(numFP)[,1],
         squaredScaledNumFP = scaledNumFP^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1])


summaryDataAll <- dataAll %>%
  group_by(participant,FP,condition,FPType,conFPDur,
           oneBackFP) %>%
  summarise(meanRT = mean(RT),
            meanAcc = mean(Acc)) %>%
  ungroup() %>%
  mutate(numFP = as.numeric(as.character(FP)),
         numOneBackFP = as.numeric(as.character(oneBackFP))) %>%
  mutate(squaredNumFP = numFP^2,
         scaledNumFP = scale(numFP)[,1],
         squaredScaledNumFP = scaledNumFP^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1])

#write_csv(data, "./Analysis/data.csv")
#write_csv(data2, "./Analysis/data2.csv")
#write_csv(summaryData, "./Analysis/summaryData.csv")
#write_csv(summaryData2, "./Analysis/summaryData2.csv")
