

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


# Remove extreme values
data <- data %>%
  filter(RT < 1.0) %>%
  filter(RT > 0.15)

# Transform RT to reduce skew
data$logRT <- log10(data$RT) # log-transform
data$invRT <- 1/data$RT


# Trimming
data2 <- data %>%
  group_by(participant) %>%
  mutate(RTzscore = compute_zscore(RT),
         logRTzscore = compute_zscore(logRT)) %>%
  filter(abs(logRTzscore) < 3) %>%
  ungroup()

# No trimming
data <- data %>%
  group_by(participant) %>%
  mutate(RTzscore=ifelse(!is.na(RT), compute_zscore(RT), NA),
         logRTzscore=ifelse(!is.na(RT), compute_zscore(logRT), NA)) %>%
  #filter(abs(logRTzscore) < 3) %>%
  ungroup()


# Average data
summaryData <- data %>%
  group_by(participant, FP, condition, FPType, conFPDur,
           orientation, prevOri, seqOri,
           oneBackFP, block, countBalance, condTestPos) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT)) %>%
  ungroup()

summaryData2 <- data2 %>%
  group_by(participant, FP, condition, FPType, conFPDur,
           orientation, prevOri, seqOri,
           oneBackFP, block, countBalance, condTestPos) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT)) %>%
  ungroup()


summaryDataAll <- dataAll %>%
  group_by(participant,FP,condition,FPType,conFPDur,
           oneBackFP) %>%
  summarise(meanRT = mean(RT),
            meanAcc = mean(Acc)) %>%
  ungroup()

#write_csv(data, "./Analysis/data.csv")
#write_csv(data2, "./Analysis/data2.csv")
#write_csv(summaryData, "./Analysis/summaryData.csv")
#write_csv(summaryData2, "./Analysis/summaryData2.csv")
