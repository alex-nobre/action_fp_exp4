

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
  mutate(across(c(participant, Counterbalance, Handedness, block, condition, orientation, foreperiod, oneBackFP), as_factor))

data$condition <- data$condition %>%
  fct_relevel(c("external", "action"))

# Create column for previous orientation and for comparison of current and previous orientations
data <- data %>%
  mutate(seqOri = ifelse(lag(orientation)==orientation, 'same', 'different'),
         prevOri = lag(orientation)) %>%
  mutate(seqOri = as.factor(seqOri),
         prevOri = as.factor(prevOri))

# Create column for trial number
data <- data %>%
  group_by(participant) %>%
  mutate(trial = seq(1,n())) %>%
  ungroup()

# Remove trials without n-1 FP values (i.e., first of each block)
data <- data %>%
  filter(!is.na(oneBackFP))

# Save data with error trials to assess accuracy
dataAcc <- data

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

###############################################################
# Add delay data
###############################################################
# delayData <- read_csv("./Analysis/delayDataAll.csv") %>%
#   mutate(across(c(participant, condition), as_factor)) %>%
#   select(-condition)
# 
# data <- inner_join(data, delayData, by = c("trial", "participant"))
# data2 <- inner_join(data2, delayData, by = c("trial", "participant"))
# dataAcc <- inner_join(dataAcc, delayData, by = c("trial", "participant"))

################################################################

# Average data
summaryData <- data %>%
  group_by(participant, foreperiod, condition,
           orientation, prevOri, seqOri,
           oneBackFP, block, Counterbalance, Handedness) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT)) %>% #,
            #meanDelay = mean(delay)) %>%
  ungroup()

summaryData2 <- data2 %>%
  group_by(participant, foreperiod, condition,
           orientation, prevOri, seqOri,
           oneBackFP, block, Counterbalance, Handedness) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT)) %>% #,
            #meanDelay = mean(delay)) %>%
  ungroup()


summaryDataAll <- dataAcc %>%
  group_by(participant, foreperiod, condition, oneBackFP, Counterbalance, Handedness) %>%
  summarise(meanRT = mean(RT),
            meanAcc = mean(Acc)) %>% #,
            #meanDelay = mean(delay)) %>%
  ungroup()

#write_csv(data, "./Analysis/data.csv")
#write_csv(data2, "./Analysis/data2.csv")
#write_csv(summaryData, "./Analysis/summaryData.csv")
#write_csv(summaryData2, "./Analysis/summaryData2.csv")
