

# Load necessary packages

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
library(viridis)

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


# Prepare data
source('./Analysis/Prepare_data.R')


# RT by FP and condition em português
rt_by_condition_port <- ggplot(data = summaryData2 %>% 
                                 group_by(participant, foreperiod, condition) %>% 
                                 summarise(meanRT = mean(meanRT)),
                               aes(x = foreperiod,
                                   y = meanRT,
                                   color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 2.7, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 2.2, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 2.0, width = 0.1, geom = "errorbar") +
  labs(title = "TR por condição e FP",
       x = "FP (s)",
       y = "TR Médio (s)",
       color = "Condição") +
  theme(plot.title = element_text(size = rel(2.0), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = rel(2.0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(1.7)),
        legend.title = element_text(size = rel(2.3)),
        legend.text = element_text(size = rel(2.1)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1)) +
  scale_color_manual(values = c("orange", "blue"), labels = c("Externa", "Ação"))

ggplot2::ggsave("./Analysis/Plots/RT_by_condition_port.tiff",
                rt_by_condition_port,
                width = 25,
                height = 16.66,
                units = "cm")

# Error rates by FP and condition em português
errors_by_condition_port <- ggplot(data = summaryDataAcc %>%
                               group_by(participant, foreperiod, condition) %>%
                               summarise(errorRate = mean(errorRate)),
                             aes(x = foreperiod,
                                 y = errorRate,
                                 color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 2.7, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 2.2, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 2.0, width = 0.1, geom = "errorbar") +
  labs(x = "FP (s)",
       y = "Proporção média de erros",
       color = "Condição",
       title = "Erros") +
  scale_color_manual(values = c("orange", "blue"),
                     label = c("Externa", "Ação")) +
  theme(plot.title = element_text(size = rel(2.0), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = rel(2.0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(1.7)),
        legend.title = element_text(size = rel(2.3)),
        legend.text = element_text(size = rel(2.1)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1),
        plot.margin = unit(c(5.5, 5.5, 5.5, 1), "pt"))
ggsave("./Analysis/Plots/error_by_condition_port.tiff",
       errors_by_condition_port,
       width = 7.7,
       height = 5.8)

# RT and errors in single panel
cond_legend <- gtable_filter(ggplot_gtable(ggplot_build(errors_by_condition_port)), "guide-box")


# Visualize
grid.arrange(rt_by_condition_port + labs(tag = "A") + 
               theme(legend.position = "none",
                     axis.text = element_text(size = rel(2.1)),
                     axis.title = element_text(size = rel(2.3)),
                     plot.title = element_text(size = rel(2.4)),
                     plot.tag = element_text(size = rel(2.5))),
             errors_by_condition_port + labs(tag = "B") +
               theme(legend.position = "none",
                     axis.text = element_text(size = rel(2.1)),
                     axis.title = element_text(size = rel(2.3)),
                     plot.title = element_text(size = rel(2.4)),
                     plot.tag = element_text(size = rel(2.5))),
             cond_legend,
             nrow = 1,
             widths = c(4/9, 4/9, 1/9))

# Save plots
comb_plots <- arrangeGrob(rt_by_condition_port + labs(tag = "A") + 
                            theme(legend.position = "none",
                                  axis.text = element_text(size = rel(2.4)),
                                  axis.title = element_text(size = rel(2.6)),
                                  plot.title = element_text(size = rel(2.7)),
                                  plot.tag = element_text(size = rel(3.4))),
                          errors_by_condition_port + labs(tag = "B") +
                            theme(legend.position = "none",
                                  axis.text = element_text(size = rel(2.4)),
                                  axis.title = element_text(size = rel(2.6)),
                                  plot.title = element_text(size = rel(2.7)),
                                  plot.tag = element_text(size = rel(3.4))),
                          cond_legend,
                          nrow = 1,
                          widths = c(4/9, 4/9, 1/9))

ggsave("./Analysis/Plots/comb_plots_portugues.jpeg",
       comb_plots,
       width = 16,
       height = 8.89)

# Sequential effects separated by condition em portugues
ggplot(data = summaryData2 %>%
         group_by(participant, foreperiod, condition, oneBackFP) %>%
         summarise(meanRT = mean(meanRT)),
       aes(x = foreperiod,
           y = meanRT,
           color = oneBackFP)) +
  geom_jitter(height = 0, width = 0.30, size = 2.7, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 2.2, aes(group = oneBackFP)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 2.0, width = 0.1, geom = "errorbar") +
  labs(title = "Efeitos sequenciais",
       x = expression(FP[n]   (s)),#x = "FP (s)",
       y = "TR Médio (s)",
       color = expression(FP[n-1])) +
  theme(plot.title = element_text(size = rel(2.7), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(2.6)),
        axis.title = element_text(size = rel(2.6)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(2.4)),
        legend.title = element_text(size = rel(2.3)),
        legend.text = element_text(size = rel(2.1)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1)) +
  facet_wrap(~condition, labeller = labeller(condition = c("external" = "Externa",
                                                           "action" = "Ação"))) + 
  scale_color_viridis(discrete = TRUE, end = 0.90)
ggsave("./Analysis/Plots/seqEff_port.png",
       width = 35,
       height = 23.32,
       units = "cm")


# Sequential effects separated by FP n-1 em portugues
ggplot(data = summaryData2 %>%
         group_by(participant, foreperiod, condition, oneBackFP) %>%
         summarise(meanRT = mean(meanRT)),
       aes(x = foreperiod,
           y = meanRT,
           color = condition)) +
  geom_jitter(height = 0, width = 0.30, size = 3.5, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 3.9, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 3.7, width = 0.1, geom = "errorbar") +
  labs(title = "Experimento 4: efeitos sequenciais",
       x = "FP (s)",
       y = "TR Médio (s)",
       color = "Condição") +
  theme(plot.title = element_text(size = rel(2.0), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(2.0)),
        axis.title = element_text(size = rel(2.0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(1.7)),
        legend.title = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.6)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1)) +
  facet_wrap(~oneBackFP,
             labeller = as_labeller(c(`1` = "FP[n-1] == 1.0",
                                      `2.8` = "FP[n-1] == 2.8"),
                                    default = label_parsed)) +
  scale_color_manual(values = c("orange", "blue"), labels = c("Externa", "Ação"))
ggsave("./Analysis/Plots/seqEff_port.tiff",
       width = 35,
       height = 23.32,
       units = "cm")