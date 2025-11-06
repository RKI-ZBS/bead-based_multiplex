####
# Densitiy plot to test for normal distribution
# Daniel Stern
# Robert Koch Institute
# Last Modified: 2025-09-05
# Version 1.0
####

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(ggthemes)
library(ggpubr)
####
# Load data
# Load data from establishment panel
load("input/dataInputComparePanels.Rdata")

# Load data from validation panel
load("input/heatmap_input.Rdata")

####
# Generate plot
# Stratified by antigen (facet_wrap)
# Histogram for establishment
# Separate plot for IgG and IgM

plot_normalized_grouped <- 
  dataInputComparePanels %>% 
  group_by(analyte) %>% 
  mutate(dataInScaled = scale(dataIn, center = TRUE, scale = TRUE)) %>% 
  ungroup() %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = dataInScaled)) +
  geom_density() +
  #  geom_density(aes(color = panel_detail)) + 
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "skyblue") +
  labs(
    x = "Scaled Datainput",
    y = "Density"
  ) +
  facet_wrap("analyte") +
  # facet_grid(panel_detail ~ analyte) +
  theme_pubclean() +
  theme(strip.background = element_blank(),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        #  legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plot_normalized_classes <-
  dataInputComparePanels %>% 
  group_by(analyte) %>% 
  mutate(dataInScaled = scale(dataIn, center = TRUE, scale = TRUE)) %>% 
  ungroup() %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = dataInScaled)) +
  geom_density() +
  #  geom_density(aes(color = panel_detail)) + 
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "skyblue") +
  labs(
    x = "Scaled Datainput",
    y = "Density"
  ) +
  # facet_wrap("analyte") +
  facet_grid(panel_detail ~ analyte) +
  theme_pubclean() +
  theme(strip.background = element_blank(),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        #  legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("output/Fig_S12.png", plot_normalized_classes,
       width = 12, height = 5, dpi = 600)

table_shapiro_analyte <-
  dataInputComparePanels %>% 
  group_by(analyte, panel_detail) %>% 
  mutate(dataInScaled = scale(dataIn)) %>%
  summarise(
    shapiro_p = shapiro.test(dataInScaled)$p.value,
    shapiro_stat = shapiro.test(dataInScaled)$statistic,
    norm = ifelse(shapiro_p > 0.05, "yes", "no"))

