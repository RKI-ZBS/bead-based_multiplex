####
# Analyse Benchmark Study
# Daniel Stern
# Robert Koch Institute
# Version 6.0
# Last modified: 2025-08-25
####

rm(list = ls(all.names = TRUE))

library(rio)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggbeeswarm)
library(stringr)

load(file = "input/dataInputFeatElim.Rdata")




baseline <- dataInput %>% 
  filter(panel == "0_baseline") %>% 
  filter(epi_cohort == "all_all")


table_performance <- dataInput_feat_select %>% 
  filter(epi_cohort == "all_all") %>% 
  mutate(establishment_cohort = case_when(establishment_cohort == "xgboost" ~ "Establishment",
                                          establishment_cohort == "xgboost_revised_data" ~ "Validation",
                                          TRUE ~ establishment_cohort)) %>% 
  filter(panel != "0_baseline") %>% 
  rbind(baseline) %>% 
  mutate(isotypes = str_replace(isotypes, "_", " "),
         panel = factor(panel, levels = c("0_baseline",
                                          "1_removed", 
                                          "2_removed",
                                          "3_removed",
                                          "4_removed",
                                          "5_removed",
                                          "6_removed", 
                                          "7_removed",
                                          "8_removed",
                                          "9_removed",
                                          "10_removed",
                                          "11_removed",
                                          "12_removed"),
                        labels = c("Baseline",
                                   "-B5-VACV", 
                                   "-A29-MPXV",
                                   "-A27-VACV", 
                                   "-H3-VACV",
                                   "-A5L-CPXV",
                                   "-ATI-C-CPXV",
                                   "-A33-VACV",
                                   "-Delta-VACV",
                                   "-B6-MPXV",
                                   "-A35-MPXV",
                                   "-ATI-N-CPXV",
                                   "-D8-VACV"),
                        ordered = TRUE)) %>% 
  group_by(panel, isotypes, establishment_cohort) %>% 
  summarise(mean_f1 = mean(f1),
            sd_f1 = sd(f1))
  



####
# Calculate pairwise decrease between different steps
# Then calculate mean difference with sd
# Export table and include in supporting information
#
output_diff <-
  dataInput_feat_select %>% 
  filter(epi_cohort == "all_all") %>% 
  filter(panel != "0_baseline") %>% 
  rbind(baseline) %>% 
  mutate(establishment_cohort = case_when(establishment_cohort == "xgboost" ~ "Combined",
                                          establishment_cohort == "xgboost_revised_data" ~ "Independent Validation",
                                          establishment_cohort == "Establishment" ~ "Combined",
                                          establishment_cohort == "Validation" ~ "Independent Validation",
                                          TRUE ~ establishment_cohort)) %>% 
  mutate(isotypes = str_replace(isotypes, "_", " "),
         panel = factor(panel, levels = c("0_baseline",
                                          "1_removed", 
                                          "2_removed",
                                          "3_removed",
                                          "4_removed",
                                          "5_removed",
                                          "6_removed", 
                                          "7_removed",
                                          "8_removed",
                                          "9_removed",
                                          "10_removed",
                                          "11_removed",
                                          "12_removed"),
                        # labels = c("Baseline",
                        #             "-B5-VACV", 
                        #            "-A29-MPXV",
                        #             "-A27-VACV", 
                        #             "-H3-VACV",
                        #             "-A5L-CPXV",
                        #             "-ATI-C-CPXV",
                        #             "-A33-VACV",
                        #             "-Delta-VACV",
                        #             "-B6-MPXV",
                        #             "-A35-MPXV",
                        #             "-ATI-N-CPXV",
                        #             "-D8-VACV"),
                        ordered = TRUE)) %>% 
  select(f1, panel, establishment_cohort, isotypes, repetition) %>% 
  pivot_wider(names_from = panel, values_from = f1) %>% 
  mutate(diff_1 = `0_baseline` - `1_removed`,
         diff_2 = `0_baseline` - `2_removed`,
         diff_3 = `0_baseline` - `3_removed`,
         diff_4 = `0_baseline` - `4_removed`,
         diff_5 = `0_baseline` - `5_removed`,
         diff_6 = `0_baseline` - `6_removed`,
         diff_7 = `0_baseline` - `7_removed`,
         diff_8 = `0_baseline` - `8_removed`,
         diff_9 = `0_baseline` - `9_removed`,
         diff_10 = `0_baseline` - `10_removed`,
         diff_11 = `0_baseline` - `11_removed`,
         diff_12 = `0_baseline` - `12_removed`) %>% 
  group_by(establishment_cohort, isotypes) %>% 
  summarise(mean_baseline = mean(`0_baseline`),
            sd_baseline = sd(`0_baseline`),
            mean_diff_1 = mean(`diff_1`),
            sd_diff_1 = sd(`diff_1`),
            mean_diff_2 = mean(`diff_2`),
            sd_diff_2 = sd(`diff_2`),
            mean_diff_3 = mean(`diff_3`),
            sd_diff_3 = sd(`diff_3`),
            mean_diff_4 = mean(`diff_4`),
            sd_diff_4 = sd(`diff_4`),
            mean_diff_5 = mean(`diff_5`),
            sd_diff_5 = sd(`diff_5`),
            mean_diff_6 = mean(`diff_6`),
            sd_diff_6 = sd(`diff_6`),
            mean_diff_7 = mean(`diff_7`),
            sd_diff_7 = sd(`diff_7`),
            mean_diff_8 = mean(`diff_8`),
            sd_diff_8 = sd(`diff_8`),
            mean_diff_9 = mean(`diff_9`),
            sd_diff_9 = sd(`diff_9`),
            mean_diff_10 = mean(`diff_10`),
            sd_diff_10 = sd(`diff_10`),
            mean_diff_11 = mean(`diff_11`),
            sd_diff_11 = sd(`diff_11`),
            mean_diff_12 = mean(`diff_12`),
            sd_diff_12 = sd(`diff_12`))



plot_diff_minpanel <-
  output_diff %>% 
  filter(establishment_cohort == "Combined") %>% 
  pivot_longer(cols = starts_with("mean"), values_to = "mean_difference", names_prefix = "mean_",
               names_to = "model_mean") %>% 
  pivot_longer(cols = starts_with("sd"), values_to = "sd_difference", names_prefix = "sd_",
               names_to = "model_sd") %>% 
  filter(model_mean == model_sd) %>% 
  filter(model_mean != "baseline") %>% 
  mutate(model_mean = as.numeric(str_remove(model_mean, "diff_"))) %>% 
  ggplot(mapping = aes(x = model_mean, y = -mean_difference, color = isotypes)) +
  geom_smooth(se = F) +
  geom_pointrange(aes(ymin = -mean_difference - sd_difference, ymax = -mean_difference + sd_difference)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -0.01, color = colorblind_pal()(8)[4]) +
  geom_hline(yintercept = -0.05, color = colorblind_pal()(8)[6]) +
  geom_vline(xintercept = 5.5,  color = colorblind_pal()(8)[4]) +
  geom_vline(xintercept = 7.5, color = colorblind_pal()(8)[6]) +
  theme_bw() +
  scale_color_manual(name = "Isotypes", values = colorblind_pal()(8)[2:8]) +
  scale_x_continuous(name = "Antigens removed", limits = c(1,12), breaks = c(1:12),
                     labels = c(
                       "-B5-VACV", 
                       "-A29-MPXV",
                       "-A27-VACV", 
                       "-H3-VACV",
                       "-A5L-CPXV",
                       "-ATI-C-CPXV",
                       "-A33-VACV",
                       "-Delta-VACV",
                       "-B6-MPXV",
                       "-A35-MPXV",
                       "-ATI-N-CPXV",
                       "-D8-VACV"))+
  scale_y_continuous(name = "F1 score decrease", breaks = c(0, -0.05, -0.10,
                                                            -0.15, -0.20, -0.25,
                                                            -0.30, -0.35), limits = c(-0.35, 0.07))+
  # facet_grid(. ~ establishment_cohort) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("output/Fig_S9_benchmark.png", plot_diff_minpanel, width = 12, height = 5, dpi = 300)

#### 
# Calculate difference for benchmark study
output_diff_benchmark <-
  dataInput %>% 
  filter(epi_cohort == "all_all") %>% 
  mutate(establishment_cohort = case_when(establishment_cohort == "Establishment" ~ "Combined",
                                          establishment_cohort == "Validation" ~ "Independent Validation",
                                          TRUE ~ establishment_cohort)) %>% 
  # filter(panel != "0_baseline") %>% 
  # rbind(baseline) %>% 
  mutate(isotypes = str_replace(isotypes, "_", " "),
         panel = factor(panel, levels = c("0_baseline",
                                          "1_ATI", 
                                          "2_D8",
                                          "3_E8",
                                          "4_ATI_E8",
                                          "5_ATI_E8_D8"),
                        labels = c("Baseline",
                                   "-ATI-N-CPXV",
                                   "-D8-VACV",
                                   "-E8-MPXV",
                                   "-ATI-N-CPXV/-E8-MPXV",
                                   "-ATI-N-CPXV/-E8-MPXV/-D8-VACV"),
                        ordered = TRUE)) %>% 
  select(f1, panel, establishment_cohort, isotypes, repetition) %>% 
  pivot_wider(names_from = panel, values_from = f1) %>% 
  mutate(diff_1 = `Baseline` - `-ATI-N-CPXV`,
         diff_2 = `Baseline` - `-D8-VACV`,
         diff_3 = `Baseline` - `-E8-MPXV`,
         diff_4 = `Baseline` - `-ATI-N-CPXV/-E8-MPXV`,
         diff_5 = `Baseline` - `-ATI-N-CPXV/-E8-MPXV/-D8-VACV`) %>% 
  group_by(establishment_cohort, isotypes) %>% 
  summarise(mean_baseline = mean(`Baseline`),
            sd_baseline = sd(`Baseline`),
            mean_diff_1 = mean(`diff_1`),
            sd_diff_1 = sd(`diff_1`),
            mean_diff_2 = mean(`diff_2`),
            sd_diff_2 = sd(`diff_2`),
            mean_diff_3 = mean(`diff_3`),
            sd_diff_3 = sd(`diff_3`),
            mean_diff_4 = mean(`diff_4`),
            sd_diff_4 = sd(`diff_4`),
            mean_diff_5 = mean(`diff_5`),
            sd_diff_5 = sd(`diff_5`))

