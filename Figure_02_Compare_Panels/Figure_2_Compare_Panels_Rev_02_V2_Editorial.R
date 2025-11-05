##
# Figure 1 - Compare Panels
# Daniel Stern
# Robert Koch Institute
# ZBS 3
# Version Final - Editorial Changes
# Last modified 2025-10-30
# Editorial changes:
# - Export source datafile for each panel of the plot
# - Export combined source data file for panels b and c and supp information
# combining information on all antigen
# - Calculate exact p values and export table
##
rm(list = ls(all.names = TRUE))

library(tidyverse)
library(ggthemes)
library(ggpubr)
library(GGally)
library(ggbeeswarm)
library(rio)
library(rstatix)


##
# Load data input 
load("input/dataInputComparePanels.Rdata")

##
# Plot panels distribution as used for ML-algorithms
# New plot: parallel plot with three groups
## Stratify for childhoodimmu
plotIgG <-
  dataInputComparePanels %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         H3, 
         A33, 
         A35,
         B5, 
         B6, 
         A5, 
         `ATI-C`, 
         `ATI-N`, 
         Delta) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  ggparcoord(columns = c(6:20), alphaLines = 0.3, groupColumn = "childhoodImmuAge",
             scale = "uniminmax", boxplot = F) +
  scale_color_manual(name = "Childhood Immunisation", values = colorblind_pal()(8)) +
  theme_pubclean() +
  coord_polar() +
  facet_grid(panel_plot ~ panel_detail, labeller = as_labeller(c("1" = "Pre", "2" = "MVA",
                                                                 "3" = "MPXV",
                                                                 "Pre" = "Pre",
                                                                 "MVA" = "MVA", "Mpox" = "MPXV"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

####
# Generate source data file: Fig. 2a
source_data_Fig2a <-
  dataInputComparePanels %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         H3, 
         A33, 
         A35,
         B5, 
         B6, 
         A5, 
         `ATI-C`, 
         `ATI-N`, 
         Delta) %>% 
  mutate(panel_detail = case_when(panel_detail == "1" ~ "Pre",
                                  panel_detail == "2" ~ "MVA",
                                  panel_detail == "3" ~ "MPXV",
                                  panel_detail == "Mpox"  ~ "MPXV",
                                  TRUE ~ panel_detail),
         panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  unique()


plotIgM <-
  dataInputComparePanels %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgM")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         H3, 
         A33, 
         A35,
         B5, 
         B6, 
         A5, 
         `ATI-C`, 
         `ATI-N`, 
         Delta) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  ggparcoord(columns = c(6:20), alphaLines = 0.3, groupColumn = "childhoodImmuAge",
             scale = "uniminmax", boxplot = F) +
  scale_color_manual(name = "Childhood Immunisation", values = colorblind_pal()(8)) +
  theme_pubclean() +
  coord_polar() +
  facet_grid(panel_plot ~ panel_detail, labeller = as_labeller(c("1" = "Pre", "2" = "MPXV",
                                                                 "3" = "Mpox",
                                                                 "Pre" = "Pre",
                                                                 "MVA" = "MVA", "Mpox" = "MPXV"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

####
# Generate source data file: Fig2b
source_data_Fig2b <-
  dataInputComparePanels %>% 
  filter(isotype %in% c("IgM")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         H3, 
         A33, 
         A35,
         B5, 
         B6, 
         A5, 
         `ATI-C`, 
         `ATI-N`, 
         Delta) %>% 
  mutate(panel_detail = case_when(panel_detail == "1" ~ "Pre",
                                  panel_detail == "2" ~ "MVA",
                                  panel_detail == "3" ~ "MPXV",
                                  panel_detail == "Mpox"  ~ "MPXV",
                                  TRUE ~ panel_detail),
         panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  unique()



plotStratifiedChildhood <-
  ggarrange(plotIgG, plotIgM, ncol = 2, align = "hv", labels = c("a", "b"))

#ggsave("output/plotStratified_Rev_02.png", plotStratifiedChildhood, width = 16, height = 6,
#       dpi = 600)


####
# Plot comparison panels
dataInIgGAll <-
  dataInputComparePanels %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  select(panel, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         analyte, dataIn) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  rename(status = panel_detail) %>% 
  group_by(analyte) %>% 
  mutate(dataInNorm = (dataIn-min(dataIn))/(max(dataIn)-min(dataIn))) %>% #min max normalization for each antigen
  ungroup() %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  unique()

####
# Export source data for figs 2
source_data_Fig2cd_FigS3_FigS4 <-
  dataInIgGAll %>% 
  unique()

# Visualize results as box plots
plotIgGSpox <-
  dataInIgGAll %>% 
  filter(analyte %in% c("E8", "A35", "B6", "ATI-N")) %>% 
  filter(isotype == "IgG") %>% 
  filter(panel %in% "Spox") %>% 
  ggplot( aes(x = status, y = dataInNorm)) +
  geom_beeswarm(alpha = 0.4, aes(color = status)) +
  geom_boxplot(outliers = FALSE, fill = NA) +
  scale_color_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  scale_x_discrete(name = "", labels = c("Pre", "MVA", "MPXV")) +
  scale_y_continuous(name = "Min Max norm.") +
  facet_grid(childhoodImmuAge ~ analyte) +
  theme_pubr() +
  stat_compare_means(aes(group = childhoodImmuAge), comparisons = list(c("Pre", "MVA"), c("MVA", "Mpox"),
                                                                       c("Pre", "Mpox")),method = "t.test",
                     label = "p.signif") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotIgGAcute <-
  dataInIgGAll %>% 
  filter(analyte %in% c("E8", "A35", "B6", "ATI-N")) %>% 
  filter(isotype == "IgG") %>% 
  filter(panel %in% "Acute") %>% 
  ggplot( aes(x = status, y = dataInNorm)) +
  geom_beeswarm(alpha = 0.4, aes(color = status)) +
  geom_boxplot(outliers = FALSE, fill = NA) +
  scale_color_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  scale_x_discrete(name = "", labels = c("Pre", "MVA", "MPXV")) +
  scale_y_continuous(name = "Min Max norm.") +
  facet_grid(childhoodImmuAge ~ analyte) +
  theme_pubr() +
  stat_compare_means(aes(group = childhoodImmuAge), comparisons = list(c("Pre", "MVA"), c("MVA", "Mpox"),
                                                                       c("Pre", "Mpox")),method = "t.test",
                     label = "p.signif") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotCombinedSerogroupsIgG <- 
  ggarrange(plotIgGAcute, plotIgGSpox, ncol = 2, align = "hv", common.legend = F,
            labels = c("c", "d"))

#### Supporting figure with other antigens
plotIgGSpoxSup <-
  dataInIgGAll %>% 
  filter(!(analyte %in% c("E8", "A35", "B6", "ATI-N"))) %>% 
  filter(isotype == "IgG") %>% 
  filter(panel %in% "Spox") %>% 
  ggplot( aes(x = status, y = dataInNorm)) +
  geom_beeswarm(alpha = 0.4, aes(color = status), corral.width = 0.7) +
  geom_boxplot(outliers = FALSE, fill = NA) +
  scale_color_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  scale_x_discrete(name = "", labels = c("Pre", "MVA", "MPXV")) +
  scale_y_continuous(name = "Min Max norm.") +
  facet_grid(childhoodImmuAge ~ analyte) +
  theme_pubr() +
  stat_compare_means(aes(group = childhoodImmuAge), comparisons = list(c("Pre", "MVA"), c("MVA", "Mpox"),
                                                                       c("Pre", "Mpox")),method = "t.test",
                     label = "p.signif") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotIgGAcuteSup <-
  dataInIgGAll %>% 
  filter(!(analyte %in% c("E8", "A35", "B6", "ATI-N"))) %>% 
  filter(isotype == "IgG") %>% 
  filter(panel %in% "Acute") %>% 
  ggplot( aes(x = status, y = dataInNorm)) +
  geom_beeswarm(alpha = 0.4, aes(color = status), corral.width = 0.7) +
  geom_boxplot(outliers = FALSE, fill = NA) +
  scale_color_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  scale_x_discrete(name = "", labels = c("Pre", "MVA", "MPXV")) +
  scale_y_continuous(name = "Min Max norm.") +
  facet_grid(childhoodImmuAge ~ analyte) +
  theme_pubr() +
  stat_compare_means(aes(group = childhoodImmuAge), comparisons = list(c("Pre", "MVA"), c("MVA", "Mpox"),
                                                                       c("Pre", "Mpox")),method = "t.test",
                     label = "p.signif") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotCombinedSerogroupsIgGSup <- 
  ggarrange(plotIgGAcuteSup, plotIgGSpoxSup, nrow = 2, align = "hv", common.legend = FALSE,
            labels = c("a", "b"))

ggsave(file = "output/Fig_S3_Rev_02.png", plotCombinedSerogroupsIgGSup,
       width = 14, height = 14, dpi = 600)

####
# Plot Supporting Figure with IgM results
plotIgMSpoxSup <-
  dataInIgGAll %>% 
  filter(isotype == "IgM") %>% 
  filter(panel %in% "Spox") %>% 
  ggplot( aes(x = status, y = dataInNorm)) +
  geom_beeswarm(alpha = 0.4, aes(color = status), corral.width = 0.7) +
  geom_boxplot(outliers = FALSE, fill = NA) +
  scale_color_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  scale_x_discrete(name = "", labels = c("Pre", "MVA", "MPXV")) +
  scale_y_continuous(name = "Min Max norm.") +
  facet_grid(childhoodImmuAge ~ analyte) +
  theme_pubr() +
  stat_compare_means(aes(group = childhoodImmuAge), comparisons = list(c("Pre", "MVA"), c("MVA", "Mpox"),
                                                                       c("Pre", "Mpox")),method = "t.test",
                     label = "p.signif") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotIgMAcuteSup <-
  dataInIgGAll %>% 
  filter(isotype == "IgM") %>% 
  filter(panel %in% "Acute") %>% 
  ggplot( aes(x = status, y = dataInNorm)) +
  geom_beeswarm(alpha = 0.4, aes(color = status), corral.width = 0.7) +
  geom_boxplot(outliers = FALSE, fill = NA) +
  scale_color_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  scale_x_discrete(name = "", labels = c("Pre", "MVA", "MPXV")) +
  scale_y_continuous(name = "Min Max norm.") +
  facet_grid(childhoodImmuAge ~ analyte) +
  theme_pubr() +
  stat_compare_means(aes(group = childhoodImmuAge), comparisons = list(c("Pre", "MVA"), c("MVA", "Mpox"),
                                                                       c("Pre", "Mpox")),method = "t.test",
                     label = "p.signif") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotCombinedSerogroupsIgMSup <- 
  ggarrange(plotIgMAcuteSup, plotIgMSpoxSup, nrow = 2, align = "hv", common.legend = FALSE,
            labels = c("a", "b"))


ggsave(file = "output/Fig_S4_Rev_02.png", plotCombinedSerogroupsIgMSup,
       width = 14, height = 18, dpi = 600)




####
# Plot ratios of different antigens 
plotRatioAcuteEpiYes <-
  dataInputComparePanels  %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  filter(childhoodImmuAge =="Yes") %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         A33, 
         A35,
         B5, 
         B6) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1/L1` = M1/L1,
         `A29/A27` = A29/A27,
         `B6/B5` = B6/B5, 
         `A35/A33` = A35/A33,
         `E8/D8` = E8/D8) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1/L1`, `A29/A27`, `A35/A33`, `B6/B5`,  `E8/D8`) %>% 
  pivot_longer(cols = c("M1/L1", "A29/A27", "A35/A33", "B6/B5", "E8/D8"), names_to = "antigens", values_to = "ratios") %>% 
  filter(antigens %in% c("A35/A33", "B6/B5", "E8/D8")) %>% 
  ggplot(mapping = aes(x = panel_detail, y = ratios)) +
  geom_beeswarm(alpha = 0.4, aes(color = panel_detail), corral.width = 0.7) +
  geom_boxplot(outliers = FALSE, fill = NA) +
  scale_color_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  theme_pubr() +
  facet_grid(. ~ antigens) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("Pre", "MVA"),
                                                                               c("MVA", "Mpox"),
                                                                               c("Pre", "Mpox"))) +
  scale_x_discrete(name = "", labels = c("Pre", "MVA", "MPXV")) +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotRatioAcuteEpiNo <-
  dataInputComparePanels %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  filter(childhoodImmuAge =="No") %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         A33, 
         A35,
         B5, 
         B6) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1/L1` = M1/L1,
         `A29/A27` = A29/A27,
         `B6/B5` = B6/B5, 
         `A35/A33` = A35/A33,
         `E8/D8` = E8/D8) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1/L1`, `A29/A27`, `A35/A33`, `B6/B5`,  `E8/D8`) %>% 
  pivot_longer(cols = c("M1/L1", "A29/A27", "A35/A33", "B6/B5", "E8/D8"), names_to = "antigens", values_to = "ratios") %>% 
  filter(antigens %in% c("A35/A33", "B6/B5", "E8/D8")) %>% 
  ggplot(mapping = aes(x = panel_detail, y = ratios)) +
  geom_beeswarm(alpha = 0.4, aes(color = panel_detail), corral.width = 0.7) +
  geom_boxplot(outliers = FALSE, fill = NA) +
  scale_color_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  theme_pubr() +
  facet_grid(. ~ antigens) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("Pre", "MVA"),
                                                                               c("MVA", "Mpox"),
                                                                               c("Pre", "Mpox"))) +
  scale_x_discrete(name = "", labels = c("Pre", "MVA", "MPXV")) +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")


plotRatios <-
  ggarrange(plotRatioAcuteEpiNo, plotRatioAcuteEpiYes, ncol = 2, align = "hv", labels = c("e", "f"))


# Plot supporting figure ratios
plotRatioAcuteEpiYesSupp <-
  dataInputComparePanels  %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  filter(childhoodImmuAge =="Yes") %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         A33, 
         A35,
         B5, 
         B6) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1/L1` = M1/L1,
         `A29/A27` = A29/A27,
         `B6/B5` = B6/B5, 
         `A35/A33` = A35/A33,
         `E8/D8` = E8/D8) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1/L1`, `A29/A27`, `A35/A33`, `B6/B5`,  `E8/D8`) %>% 
  pivot_longer(cols = c("M1/L1", "A29/A27", "A35/A33", "B6/B5", "E8/D8"), names_to = "antigens", values_to = "ratios") %>% 
  filter(!(antigens %in% c("A35/A33", "B6/B5", "E8/D8"))) %>% 
  ggplot(mapping = aes(x = panel_detail, y = ratios)) +
  geom_beeswarm(alpha = 0.4, aes(color = panel_detail), corral.width = 0.7) +
  geom_boxplot(outliers = FALSE, fill = NA) +
  scale_color_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  theme_pubr() +
  facet_grid(. ~ antigens) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("Pre", "MVA"),
                                                                               c("MVA", "Mpox"),
                                                                               c("Pre", "Mpox"))) +
  scale_x_discrete(name = "", labels = c("Pre", "MVA", "MPXV")) +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotRatioAcuteEpiNoSupp <-
  dataInputComparePanels %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  filter(childhoodImmuAge =="No") %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         A33, 
         A35,
         B5, 
         B6) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1/L1` = M1/L1,
         `A29/A27` = A29/A27,
         `B6/B5` = B6/B5, 
         `A35/A33` = A35/A33,
         `E8/D8` = E8/D8) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1/L1`, `A29/A27`, `A35/A33`, `B6/B5`,  `E8/D8`) %>% 
  pivot_longer(cols = c("M1/L1", "A29/A27", "A35/A33", "B6/B5", "E8/D8"), names_to = "antigens", values_to = "ratios") %>% 
  filter(!(antigens %in% c("A35/A33", "B6/B5", "E8/D8"))) %>% 
  ggplot(mapping = aes(x = panel_detail, y = ratios)) +
  geom_beeswarm(alpha = 0.4, aes(color = panel_detail), corral.width = 0.7) +
  geom_boxplot(outliers = FALSE, fill = NA) +
  scale_color_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  theme_pubr() +
  facet_grid(. ~ antigens) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("Pre", "MVA"),
                                                                               c("MVA", "Mpox"),
                                                                               c("Pre", "Mpox"))) +
  scale_x_discrete(name = "", labels = c("Pre", "MVA", "MPXV")) +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")


plotRatiosSupp <-
  ggarrange(plotRatioAcuteEpiNoSupp, plotRatioAcuteEpiYesSupp, ncol = 2, align = "hv", labels = c("a", "b"))

ggsave(file = "output/Fig_S5_Rev_02.png", plotRatiosSupp, width = 8, height = 5,
       dpi = 600)

####
# Generate source data file: plot ratios Fig2ef, FigS5
source_data_Fig2ef_FigS5 <-
  plotRatioAcuteEpiYes <-
  dataInputComparePanels  %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  # filter(childhoodImmuAge =="Yes") %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         A33, 
         A35,
         B5, 
         B6) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1/L1` = M1/L1,
         `A29/A27` = A29/A27,
         `B6/B5` = B6/B5, 
         `A35/A33` = A35/A33,
         `E8/D8` = E8/D8) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1/L1`, `A29/A27`, `A35/A33`, `B6/B5`,  `E8/D8`) %>% 
  pivot_longer(cols = c("M1/L1", "A29/A27", "A35/A33", "B6/B5", "E8/D8"), names_to = "antigens", values_to = "ratios")


plotFigLower <-
  ggarrange( plotCombinedSerogroupsIgG,
             plotRatios, nrow = 2, align = "hv",
             heights = c(2,1.4),
             common.legend = T)
plotFig2New <- 
  ggarrange(plotStratifiedChildhood, plotFigLower,
            nrow = 2, align = "hv",
            common.legend = F, heights = c(1,1.7))

# Save figure 2
ggsave("output/Figure_2_Rev_02_V2.png", plotFig2New, width = 14, height = 12,
       dpi = 600)
ggsave("output/Figure_2_Rev_02_V2.pdf", plotFig2New, width = 14, height = 12,
       dpi = 600)


####
# Generate table with multiplex comparisons providing exact p-values and 
# test statistics for multiple comparisons presented in Fig2, S3, and S4
table_p_values_Fig2cd_FigS3_FigS4 <-
  source_data_Fig2cd_FigS3_FigS4 %>% 
  group_by(isotype, analyte, panel, childhoodImmuAge) %>% 
  select(status, dataInNorm) %>% 
  rstatix::pairwise_t_test(dataInNorm ~ status, comparisons = list(c("Pre", "MVA"),
                                                                   c("MVA", "Mpox"),
                                                                   c("Pre", "Mpox")))

table_p_values_Fig2ef_FigS5 <-
  source_data_Fig2ef_FigS5 %>% 
  group_by(isotype, panel_plot, childhoodImmuAge, antigens) %>% 
  select(panel_detail, ratios) %>% 
  rstatix::pairwise_t_test(ratios ~ panel_detail, comparisons = list(c("Pre", "MVA"),
                                                                   c("MVA", "Mpox"),
                                                                   c("Pre", "Mpox")))



####
# Export source data files
export(file = "output/Source_Data_Fig2.xlsx", list(Fig2a = source_data_Fig2a,
                                                   Fig2b = source_data_Fig2b,
                                                   Fig2cd_FigS3_FigS4 = source_data_Fig2cd_FigS3_FigS4,
                                                   Fig2ef_FigS5 = source_data_Fig2ef_FigS5,
                                                   Fig2cd_FigS3_FigS4_p_values = table_p_values_Fig2cd_FigS3_FigS4,
                                                   Fig2ef_FigS5_p_values = table_p_values_Fig2ef_FigS5))
