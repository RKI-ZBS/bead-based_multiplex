####
# Analyse Reproducibility
# Daniel Stern
# RKI
# 2025-08-01
# Version 1.0
####

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(caret)
library(ggpubr)
library(ggthemes)
library(broom)

# Load data
load("input/dataInRep.Rdata")

# Load first experiment for repeatability
load("input/dataInRepSpoxFiltered.Rdata")

# Harmonize and prepare dataframes for joining
dataInValidation <-
  dataInRep %>% 
  select(sampleID_metadata, experiment, panel, analyte, isotype, dataIn, 
         serostatus, serostatus_cat) %>% 
  filter(experiment %in% c("Rev_01_HSA", "Rev_01_Rep_HSA")) %>% 
  mutate(experiment = if_else(grepl("Rep", experiment), "1", "2"),
         panel = "Validation")

dataInSpox <- 
  dataInRepSpoxFiltered %>% 
  select(sampleID_metadata, experiment = panel, panel = panel_detail, analyte, isotype, 
         dataIn, serostatus, serostatus_cat) %>% 
  mutate(experiment = if_else(grepl("Rep", experiment), "1", "2"),
         panel = "Epi")

# Combine both dataframe for a combined analysis of all repliced measurement
dataInCombined <-
  rbind(dataInValidation, dataInSpox) %>% 
  filter(sampleID_metadata != "N_27287ad4") %>% # Remove duplicated measurement
  filter(sampleID_metadata != "S_5c13d8fe") %>% # Remove outlier for epi measurement
  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                              "D8L", "E8", "H3L",
                                              "L1R", "M1",
                                              "A33R", "A35R",
                                              "B5R", "B6",
                                              "A5L", "ATI-C", "ATI-N", "VACV", "Delta"),
                          labels = c("A27", "A29",
                                     "D8", "E8", "H3",
                                     "L1", "M1",
                                     "A33", "A35",
                                     "B5", "B6",
                                     "A5", "ATI-C", "ATI-N", "VACV", "Delta"), ordered = TRUE))

# Generate wide dataframe
dataInCombinedWide <-
  dataInCombined %>% 
  pivot_wider(names_from = experiment, values_from = c(dataIn, serostatus, serostatus_cat)) %>% 
  filter(!is.na(dataIn_1)) %>% 
  filter(!is.na(dataIn_2))

# Identify outlier
# dataInCombinedWide %>% 
#  filter(dataIn_1 > 2 & dataIn_2 < 2 & analyte == "VACV")


# Plot Repsoducibilty Scatter Plot
plot_reproducibilityIgG <-
  dataInCombinedWide %>% 
  filter(isotype == "IgG") %>% 
  filter(analyte != "VACV") %>% 
  ggplot(mapping = aes(x =dataIn_1, y = dataIn_2, color = panel)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +  # Gerade mit Steigung 1 hinzufügen
  facet_grid(panel ~ analyte) +
  theme_minimal() +
  scale_color_manual(name = "Sero Cohort", values = colorblind_pal()(8)[2:8]) +
  scale_x_continuous(name = "Multiplex quant.",  breaks = c(0,2,4,6), limits = c(0,6)) +
  scale_y_continuous(name = "Multiplex quant. rep.",  breaks = c(0,2,4,6), limits = c(0,6)) +
  theme(strip.background = element_blank(),
        legend.position = "top",
        #  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plot_reproducibilityIgM <-
  dataInCombinedWide %>% 
  filter(isotype == "IgM") %>% 
  filter(analyte != "VACV") %>% 
  ggplot(mapping = aes(x =dataIn_1, y = dataIn_2, color = panel)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +  # Gerade mit Steigung 1 hinzufügen
  facet_grid(panel ~ analyte) +
  theme_minimal() +
  scale_color_manual(name = "Sero Cohort", values = colorblind_pal()(8)[2:8]) +
  scale_x_continuous(name = "Multiplex quant.", breaks = c(0,2,4,6), limits = c(0,6)) +
  scale_y_continuous(name = "Multiplex quant. rep.",  breaks = c(0,2,4,6), limits = c(0,6)) +
  theme(strip.background = element_blank(),
        legend.position = "top",
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

####
# Calculate linear regression
linear_regression <-
  dataInCombinedWide%>% 
  group_by(panel, analyte, isotype) %>%
  do(broom::glance(lm(dataIn_1 ~ dataIn_2, data = .)))

plot_linear_regression_IgG_Epi <-
  linear_regression %>% 
  filter(isotype == "IgG") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel == "Epi") %>% 
  group_by(panel, isotype) %>% 
  mutate(analyte = factor(analyte, levels = analyte[order(r.squared, decreasing = TRUE)])) %>% 
  ggplot(mapping = aes(x = analyte, y = r.squared, color = panel)) +
  geom_point(size = 3) +
 # facet_grid(. ~ isotype) +
  scale_color_manual(name = "Sero Cohort", values = colorblind_pal()(8)[2:8]) +
  scale_y_continuous(name = "R Square", limits = c(0,1)) + 
  scale_x_discrete(name = "") +
  theme_minimal() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plot_linear_regression_IgG_Validation <-
  linear_regression %>% 
  filter(isotype == "IgG") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel == "Validation") %>% 
  group_by(panel, isotype) %>% 
  mutate(analyte = factor(analyte, levels = analyte[order(r.squared, decreasing = TRUE)])) %>% 
  ggplot(mapping = aes(x = analyte, y = r.squared, color = panel)) +
  geom_point(size = 3) +
 # facet_grid(. ~ isotype) +
  scale_color_manual(name = "Sero Cohort", values = colorblind_pal()(8)[3:8]) +
  scale_y_continuous(name = "R Square", limits = c(0,1)) + 
  scale_x_discrete(name = "") +
  theme_minimal() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plot_linear_regression_IgM_Epi <-
  linear_regression %>% 
  filter(isotype == "IgM") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel == "Epi") %>% 
  group_by(panel, isotype) %>% 
  mutate(analyte = factor(analyte, levels = analyte[order(r.squared, decreasing = TRUE)])) %>% 
  ggplot(mapping = aes(x = analyte, y = r.squared, color = panel)) +
  geom_point(size = 3) +
  #facet_grid(. ~ isotype) +
  scale_color_manual(name = "Sero Cohort", values = colorblind_pal()(8)[2:8]) +
  scale_y_continuous(name = "R Square", limits = c(0,1)) + 
  scale_x_discrete(name = "") +
  theme_minimal() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plot_linear_regression_IgM_Validation <-
  linear_regression %>% 
  filter(isotype == "IgM") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel == "Validation") %>% 
  group_by(panel, isotype) %>% 
  mutate(analyte = factor(analyte, levels = analyte[order(r.squared, decreasing = TRUE)])) %>% 
  ggplot(mapping = aes(x = analyte, y = r.squared, color = panel)) +
  geom_point(size = 3) +
 # facet_grid(. ~ isotype) +
  scale_color_manual(name = "Sero Cohort", values = colorblind_pal()(8)[3:8]) +
  scale_y_continuous(name = "R Square", limits = c(0,1)) + 
  scale_x_discrete(name = "") +
 theme_minimal() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



####
# Test agreement between borderline measurement
dataInput_agreement_Epi_Clear <- 
  dataInCombinedWide %>% 
  filter(panel == "Epi") %>% 
  filter(!(grepl("borderline", serostatus_cat_1) | (grepl("borderline", serostatus_cat_2))))

epi_clear_conf <-
  confusionMatrix(data = factor(dataInput_agreement_Epi_Clear$serostatus_2),
                  reference = factor(dataInput_agreement_Epi_Clear$serostatus_1),
                  positive = "1", mode = "everything")

dataInput_agreement_Validation_Clear <- 
  dataInCombinedWide %>% 
  filter(panel == "Validation") %>% 
  filter(!(grepl("borderline", serostatus_cat_1) | (grepl("borderline", serostatus_cat_2))))

validation_clear_conf <-
  confusionMatrix(data = factor(dataInput_agreement_Validation_Clear$serostatus_2),
                  reference = factor(dataInput_agreement_Validation_Clear$serostatus_1),
                  positive = "1", mode = "everything")


dataInput_agreement_Epi_All <- 
  dataInCombinedWide %>% 
  filter(panel == "Epi")

epi_all_conf <-
  confusionMatrix(data = factor(dataInput_agreement_Epi_All$serostatus_2),
                  reference = factor(dataInput_agreement_Epi_All$serostatus_1),
                  positive = "1", mode = "everything")

dataInput_agreement_Validation_All <- 
  dataInCombinedWide %>% 
  filter(panel == "Validation") 

validataion_all_conf <-
  confusionMatrix(data = factor(dataInput_agreement_Validation_All$serostatus_2),
                  reference = factor(dataInput_agreement_Validation_All$serostatus_1),
                  positive = "1", mode = "everything")


confusion_function <- function(analyteIn, isotypeIn, panelIn){
  dataInput_agreement <- 
    dataInCombinedWide %>% 
    filter(analyte == analyteIn) %>% 
    filter(isotype == isotypeIn) %>% 
    filter(panel == panelIn)
  
  output <- 
    confusionMatrix(data = factor(dataInput_agreement$serostatus_2),
                    reference = factor(dataInput_agreement$serostatus_1),
                    positive = "1", mode = "everything")
  return(output)
}

conf_A27_IgG <-
  confusion_function("A27", "IgG", "Validation")
conf_A27_IgM <-
  confusion_function("A27", "IgM", "Validation")
conf_A29_IgG <-
  confusion_function("A29", "IgG", "Validation")
conf_A29_IgM <-
  confusion_function("A29", "IgM", "Validation")

conf_L1_IgG <-
  confusion_function("L1", "IgG", "Validation")
conf_L1_IgM <-
  confusion_function("L1", "IgM", "Validation")
conf_M1_IgG <-
  confusion_function("M1", "IgG", "Validation")
conf_M1_IgM <-
  confusion_function("M1", "IgM", "Validation")

conf_D8_IgG <-
  confusion_function("D8", "IgG", "Validation")
conf_D8_IgM <-
  confusion_function("D8", "IgM", "Validation")
conf_E8_IgG <-
  confusion_function("E8", "IgG", "Validation")
conf_E8_IgM <-
  confusion_function("E8", "IgM", "Validation")


conf_A33_IgG <-
  confusion_function("A33", "IgG", "Validation")
conf_A33_IgM <-
  confusion_function("A33", "IgM", "Validation")
conf_A35_IgG <-
  confusion_function("A35", "IgG", "Validation")
conf_A35_IgM <-
  confusion_function("A35", "IgM", "Validation")


conf_B5_IgG <-
  confusion_function("B5", "IgG", "Validation")
conf_B5_IgM <-
  confusion_function("B5", "IgM", "Validation")
conf_B6_IgG <-
  confusion_function("B6", "IgG", "Validation")
conf_B6_IgM <-
  confusion_function("B6", "IgM", "Validation")


conf_H3_IgG <-
  confusion_function("H3", "IgG", "Validation")
conf_H3_IgM <-
  confusion_function("H3", "IgM", "Validation")
conf_A5_IgG <-
  confusion_function("A5", "IgG", "Validation")
conf_A5_IgM <-
  confusion_function("A5", "IgM", "Validation")


conf_ATI.N_IgG <-
  confusion_function("ATI-N", "IgG", "Validation")
conf_ATI.N_IgM <-
  confusion_function("ATI-N", "IgM", "Validation")
conf_ATI.C_IgG <-
  confusion_function("ATI-C", "IgG", "Validation")
conf_ATI.C_IgM <-
  confusion_function("ATI-C", "IgM", "Validation")

conf_Delta_IgG <-
  confusion_function("Delta", "IgG", "Validation")
conf_Delta_IgM <-
  confusion_function("Delta", "IgM", "Validation")



####
# Generate table with F1-value output
names_conf <- ls(pattern = "^conf_")


# Calculate mean and confidence intervalls
get_stats <- function(dataIn, dataInName) {
  valuesIn <- (dataIn[[4]][7])
  namesIn <- names(dataIn[[4]][7])
  
  output <- tibble(namesIn, valuesIn) %>% 
    pivot_wider(names_from = namesIn, values_from = valuesIn) %>% 
    mutate(input = str_remove(dataInName, "conf_"))
  return(output)
}


# Write data into one table
perf_matrices <- list()
for(i in 1:length(names_conf)){
  output <- get_stats(get(names_conf[i]), names_conf[i])
  perf_matrices[[i]] <- output
}

perf_dataframe <-
  do.call(rbind, perf_matrices) %>% 
  separate_wider_delim(input, delim = "_", names = c("analyte", "isotype")) %>% 
  mutate(analyte = str_replace(analyte, "\\.", "-"),
         panel = "Validation")

confusion_function_spox <- function(analyteIn, isotypeIn){
  dataInput_agreement <- 
    dataInRepSpoxFiltered%>% 
    filter(!is.na(panel)) %>% 
    filter(!is.na(dataIn)) %>% 
    select(-panel_detail) %>% 
    pivot_wider(names_from = panel, values_from = c(dataIn, serostatus_cat, serostatus)) %>% 
    mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                                "D8L", "E8", "H3L",
                                                "L1R", "M1",
                                                "A33R", "A35R",
                                                "B5R", "B6",
                                                "A5L", "ATI-C", "ATI-N", "VACV", "Delta"),
                            labels = c("A27", "A29",
                                       "D8", "E8", "H3",
                                       "L1", "M1",
                                       "A33", "A35",
                                       "B5", "B6",
                                       "A5", "ATI-C", "ATI-N", "VACV", "Delta"), ordered = TRUE)) %>% 
    filter(analyte == analyteIn) %>% 
    filter(isotype == isotypeIn)
  
  output <- 
    confusionMatrix(data = factor(dataInput_agreement$serostatus_SPox),
                    reference = factor(dataInput_agreement$serostatus_SPox_Rep),
                    positive = "1", mode = "everything")
  return(output)
}


spox_conf_A27_IgG <-
  confusion_function_spox("A27", "IgG")
spox_conf_A27_IgM <-
  confusion_function_spox("A27", "IgM")
spox_conf_A29_IgG <-
  confusion_function_spox("A29", "IgG")
spox_conf_A29_IgM <-
  confusion_function_spox("A29", "IgM")

spox_conf_L1_IgG <-
  confusion_function_spox("L1", "IgG")
spox_conf_L1_IgM <-
  confusion_function_spox("L1", "IgM")
spox_conf_M1_IgG <-
  confusion_function_spox("M1", "IgG")
spox_conf_M1_IgM <-
  confusion_function_spox("M1", "IgM")

spox_conf_D8_IgG <-
  confusion_function_spox("D8", "IgG")
spox_conf_D8_IgM <-
  confusion_function_spox("D8", "IgM")
spox_conf_E8_IgG <-
  confusion_function_spox("E8", "IgG")
spox_conf_E8_IgM <-
  confusion_function_spox("E8", "IgM")


spox_conf_A33_IgG <-
  confusion_function_spox("A33", "IgG")
spox_conf_A33_IgM <-
  confusion_function_spox("A33", "IgM")
spox_conf_A35_IgG <-
  confusion_function_spox("A35", "IgG")
spox_conf_A35_IgM <-
  confusion_function_spox("A35", "IgM")


spox_conf_B5_IgG <-
  confusion_function_spox("B5", "IgG")
spox_conf_B5_IgM <-
  confusion_function_spox("B5", "IgM")
spox_conf_B6_IgG <-
  confusion_function_spox("B6", "IgG")
spox_conf_B6_IgM <-
  confusion_function_spox("B6", "IgM")


spox_conf_H3_IgG <-
  confusion_function_spox("H3", "IgG")
spox_conf_H3_IgM <-
  confusion_function_spox("H3", "IgM")
spox_conf_A5_IgG <-
  confusion_function_spox("A5", "IgG")
spox_conf_A5_IgM <-
  confusion_function_spox("A5", "IgM")


spox_conf_ATI.N_IgG <-
  confusion_function_spox("ATI-N", "IgG")
spox_conf_ATI.N_IgM <-
  confusion_function_spox("ATI-N", "IgM")
spox_conf_ATI.C_IgG <-
  confusion_function_spox("ATI-C", "IgG")
spox_conf_ATI.C_IgM <-
  confusion_function_spox("ATI-C", "IgM")

spox_conf_Delta_IgG <-
  confusion_function_spox("Delta", "IgG")
spox_conf_Delta_IgM <-
  confusion_function_spox("Delta", "IgM")



####
# Generate table with F1-value output
names_spox <- ls(pattern = "^spox_")
#precision_names <- ls(pattern = "^precision_model")
#recall_names <- ls(pattern = "^recall_model")

# Calculate mean and confidence intervalls
get_stats_spox <- function(dataIn, dataInName) {
  valuesIn <- (dataIn[[4]][7])
  namesIn <- names(dataIn[[4]][7])
  # var_name <- deparse(substitute(dataIn))
  
  output <- tibble(namesIn, valuesIn) %>% 
    pivot_wider(names_from = namesIn, values_from = valuesIn) %>% 
    mutate(input = str_remove(dataInName, "spox_conf_"))
  return(output)
}

# Write data into one table
perf_matrices_spox <- list()
for(i in 1:length(names_spox)){
  output <- get_stats_spox(get(names_spox[i]), names_spox[i])
  perf_matrices_spox[[i]] <- output
}

perf_dataframe_spox <-
  do.call(rbind, perf_matrices_spox) %>% 
  separate_wider_delim(input, delim = "_", names = c("analyte", "isotype")) %>% 
  mutate(analyte = str_replace(analyte, "\\.", "-"),
         panel = "Epi")




plot_f1_epi_IgG <-
  perf_dataframe %>% 
  rbind(perf_dataframe_spox) %>% 
  filter(panel == "Epi") %>% 
  filter(isotype == "IgG") %>% 
  filter(analyte != "VACV") %>% 
  mutate(analyte = factor(analyte, levels = analyte[order(F1, decreasing = TRUE)])) %>% 
  ggplot(mapping = aes(x = analyte, y = F1, color = panel)) +
  geom_point(size = 3) +
 # facet_grid(. ~isotype, scales = "free") +
  scale_color_manual(name = "Sero Cohort", values = colorblind_pal()(8)[2:8]) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = "F1",
       x = "") +
  theme_minimal() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plot_f1_validation_IgG <-
  perf_dataframe %>% 
  rbind(perf_dataframe_spox) %>% 
  filter(panel == "Validation") %>% 
  filter(isotype == "IgG") %>% 
  filter(analyte != "VACV") %>% 
  mutate(analyte = factor(analyte, levels = analyte[order(F1, decreasing = TRUE)])) %>% 
  ggplot(mapping = aes(x = analyte, y = F1, color = panel)) +
  geom_point(size = 3) +
 # facet_grid(. ~isotype, scales = "free") +
  scale_color_manual(name = "Sero Cohort", values = colorblind_pal()(8)[3:8]) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = "F1",
       x = "") +
  theme_minimal() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plot_f1_epi_IgM <-
  perf_dataframe %>% 
  rbind(perf_dataframe_spox) %>% 
  filter(panel == "Epi") %>% 
  filter(isotype == "IgM") %>% 
  filter(analyte != "VACV") %>% 
  mutate(analyte = factor(analyte, levels = analyte[order(F1, decreasing = TRUE)])) %>% 
  ggplot(mapping = aes(x = analyte, y = F1, color = panel)) +
  geom_point(size = 3) +
 # facet_grid(. ~isotype, scales = "free") +
  scale_color_manual(name = "Sero Cohort", values = colorblind_pal()(8)[2:8]) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = "F1",
       x = "") +
  theme_minimal() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plot_f1_validation_IgM <-
  perf_dataframe %>% 
  rbind(perf_dataframe_spox) %>% 
  filter(panel == "Validation") %>% 
  filter(isotype == "IgM") %>% 
  filter(analyte != "VACV") %>% 
  mutate(analyte = factor(analyte, levels = analyte[order(F1, decreasing = TRUE)])) %>% 
  ggplot(mapping = aes(x = analyte, y = F1, color = panel)) +
  geom_point(size = 3) +
#  facet_grid(. ~isotype, scales = "free") +
  scale_color_manual(name = "Sero Cohort", values = colorblind_pal()(8)[3:8]) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = "F1",
       x = "") +
  theme_minimal() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))






plot_combined <-
  ggarrange(ggarrange(plot_reproducibilityIgG, plot_reproducibilityIgM,
                      ncol = 2, align = "hv", common.legend = TRUE, labels = "auto"),
            ggarrange(plot_linear_regression_IgG_Epi, plot_linear_regression_IgM_Epi,
                      plot_linear_regression_IgG_Validation, plot_linear_regression_IgM_Validation,
                      plot_f1_epi_IgG, plot_f1_epi_IgM,
                      plot_f1_validation_IgG, plot_f1_validation_IgM,
                      ncol = 2,
                      nrow = 4, 
                      align = "hv", common.legend = F, labels = c("c", "d", 
                                                                     "", "",
                                                                     "g", "h", "", "")),
            nrow = 2, common.legend = TRUE, heights = c(1,2))

ggsave("output/SFig_reproducibility.png", plot_combined,
       width = 13, height = 10, dpi = 600)

length(unique(dataInCombined$sampleID_metadata)) # 168 total # 129
length(unique(dataInRepSpoxFiltered$sampleID_metadata)) # 40 Spox
length(unique(dataInValidation$sampleID_metadata)) # 130 Validation

dataInCombined %>% 
  group_by(panel) %>% 
  select(sampleID_metadata) %>% 
  unique() %>% 
  count()
