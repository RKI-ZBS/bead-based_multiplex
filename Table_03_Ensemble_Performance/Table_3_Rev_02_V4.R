####
# Generate Fig 6 and sup fig s11 - Performance parameters for ML predictions
# Daniel Stern
# RKI
# Version 4.0
# Last modified 2025-08-24

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(rio)
library(ggthemes)
library(ggpubr)

load("input/table_performance_establishment_Rev_04.Rdata")
load("input/table_performance_validation_Rev_02.Rdata")
load("input/f1_stats_export.Rdata")
load("input/table_GBC_class_establishment_Rev_04.Rdata")
load("input/table_GBC_class_validation_Rev_02.Rdata")

# Generate on dataframe with data from both panels
table_performance_establishment_combined <- 
  table_performance_establishment %>% 
  mutate(Cohort = "Combined",
         Algorithm = case_when(grepl("GBC", Model) ~ "GBC",
                               grepl("LDA", Model) ~ "LDA",
                               grepl("RF", Model) ~ "RF",
                               grepl("mean", Model) ~ "Ensemble mean",
                               grepl("serostatus", Model) ~ "Ensemble serostatus"),
         filtered = case_when(grepl("filtered", Model, ignore.case = TRUE) ~ "confidence > 0.5",
                              TRUE ~ "none"),
         n = case_when(Algorithm == "Ensemble serostatus" & filtered == "confidence > 0.5" ~ 1155,
                       Algorithm == "Ensemble mean" & filtered == "confidence > 0.5" ~ 1129,
                       Algorithm == "LDA" & filtered == "confidence > 0.5" ~ 1140,
                       Algorithm == "GBC" & filtered == "confidence > 0.5" ~ 1156,
                       Algorithm == "RF" & filtered == "confidence > 0.5" ~ 1173,
                       filtered == "none" ~ 1260),
         Class = case_when(grepl("MPXV", Model) ~ "MPXV",
                           grepl("MVA", Model) ~ "MVA",
                           grepl("Pre", Model) ~ "Pre",
                           TRUE ~ "Macro Avg."))


table_performance_validation_combined <- 
  table_performance_validation %>% 
  mutate(Cohort = "Independent Validation",
         Algorithm = case_when(grepl("GBC", Model) ~ "GBC",
                               grepl("LDA", Model) ~ "LDA",
                               grepl("RF", Model) ~ "RF",
                               grepl("mean", Model) ~ "Ensemble mean",
                               grepl("serostatus", Model) ~ "Ensemble serostatus"),
         filtered = case_when(grepl("filtered", Model, ignore.case = TRUE) ~ "confidence > 0.5",
                              TRUE ~ "none"),
         n = case_when(filtered == "confidence > 0.5" ~ 124,
                       filtered == "none" ~ 143),
         Class = case_when(grepl("MPXV", Model) ~ "MPXV",
                           grepl("MVA", Model) ~ "MVA",
                           grepl("Pre", Model) ~ "Pre",
                           TRUE ~ "Macro Avg.")) 

table_performance <-
  rbind(table_performance_establishment_combined, table_performance_validation_combined) %>% 
  mutate(Algorithm = factor(Algorithm, levels = c("LDA", "RF", "GBC",
                                                  "Ensemble mean", "Ensemble serostatus"),
                            ordered = TRUE),
         filtered = factor(filtered, levels = c("none", "confidence > 0.5"),
                           ordered = TRUE)) %>% 
  arrange(filtered, Cohort, Algorithm) %>% 
  select(Cohort, Algorithm, Class, filtered, n, F1 = F1_formatted, Precision = precision_formatted, 
         Recall = recall_formatted, Sensitivity = sens_formatted, Specificity = spec_formatted)

# Filter new table 3 with only confident prediction
table_3 <-
  table_performance %>% 
  filter(Algorithm != "Ensemble mean") %>% 
  filter(Class == "Macro Avg.") #%>% 
 # arrange(Algorithm)

# Filter new supporting table 11 with class wise predictions on independent validation
table_s11 <-
  table_performance %>% 
  filter(Class != "Macro Avg.") %>% 
  filter(Algorithm != "Ensemble mean") %>% 
  arrange(Algorithm)


# Export new table 3 with only confident predictions
export(table_3, "output/table_3.xlsx")

# Export new supporting table iwth all predictions
export(table_s11, "output/table_s11.xlsx")

