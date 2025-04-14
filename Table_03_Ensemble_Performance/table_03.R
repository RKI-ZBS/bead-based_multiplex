####
# Generate Table 3 - Performance parameters for ML predictions
# Daniel Stern
# RKI
# Version 1.0
# Last modified 2025-03-31

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(rio)

load("input/table_performance_establishment.Rdata")
load("input/table_performance_validatation.Rdata")

table_performance <-
  rbind(table_performance_establishment, table_performance_validatation) %>% 
  mutate(F1 = round(F1, digits = 2),
         Sensitivity = round(Sensitivity, digits = 2),
         Specificity = round(Specificity, digits = 2))
export(table_performance, "output/table_3.xlsx")
