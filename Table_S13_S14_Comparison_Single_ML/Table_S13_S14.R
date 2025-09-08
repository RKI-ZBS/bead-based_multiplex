####
# Compare ML analysis with single parameter analysis
# Prepare dataframes for analysis 
# Rev 02
# Daniel Stern RKI
# Version 4.0
# Last modified 2025-08-27
# Changes to version 3.0
# - Implemented predictions based on D8 and B6 (seropositivity)
# - Ratio A35/A33 (infected) and combination from B6 and ratio A35/A33 (multi)
####

####
# Prepare workspace: clean environment
rm(list = ls(all.names = TRUE))

####
# Load necessary libraries
library(rio)
library(tidyverse)
# library(caret)
library(cvms)
library(ggpubr)
library(ggthemes)
library(yardstick)

set.seed(123)

####
# Load data
# 1) Load ML predictions for complete dataset:
# data dataInput_ML_predctions.Rdata from folder 00_09
load("input/dataInput_ML_predctions.Rdata")
dataInput_ML_predictions <- 
  dataInput_ML_predctions %>% 
  mutate(max_col = if_else(max_col == "MPXV", "Mpox", max_col),
         LDA_max_col = if_else(LDA_max_col == "MPXV", "Mpox", LDA_max_col),
         GBC_max_col = if_else(GBC_max_col == "MPXV", "Mpox", GBC_max_col),
         RF_max_col = if_else(RF_max_col == "MPXV", "Mpox", RF_max_col),
         mean_max_col = if_else(mean_max_col == "MPXV", "Mpox", mean_max_col))

rm(dataInput_ML_predctions)

# 2) Load dataset with quantified values for complete dataset:
# dataInput_Analyte_Threshold.Rdata from folder 00_09
load("input/dataInput_Analyte_Threshold.Rdata")

# 3) Load heatmap input to determine performance for combined ATI-N/E8 input
# on the validation dataset
load("input/heatmap_input.Rdata")

# 4) Load threshold data
load("input/treshold_population_Rev.RData")

# 5) Load predictions for GBC with most important antigens removed
# Data generated in folder 00_13 
load("input/dataInBenchmark_ensemble.Rdata")


cutoff_E8 <- 
  threshold_population %>% 
  filter(antigene == "E8L" & CI == "median") %>% 
  pull()

cutoff_ATI <- 
  threshold_population %>% 
  filter(antigene == "ATI-N" & CI == "median") %>% 
  pull()

cutoff_B6 <- 
  threshold_population %>% 
  filter(antigene == "B6R" & CI == "median") %>% 
  pull()

cutoff_A35A33 <- 
  threshold_population %>% 
  filter(antigene == "ratioA35A33_Infected" & CI == "median") %>% 
  pull()


# 5) Load predictions validation ML 
load("input/ensemblePrediction_Rev_02.Rdata")
validation_GBC_prediction <- 
  ensemblePrediction %>% 
  select(GBC_real = real, GBC_pred = GBC_max_col) %>% 
  mutate(GBC_real = if_else(GBC_real == "MPXV", "Mpox", GBC_real),
         GBC_pred = if_else(GBC_pred == "MPXV", "Mpox", GBC_pred))

validation_GBC_infected <- 
  ensemblePrediction %>% 
  select(GBC_real = real, GBC_pred = GBC_max_col) %>% 
  mutate(GBC_inf_real = if_else(GBC_real == "MPXV", "Mpox", "Non infected"),
         GBC_inf_pred = if_else(GBC_pred == "MPXV", "Mpox", "Non infected"))

validation_GBC_serostatus <- 
  ensemblePrediction %>% 
  select(GBC_real = real, GBC_pred = GBC_max_col) %>% 
  mutate(GBC_sero_real = if_else(GBC_real %in% c("MPXV", "MVA"), "positive", "negative"),
         GBC_sero_pred = if_else(GBC_pred %in% c("MPXV", "MVA"), "positive", "negative"))


####
# Analyse performance for prediction on two outcomes based on ATI-N or E8
# Recode GBC algorithm to mirror the outcome


# Prediction for Mpox based on single antigen
dataInput_MPXV_ATI_N <-
  dataInput_Analyte_Threshold %>% 
  filter(analyte == "ATI-N") %>% 
  filter(panel_detail %in% c("Pre", "MVA", "Mpox")) %>% 
  mutate(prediction_ATI = if_else(dataIn > Infected, "Mpox", "Non infected"),
         ground_truth_ATI = if_else(panel_detail == "Mpox", "Mpox", "Non infected")) %>% 
  unique()

dataInput_Positive_E8 <-
  dataInput_Analyte_Threshold %>% 
  filter(analyte == "E8") %>% 
  filter(panel_detail %in% c("Pre", "MVA", "Mpox")) %>% 
  mutate(prediction_pos_E8 = if_else(dataIn > Serostatus, "Positive", "Negative"),
         ground_truth_pos_E8 = if_else(panel_detail %in% c("Mpox", "MVA"), "Positive", "Negative")) %>% 
  select(sampleID_metadata, prediction_pos_E8, ground_truth_pos_E8)

dataInput_Positive_B6 <-
  dataInput_Analyte_Threshold %>% 
  filter(analyte == "B6") %>% 
  filter(panel_detail %in% c("Pre", "MVA", "Mpox")) %>% 
  mutate(prediction_pos_B6 = if_else(dataIn > Serostatus, "Positive", "Negative"),
         ground_truth_pos_B6 = if_else(panel_detail %in% c("Mpox", "MVA"), "Positive", "Negative")) %>% 
  select(sampleID_metadata, prediction_pos_B6, ground_truth_pos_B6)


# Predictions for Mpox based on ratio A35/A33
dataInput_MPXV_ratioA35A33 <-
  dataInput_Analyte_Threshold %>% 
  filter(analyte == "A35_A33") %>% 
  filter(panel_detail %in% c("Pre", "MVA", "Mpox")) %>% 
  mutate(prediction_A35A33 = if_else(dataIn > Infected, "Mpox", "Non infected"),
         ground_truth_A35A33 = if_else(panel_detail == "Mpox", "Mpox", "Non infected")) %>% 
  select(sampleID_metadata, prediction_A35A33, ground_truth_A35A33)

# Prediction for Mpox based on GBC
dataInput_MPXV_ML <-
  dataInput_ML_predictions %>% 
  mutate(prediction_GBC = if_else(GBC_max_col =="Mpox", "Mpox", "Non infected"),
         prediction_RF = if_else(RF_max_col =="Mpox", "Mpox", "Non infected"),
         prediction_LDA = if_else(LDA_max_col =="Mpox", "Mpox", "Non infected"),
         ground_truth_ML = if_else(real =="Mpox", "Mpox", "Non infected")) %>% 
  unique()

# Prediction for seropositive based on GBC
dataInput_Positive_ML <-
  dataInput_ML_predictions %>% 
  mutate(prediction_pos_GBC = if_else(GBC_max_col %in% c("Mpox", "MVA"), "Positive", "Negative"),
         prediction_pos_RF = if_else(RF_max_col %in% c("Mpox", "MVA"), "Positive", "Negative"),
         prediction_pos_LDA = if_else(LDA_max_col %in% c("Mpox", "MVA"),"Positive", "Negative"),
         ground_truth_pos_ML = if_else(real %in% c("Mpox", "MVA"), "Positive", "Negative")) %>% 
  unique()

# Combine in one dataframe
dataInput_ATI <-
  dataInput_MPXV_ATI_N %>% 
  left_join(dataInput_MPXV_ML, by = "sampleID_metadata") %>% 
  left_join(dataInput_Positive_E8, by = "sampleID_metadata") %>% 
  left_join(dataInput_Positive_ML, by = "sampleID_metadata") %>% 
  left_join(dataInput_Positive_B6, by = "sampleID_metadata") %>% 
  left_join(dataInput_MPXV_ratioA35A33, by = "sampleID_metadata") 


# Set up prediction for bootstrap analysis
# Prediction for Mpox based on ML algorithm
truth      <- dataInput_ATI$ground_truth_ATI    
truth_pos      <- dataInput_ATI$ground_truth_pos_E8 
pred_mod_ATI  <- dataInput_ATI$prediction_ATI   
pred_mod_A35A33 <- dataInput_ATI$prediction_A35A33
pred_mod_GBC  <- dataInput_ATI$prediction_GBC      
pred_mod_pos_E8 <-dataInput_ATI$prediction_pos_E8
pred_mod_pos_B6 <- dataInput_ATI$prediction_pos_B6
pred_mod_pos_GBC <- dataInput_ATI$prediction_pos_GBC

n <- length(truth)
n_boot <- 2000      

# Define variables for differences
diffs_ATI_GBC <- numeric(n_boot)
diffs_A35A33_GBC <- numeric(n_boot)
diffs_f1_pos_E8_GBC <- numeric(n_boot)
diffs_f1_pos_B6_GBC <- numeric(n_boot)

# Define variables for f1 values, precision and recall
f1_model_ATI <- numeric(n_boot)
f1_model_A35A33 <- numeric(n_boot)
f1_model_pos_E8 <- numeric(n_boot)
f1_model_pos_B6 <- numeric(n_boot)
f1_model_GBC <- numeric(n_boot)
f1_model_pos_GBC <- numeric(n_boot)

# run bootstrap analysis
for(i in 1:n_boot) {
  idx <- sample(1:n, size = n, replace = TRUE)
  # f1_ATI <- F_meas(data = factor(pred_mod_ATI[idx]), reference = factor(truth[idx]))
  f1_ATI <- yardstick::f_meas_vec(estimate = factor(pred_mod_ATI[idx],
                                                    levels = c("Non infected", "Mpox")),
                                  truth = factor(truth[idx],
                                                 levels = c("Non infected", "Mpox")), estimator = "binary",
                                  event_level = "second")
  f1_A35A33 <- yardstick::f_meas_vec(estimate = factor(pred_mod_A35A33[idx],
                                                       levels = c("Non infected", "Mpox")),
                                     truth = factor(truth[idx],
                                                    levels = c("Non infected", "Mpox")), estimator = "binary",
                                     event_level = "second")
  
  #  f1_pos_E8 <- F_meas(data = factor(pred_mod_pos_E8[idx]), reference = factor(truth_pos[idx]))
  f1_pos_E8 <- yardstick::f_meas_vec(estimate = factor(pred_mod_pos_E8[idx],
                                                       levels = c("Negative", "Positive")),
                                     truth = factor(truth_pos[idx],
                                                    levels = c("Negative", "Positive")), estimator = "binary",
                                     event_level = "second")
  f1_pos_B6 <- yardstick::f_meas_vec(estimate = factor(pred_mod_pos_B6[idx],
                                                       levels = c("Negative", "Positive")),
                                     truth = factor(truth_pos[idx],
                                                    levels = c("Negative", "Positive")), estimator = "binary",
                                     event_level = "second")
  #  f1_pos_GBC <- F_meas(data = factor(pred_mod_pos_GBC[idx]), reference = factor(truth_pos[idx]))
  f1_pos_GBC <- yardstick::f_meas_vec(estimate = factor(pred_mod_pos_GBC[idx],
                                                        levels = c("Negative", "Positive")),
                                      truth = factor(truth_pos[idx],
                                                     levels = c("Negative", "Positive")), estimator = "binary",
                                      event_level = "second")
  #  f1_GBC <- F_meas(data = factor(pred_mod_GBC[idx]), reference = factor(truth[idx]))
  f1_GBC <- yardstick::f_meas_vec(estimate = factor(pred_mod_GBC[idx],
                                                    levels = c("Non infected", "Mpox")),
                                  truth = factor(truth[idx],
                                                 levels = c("Non infected", "Mpox")), estimator = "binary",
                                  event_level = "second")
  f1_model_ATI[i] <- f1_ATI
  f1_model_A35A33[i] <- f1_A35A33
  f1_model_pos_E8[i] <- f1_pos_E8
  f1_model_pos_B6[i] <- f1_pos_B6
  f1_model_GBC[i] <- f1_GBC
  f1_model_pos_GBC[i] <- f1_pos_GBC
  diffs_ATI_GBC[i] <- f1_ATI - f1_GBC
  diffs_f1_pos_E8_GBC[i] <- f1_pos_E8 - f1_pos_GBC
}


# Predict based on both E8 and ATI-N on establishment panel
dataInput_Combined_ATI_E8 <-
  dataInput_MPXV_ATI_N %>% 
  left_join(dataInput_Positive_E8, by = "sampleID_metadata") %>% 
  mutate(prediction_combined = case_when(prediction_ATI == "Mpox" & prediction_pos_E8 == "Positive" ~ "Mpox",
                                         prediction_ATI == "Non infected" & prediction_pos_E8 == "Positive" ~ "MVA",
                                         prediction_ATI == "Non infected" & prediction_pos_E8 == "Negative" ~ "Pre",
                                         prediction_ATI == "Mpox" & prediction_pos_E8 == "Negative" ~ "Pre")) %>% 
  left_join(dataInput_ML_predictions, by = "sampleID_metadata") %>% 
  filter(!is.na(prediction_combined)) 

# Prediction outcome on validation panel based on both E8 and ATI-N 
prediction_validation_ATI_E8_combined <-
  heatmap_input_IgG %>% 
  select(panel, real, `ATI-N`, E8) %>% 
  mutate(real = if_else(real == "MPXV", "Mpox", real),
         cutoff_ATI = cutoff_ATI, cutoff_E8 = cutoff_E8,
         prediction_ATI = if_else(`ATI-N` < cutoff_ATI, "Non infected", "Mpox"),
         prediction_pos_E8 = if_else(E8 < cutoff_E8, "Negative", "Positive"),
         prediction_combined = case_when(prediction_ATI == "Mpox" & prediction_pos_E8 == "Positive" ~ "Mpox",
                                         prediction_ATI == "Non infected" & prediction_pos_E8 == "Positive" ~ "MVA",
                                         prediction_ATI == "Non infected" & prediction_pos_E8 == "Negative" ~ "Pre",
                                         prediction_ATI == "Mpox" & prediction_pos_E8 == "Negative" ~ "Pre"))

prediction_validation_ATI <-
  heatmap_input_IgG %>% 
  select(panel, real, `ATI-N`) %>% 
  mutate(real = if_else(real == "MPXV", "Mpox", "Non infected"),
         cutoff_ATI = cutoff_ATI, 
         prediction_ATI = if_else(`ATI-N` < cutoff_ATI, "Non infected", "Mpox"))

prediction_validation_E8 <-
  heatmap_input_IgG %>% 
  select(panel, real, E8) %>% 
  mutate(real_serostatus = if_else(real %in% c("MPXV", "MVA"), "positive", "negative"),
         cutoff_E8 = cutoff_E8, 
         prediction_serostatus = if_else(E8 < cutoff_E8, "negative", "positive"))



# Prediction based on both B6 and A35_A33 on establishment panel
dataInput_Combined_B6_A35A33 <-
  dataInput_MPXV_ratioA35A33 %>% 
  left_join(dataInput_Positive_B6, by = "sampleID_metadata") %>% 
  mutate(prediction_combined = case_when(prediction_A35A33 == "Mpox" & prediction_pos_B6 == "Positive" ~ "Mpox",
                                         prediction_A35A33 == "Non infected" & prediction_pos_B6 == "Positive" ~ "MVA",
                                         prediction_A35A33 == "Non infected" & prediction_pos_B6 == "Negative" ~ "Pre",
                                         prediction_A35A33 == "Mpox" & prediction_pos_B6 == "Negative" ~ "Pre")) %>% 
  left_join(dataInput_ML_predictions, by = "sampleID_metadata") %>% 
  filter(!is.na(prediction_combined)) 


# Prediction based on both B6 and A35_A33 on validation panel
prediction_validation_B6_A35A33_combined <-
  heatmap_input_IgG %>% 
  select(panel, real, B6, A35, A33) %>% 
  mutate(real = if_else(real == "MPXV", "Mpox", real),
         cutoff_A35A33 = cutoff_A35A33,
         cutoff_B6 = cutoff_B6,
         A35_A33 = A35/A33,
         prediction_A35A33 = if_else(A35_A33 < cutoff_A35A33, "Non infected", "Mpox"),
         prediction_pos_B6 = if_else(B6 < cutoff_B6, "Negative", "Positive"),
         prediction_combined = case_when(prediction_A35A33 == "Mpox" & prediction_pos_B6 == "Positive" ~ "Mpox",
                                         prediction_A35A33 == "Non infected" & prediction_pos_B6 == "Positive" ~ "MVA",
                                         prediction_A35A33 == "Non infected" & prediction_pos_B6 == "Negative" ~ "Pre",
                                         prediction_A35A33 == "Mpox" & prediction_pos_B6 == "Negative" ~ "Pre"))

prediction_validation_A35A33 <-
  heatmap_input_IgG %>% 
  select(panel, real, A35, A33) %>% 
  mutate(real = if_else(real == "MPXV", "Mpox", "Non infected"),
         cutoff_A35A33 = cutoff_A35A33,
         A35_A33 = A35/A33,
         prediction_A35A33 = if_else(A35_A33 <  cutoff_A35A33, "Non infected", "Mpox"))

prediction_validation_B6 <-
  heatmap_input_IgG %>% 
  select(panel, real, B6) %>% 
  mutate(real_serostatus = if_else(real %in% c("MPXV", "MVA"), "positive", "negative"),
         cutoff_B6 = cutoff_B6, 
         prediction_serostatus = if_else(B6 < cutoff_B6, "negative", "positive"))






####
# Boolean on validation panel
n_val <- length(prediction_validation_ATI_E8_combined$prediction_combined)
n_boot <- 2000 
pred_mod_ATI_E8_combined_val  <- prediction_validation_ATI_E8_combined$prediction_combined
truth_combined_val      <- prediction_validation_ATI_E8_combined$real
pred_mod_GBC_val <- validation_GBC_prediction$GBC_pred
truth_GBC_val <- validation_GBC_prediction$GBC_real

f1_model_ATI_E8_combined_val <- numeric(n_boot)
f1_model_B6_A35A33_combined_val <- numeric(n_boot)
f1_model_GBC_val <- numeric(n_boot)
f1_model_ATI_val <- numeric(n_boot)
f1_model_A35A33_val <- numeric(n_boot)
f1_model_E8_val <- numeric(n_boot)
f1_model_B6_val <- numeric(n_boot)
f1_model_GBC_inf_val <- numeric(n_boot)
f1_model_GBC_sero_val <- numeric(n_boot)


for(i in 1:n_boot) {
  idx <- sample(1:n_val, size = n_val, replace = TRUE)
  f1_ATI_E8_combined_val <- yardstick::f_meas_vec(estimate = factor(prediction_validation_ATI_E8_combined$prediction_combined[idx],
                                                                    levels = c("Pre", "MVA", "Mpox")),
                                                  truth = factor(prediction_validation_ATI_E8_combined$real[idx],
                                                                 levels = c("Pre", "MVA", "Mpox")), estimator = "macro")
  
  f1_B6_A35A33_combined_val <- yardstick::f_meas_vec(estimate = factor(prediction_validation_B6_A35A33_combined$prediction_combined[idx],
                                                                       levels = c("Pre", "MVA", "Mpox")),
                                                     truth = factor(prediction_validation_B6_A35A33_combined$real[idx],
                                                                    levels = c("Pre", "MVA", "Mpox")), estimator = "macro")
  
  f1_GBC_val <- yardstick::f_meas_vec(estimate = factor(validation_GBC_prediction$GBC_pred[idx],
                                                        levels = c("Pre", "MVA", "Mpox")),
                                      truth = factor(validation_GBC_prediction$GBC_real[idx],
                                                     levels = c("Pre", "MVA", "Mpox")), estimator = "macro")
  
  f1_ATI_val <- yardstick::f_meas_vec(estimate = factor(prediction_validation_ATI$prediction_ATI[idx],
                                                        levels = c("Non infected", "Mpox")),
                                      truth = factor(prediction_validation_ATI$real[idx],
                                                     levels = c("Non infected", "Mpox")), estimator = "binary", event_level = "second")
  
  
  f1_A35A33_val <- yardstick::f_meas_vec(estimate = factor(prediction_validation_A35A33$prediction_A35A33[idx],
                                                           levels = c("Non infected", "Mpox")),
                                         truth = factor(prediction_validation_A35A33$real[idx],
                                                        levels = c("Non infected", "Mpox")), estimator = "binary", event_level = "second")
  
  f1_E8_val <- yardstick::f_meas_vec(estimate = factor(prediction_validation_E8$prediction_serostatus[idx],
                                                       levels =  c("negative", "positive")),
                                     truth = factor(prediction_validation_E8$real_serostatus[idx],
                                                    levels =  c("negative", "positive")), estimator = "binary", event_level = "second")
  
  f1_B6_val <- yardstick::f_meas_vec(estimate = factor(prediction_validation_B6$prediction_serostatus[idx],
                                                       levels =  c("negative", "positive")),
                                     truth = factor(prediction_validation_B6$real_serostatus[idx],
                                                    levels =  c("negative", "positive")), estimator = "binary", event_level = "second")
  
  
  f1_GBC_inf_val <- yardstick::f_meas_vec(estimate = factor(validation_GBC_infected$GBC_inf_pred[idx],
                                                            levels = c("Non infected", "Mpox")),
                                          truth = factor(validation_GBC_infected$GBC_inf_real[idx],
                                                         levels = c("Non infected", "Mpox")), estimator = "binary",
                                          event_level = "second")
  f1_GBC_sero_val <- yardstick::f_meas_vec(estimate = factor(validation_GBC_serostatus$GBC_sero_pred[idx],
                                                             levels = c("negative", "positive")),
                                           truth = factor(validation_GBC_serostatus$GBC_sero_real[idx],
                                                          levels = c("negative", "positive")), estimator = "binary",
                                           event_level = "second")
  
  
  
  f1_model_ATI_E8_combined_val[i] <-  f1_ATI_E8_combined_val
  f1_model_B6_A35A33_combined_val[i] <-  f1_B6_A35A33_combined_val
  
  f1_model_GBC_val[i] <-  f1_GBC_val
  
  f1_model_ATI_val[i] <-    f1_ATI_val
  f1_model_A35A33_val[i] <-    f1_A35A33_val
  f1_model_GBC_inf_val[i] <-   f1_GBC_inf_val
  
  f1_model_E8_val[i] <-  f1_E8_val
  f1_model_B6_val[i] <-  f1_B6_val
  f1_model_GBC_sero_val[i] <- f1_GBC_sero_val
}

####

# Prediction for Mpox based on ML algorithm on establishment panel
truth_combined      <- dataInput_Combined_ATI_E8$panel_detail   
pred_mod_ATI_E8_combined  <- dataInput_Combined_ATI_E8$prediction_combined  
pre_mod_B6_A35A33_combined <- dataInput_Combined_B6_A35A33$prediction_combined

pred_mod_GBC_combined  <- dataInput_Combined_ATI_E8$GBC_max_col 

diffs_f1_Combined_GBC <- numeric(n_boot)
diffs_prec_Combined_GBC <- numeric(n_boot)
diffs_rec_Combined_GBC <- numeric(n_boot)

f1_model_ATI_E8_combined <- numeric(n_boot)
prec_model_ATI_E8_combined <- numeric(n_boot)
rec_model_ATI_E8_combined <- numeric(n_boot)

f1_model_B6_A35A33_combined <- numeric(n_boot)
prec_model_B6_A35A33_combined <- numeric(n_boot)
rec_model_B6_A35A33_combined <- numeric(n_boot)

f1_model_GBC_combined <- numeric(n_boot)
prec_model_GBC_combined <- numeric(n_boot)
rec_model_GBC_combined <- numeric(n_boot)


for(i in 1:n_boot) {
  idx <- sample(1:n, size = n, replace = TRUE)
  f1_ATI_E8_combined <- yardstick::f_meas_vec(estimate = factor(dataInput_Combined_ATI_E8$prediction_combined[idx], levels = c("Pre", "MVA", "Mpox")), truth = factor(dataInput_Combined_ATI_E8$panel_detail[idx], levels = c("Pre", "MVA", "Mpox")), estimator = "macro")
  f1_B6_A35A33_combined <- yardstick::f_meas_vec(estimate = factor(dataInput_Combined_B6_A35A33$prediction_combined[idx], levels = c("Pre", "MVA", "Mpox")), truth = factor(dataInput_Combined_ATI_E8$panel_detail[idx], levels = c("Pre", "MVA", "Mpox")), estimator = "macro")
  f1_Combined_GBC <- yardstick::f_meas_vec(estimate = factor(dataInput_Combined_ATI_E8$GBC_max_col[idx], levels = c("Pre", "MVA", "Mpox")), truth = factor(dataInput_Combined_ATI_E8$panel_detail[idx], levels = c("Pre", "MVA", "Mpox")), estimator = "macro")
  f1_model_ATI_E8_combined[i] <-  f1_ATI_E8_combined
  f1_model_B6_A35A33_combined[i] <-  f1_B6_A35A33_combined
  f1_model_GBC_combined[i] <-  f1_Combined_GBC
  diffs_f1_Combined_GBC[i] <- f1_ATI_E8_combined - f1_Combined_GBC
}

# Mittelwert und 95%-Konfidenzintervall
mean_diffs_f1_Combined_GBC <- mean( diffs_f1_Combined_GBC)
ci_f1_Combined_GBC <- quantile(diffs_f1_Combined_GBC, c(0.025, 0.975))

mean(f1_model_GBC_combined)
quantile(f1_model_GBC_combined, c(0.025, 0.975))

mean(f1_model_ATI_E8_combined)
quantile(f1_model_ATI_E8_combined, c(0.025, 0.975))

mean(f1_model_B6_A35A33_combined)
quantile(f1_model_B6_A35A33_combined, c(0.025, 0.975))

####
# Predictions based on GBC ensemble with baseline model and single antigens 
# removed
# 1) Recode ground truth to binary infected and binary serostatus
recode_function <- function(dataIn){
  output <- dataIn %>% 
    select(real = real, pred = GBC_max_col) %>% 
    mutate(inf_real = if_else(real == "MPXV", "Mpox", "Non infected"),
           inf_pred = if_else(pred == "MPXV", "Mpox", "Non infected"),
           sero_real = if_else(real %in% c("MPXV", "MVA"), "positive", "negative"),
           sero_pred = if_else(pred %in% c("MPXV", "MVA"), "positive", "negative"))
  return(output)
}

bootIn_0_comb <- recode_function(ensemblePrediction_0_combined)
bootIn_0_vali <- recode_function(ensemblePrediction_0_validation)

bootIn_1_comb <- recode_function(ensemblePrediction_1_combined)
bootIn_1_vali <- recode_function(ensemblePrediction_1_validation)

bootIn_2_comb <- recode_function(ensemblePrediction_2_combined)
bootIn_2_vali <- recode_function(ensemblePrediction_2_validation)

bootIn_3_comb <- recode_function(ensemblePrediction_3_combined)
bootIn_3_vali <- recode_function(ensemblePrediction_3_validation)

bootIn_4_comb <- recode_function(ensemblePrediction_4_combined)
bootIn_4_vali <- recode_function(ensemblePrediction_4_validation)

bootIn_5_comb <- recode_function(ensemblePrediction_5_combined)
bootIn_5_vali <- recode_function(ensemblePrediction_5_validation)

# 2) Define variables need for predictions across all variables
f1_model_0_comb_multi <- numeric(n_boot)
f1_model_0_comb_infected <- numeric(n_boot)
f1_model_0_comb_serostatus <- numeric(n_boot)

f1_model_0_vali_multi <- numeric(n_boot)
f1_model_0_vali_infected <- numeric(n_boot)
f1_model_0_vali_serostatus <- numeric(n_boot)

f1_model_1_comb_multi <- numeric(n_boot)
f1_model_1_comb_infected <- numeric(n_boot)
f1_model_1_comb_serostatus <- numeric(n_boot)

f1_model_1_vali_multi <- numeric(n_boot)
f1_model_1_vali_infected <- numeric(n_boot)
f1_model_1_vali_serostatus <- numeric(n_boot)

f1_model_2_comb_multi <- numeric(n_boot)
f1_model_2_comb_infected <- numeric(n_boot)
f1_model_2_comb_serostatus <- numeric(n_boot)

f1_model_2_vali_multi <- numeric(n_boot)
f1_model_2_vali_infected <- numeric(n_boot)
f1_model_2_vali_serostatus <- numeric(n_boot)

f1_model_3_comb_multi <- numeric(n_boot)
f1_model_3_comb_infected <- numeric(n_boot)
f1_model_3_comb_serostatus <- numeric(n_boot)

f1_model_3_vali_multi <- numeric(n_boot)
f1_model_3_vali_infected <- numeric(n_boot)
f1_model_3_vali_serostatus <- numeric(n_boot)

f1_model_4_comb_multi <- numeric(n_boot)
f1_model_4_comb_infected <- numeric(n_boot)
f1_model_4_comb_serostatus <- numeric(n_boot)

f1_model_4_vali_multi <- numeric(n_boot)
f1_model_4_vali_infected <- numeric(n_boot)
f1_model_4_vali_serostatus <- numeric(n_boot)

f1_model_5_comb_multi <- numeric(n_boot)
f1_model_5_comb_infected <- numeric(n_boot)
f1_model_5_comb_serostatus <- numeric(n_boot)

f1_model_5_vali_multi <- numeric(n_boot)
f1_model_5_vali_infected <- numeric(n_boot)
f1_model_5_vali_serostatus <- numeric(n_boot)

# 3) Bootstrap analysis on combined cohort
for(i in 1:n_boot) {
  idx <- sample(1:n, size = n, replace = TRUE)
  f1_0_comb_multi <- yardstick::f_meas_vec(estimate = factor(bootIn_0_comb$pred[idx], levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(bootIn_0_comb$real[idx], levels = c("Pre", "MVA", "MPXV")),
                                           estimator = "macro")
  f1_0_comb_infected <- yardstick::f_meas_vec(estimate = factor(bootIn_0_comb$inf_pred[idx], levels = c("Non infected", "Mpox")),
                                              truth = factor(bootIn_0_comb$inf_real[idx], levels = c("Non infected", "Mpox")),
                                                             estimator = "binary",
                                                             event_level = "second")
  f1_0_comb_serostatus <- yardstick::f_meas_vec(estimate = factor(bootIn_0_comb$sero_pred[idx], levels = c("negative", "positive")),
                                                truth = factor(bootIn_0_comb$sero_real[idx], levels = c("negative", "positive")),
                                                estimator = "binary",
                                                event_level = "second")
  
  f1_model_0_comb_multi[i] <-   f1_0_comb_multi
  f1_model_0_comb_infected[i] <- f1_0_comb_infected
  f1_model_0_comb_serostatus[i] <- f1_0_comb_serostatus
  
  # 1: ATI-N removed
  f1_1_comb_multi <- yardstick::f_meas_vec(estimate = factor(bootIn_1_comb$pred[idx], levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(bootIn_1_comb$real[idx], levels = c("Pre", "MVA", "MPXV")),
                                           estimator = "macro")
  f1_1_comb_infected <- yardstick::f_meas_vec(estimate = factor(bootIn_1_comb$inf_pred[idx], levels = c("Non infected", "Mpox")),
                                              truth = factor(bootIn_1_comb$inf_real[idx], levels = c("Non infected", "Mpox")),
                                              estimator = "binary",
                                              event_level = "second")
  f1_1_comb_serostatus <- yardstick::f_meas_vec(estimate = factor(bootIn_1_comb$sero_pred[idx], levels = c("negative", "positive")),
                                                truth = factor(bootIn_1_comb$sero_real[idx], levels = c("negative", "positive")),
                                                estimator = "binary",
                                                event_level = "second")
  
  f1_model_1_comb_multi[i] <-   f1_1_comb_multi
  f1_model_1_comb_infected[i] <- f1_1_comb_infected
  f1_model_1_comb_serostatus[i] <- f1_1_comb_serostatus
  
  
  # 2: ATI-N and E8 removed
  f1_2_comb_multi <- yardstick::f_meas_vec(estimate = factor(bootIn_2_comb$pred[idx], levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(bootIn_2_comb$real[idx], levels = c("Pre", "MVA", "MPXV")),
                                           estimator = "macro")
  f1_2_comb_infected <- yardstick::f_meas_vec(estimate = factor(bootIn_2_comb$inf_pred[idx], levels = c("Non infected", "Mpox")),
                                              truth = factor(bootIn_2_comb$inf_real[idx], levels = c("Non infected", "Mpox")),
                                              estimator = "binary",
                                              event_level = "second")
  f1_2_comb_serostatus <- yardstick::f_meas_vec(estimate = factor(bootIn_2_comb$sero_pred[idx], levels = c("negative", "positive")),
                                                truth = factor(bootIn_2_comb$sero_real[idx], levels = c("negative", "positive")),
                                                estimator = "binary",
                                                event_level = "second")
  
  f1_model_2_comb_multi[i] <-   f1_2_comb_multi
  f1_model_2_comb_infected[i] <- f1_2_comb_infected
  f1_model_2_comb_serostatus[i] <- f1_2_comb_serostatus
  
  
  # 3: ATI-N and E8 and D8 removed
  f1_3_comb_multi <- yardstick::f_meas_vec(estimate = factor(bootIn_3_comb$pred[idx], levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(bootIn_3_comb$real[idx], levels = c("Pre", "MVA", "MPXV")),
                                           estimator = "macro")
  f1_3_comb_infected <- yardstick::f_meas_vec(estimate = factor(bootIn_3_comb$inf_pred[idx], levels = c("Non infected", "Mpox")),
                                              truth = factor(bootIn_3_comb$inf_real[idx], levels = c("Non infected", "Mpox")),
                                              estimator = "binary",
                                              event_level = "second")
  f1_3_comb_serostatus <- yardstick::f_meas_vec(estimate = factor(bootIn_3_comb$sero_pred[idx], levels = c("negative", "positive")),
                                                truth = factor(bootIn_3_comb$sero_real[idx], levels = c("negative", "positive")),
                                                estimator = "binary",
                                                event_level = "second")
  
  f1_model_3_comb_multi[i] <-   f1_3_comb_multi
  f1_model_3_comb_infected[i] <- f1_3_comb_infected
  f1_model_3_comb_serostatus[i] <- f1_3_comb_serostatus
  
  
  # 4: D8 removed
  f1_4_comb_multi <- yardstick::f_meas_vec(estimate = factor(bootIn_4_comb$pred[idx], levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(bootIn_4_comb$real[idx], levels = c("Pre", "MVA", "MPXV")),
                                           estimator = "macro")
  f1_4_comb_infected <- yardstick::f_meas_vec(estimate = factor(bootIn_4_comb$inf_pred[idx], levels = c("Non infected", "Mpox")),
                                              truth = factor(bootIn_4_comb$inf_real[idx], levels = c("Non infected", "Mpox")),
                                              estimator = "binary",
                                              event_level = "second")
  f1_4_comb_serostatus <- yardstick::f_meas_vec(estimate = factor(bootIn_4_comb$sero_pred[idx], levels = c("negative", "positive")),
                                                truth = factor(bootIn_4_comb$sero_real[idx], levels = c("negative", "positive")),
                                                estimator = "binary",
                                                event_level = "second")
  
  f1_model_4_comb_multi[i] <-   f1_4_comb_multi
  f1_model_4_comb_infected[i] <- f1_4_comb_infected
  f1_model_4_comb_serostatus[i] <- f1_4_comb_serostatus
  
  
  
  # 5: E8 removed
  f1_5_comb_multi <- yardstick::f_meas_vec(estimate = factor(bootIn_5_comb$pred[idx], levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(bootIn_5_comb$real[idx], levels = c("Pre", "MVA", "MPXV")),
                                           estimator = "macro")
  f1_5_comb_infected <- yardstick::f_meas_vec(estimate = factor(bootIn_5_comb$inf_pred[idx], levels = c("Non infected", "Mpox")),
                                              truth = factor(bootIn_5_comb$inf_real[idx], levels = c("Non infected", "Mpox")),
                                              estimator = "binary",
                                              event_level = "second")
  f1_5_comb_serostatus <- yardstick::f_meas_vec(estimate = factor(bootIn_5_comb$sero_pred[idx], levels = c("negative", "positive")),
                                                truth = factor(bootIn_5_comb$sero_real[idx], levels = c("negative", "positive")),
                                                estimator = "binary",
                                                event_level = "second")
  
  f1_model_5_comb_multi[i] <-   f1_5_comb_multi
  f1_model_5_comb_infected[i] <- f1_5_comb_infected
  f1_model_5_comb_serostatus[i] <- f1_5_comb_serostatus
}

# 4) Bootstrap analysis on indendent validation cohort
for(i in 1:n_boot) {
  idx <- sample(1:n_val, size = n_val, replace = TRUE)
  f1_0_vali_multi <- yardstick::f_meas_vec(estimate = factor(bootIn_0_vali$pred[idx], levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(bootIn_0_vali$real[idx], levels = c("Pre", "MVA", "MPXV")),
                                           estimator = "macro")
  f1_0_vali_infected <- yardstick::f_meas_vec(estimate = factor(bootIn_0_vali$inf_pred[idx], levels = c("Non infected", "Mpox")),
                                              truth = factor(bootIn_0_vali$inf_real[idx], levels = c("Non infected", "Mpox")),
                                              estimator = "binary",
                                              event_level = "second")
  f1_0_vali_serostatus <- yardstick::f_meas_vec(estimate = factor(bootIn_0_vali$sero_pred[idx], levels = c("negative", "positive")),
                                                truth = factor(bootIn_0_vali$sero_real[idx], levels = c("negative", "positive")),
                                                estimator = "binary",
                                                event_level = "second")
  
  f1_model_0_vali_multi[i] <-   f1_0_vali_multi
  f1_model_0_vali_infected[i] <- f1_0_vali_infected
  f1_model_0_vali_serostatus[i] <- f1_0_vali_serostatus
  
  # 1: ATI-N removed
  f1_1_vali_multi <- yardstick::f_meas_vec(estimate = factor(bootIn_1_vali$pred[idx], levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(bootIn_1_vali$real[idx], levels = c("Pre", "MVA", "MPXV")),
                                           estimator = "macro")
  f1_1_vali_infected <- yardstick::f_meas_vec(estimate = factor(bootIn_1_vali$inf_pred[idx], levels = c("Non infected", "Mpox")),
                                              truth = factor(bootIn_1_vali$inf_real[idx], levels = c("Non infected", "Mpox")),
                                              estimator = "binary",
                                              event_level = "second")
  f1_1_vali_serostatus <- yardstick::f_meas_vec(estimate = factor(bootIn_1_vali$sero_pred[idx], levels = c("negative", "positive")),
                                                truth = factor(bootIn_1_vali$sero_real[idx], levels = c("negative", "positive")),
                                                estimator = "binary",
                                                event_level = "second")
  
  f1_model_1_vali_multi[i] <-   f1_1_vali_multi
  f1_model_1_vali_infected[i] <- f1_1_vali_infected
  f1_model_1_vali_serostatus[i] <- f1_1_vali_serostatus
  
  
  # 2: ATI-N and E8 removed
  f1_2_vali_multi <- yardstick::f_meas_vec(estimate = factor(bootIn_2_vali$pred[idx], levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(bootIn_2_vali$real[idx], levels = c("Pre", "MVA", "MPXV")),
                                           estimator = "macro")
  f1_2_vali_infected <- yardstick::f_meas_vec(estimate = factor(bootIn_2_vali$inf_pred[idx], levels = c("Non infected", "Mpox")),
                                              truth = factor(bootIn_2_vali$inf_real[idx], levels = c("Non infected", "Mpox")),
                                              estimator = "binary",
                                              event_level = "second")
  f1_2_vali_serostatus <- yardstick::f_meas_vec(estimate = factor(bootIn_2_vali$sero_pred[idx], levels = c("negative", "positive")),
                                                truth = factor(bootIn_2_vali$sero_real[idx], levels = c("negative", "positive")),
                                                estimator = "binary",
                                                event_level = "second")
  
  f1_model_2_vali_multi[i] <-   f1_2_vali_multi
  f1_model_2_vali_infected[i] <- f1_2_vali_infected
  f1_model_2_vali_serostatus[i] <- f1_2_vali_serostatus
  
  
  # 3: ATI-N and E8 and D8 removed
  f1_3_vali_multi <- yardstick::f_meas_vec(estimate = factor(bootIn_3_vali$pred[idx], levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(bootIn_3_vali$real[idx], levels = c("Pre", "MVA", "MPXV")),
                                           estimator = "macro")
  f1_3_vali_infected <- yardstick::f_meas_vec(estimate = factor(bootIn_3_vali$inf_pred[idx], levels = c("Non infected", "Mpox")),
                                              truth = factor(bootIn_3_vali$inf_real[idx], levels = c("Non infected", "Mpox")),
                                              estimator = "binary",
                                              event_level = "second")
  f1_3_vali_serostatus <- yardstick::f_meas_vec(estimate = factor(bootIn_3_vali$sero_pred[idx], levels = c("negative", "positive")),
                                                truth = factor(bootIn_3_vali$sero_real[idx], levels = c("negative", "positive")),
                                                estimator = "binary",
                                                event_level = "second")
  
  f1_model_3_vali_multi[i] <-   f1_3_vali_multi
  f1_model_3_vali_infected[i] <- f1_3_vali_infected
  f1_model_3_vali_serostatus[i] <- f1_3_vali_serostatus
  
  
  # 4: D8 removed
  f1_4_vali_multi <- yardstick::f_meas_vec(estimate = factor(bootIn_4_vali$pred[idx], levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(bootIn_4_vali$real[idx], levels = c("Pre", "MVA", "MPXV")),
                                           estimator = "macro")
  f1_4_vali_infected <- yardstick::f_meas_vec(estimate = factor(bootIn_4_vali$inf_pred[idx], levels = c("Non infected", "Mpox")),
                                              truth = factor(bootIn_4_vali$inf_real[idx], levels = c("Non infected", "Mpox")),
                                              estimator = "binary",
                                              event_level = "second")
  f1_4_vali_serostatus <- yardstick::f_meas_vec(estimate = factor(bootIn_4_vali$sero_pred[idx], levels = c("negative", "positive")),
                                                truth = factor(bootIn_4_vali$sero_real[idx], levels = c("negative", "positive")),
                                                estimator = "binary",
                                                event_level = "second")
  
  f1_model_4_vali_multi[i] <-   f1_4_vali_multi
  f1_model_4_vali_infected[i] <- f1_4_vali_infected
  f1_model_4_vali_serostatus[i] <- f1_4_vali_serostatus
  
  
  
  # 5: E8 removed
  f1_5_vali_multi <- yardstick::f_meas_vec(estimate = factor(bootIn_5_vali$pred[idx], levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(bootIn_5_vali$real[idx], levels = c("Pre", "MVA", "MPXV")),
                                           estimator = "macro")
  f1_5_vali_infected <- yardstick::f_meas_vec(estimate = factor(bootIn_5_vali$inf_pred[idx], levels = c("Non infected", "Mpox")),
                                              truth = factor(bootIn_5_vali$inf_real[idx], levels = c("Non infected", "Mpox")),
                                              estimator = "binary",
                                              event_level = "second")
  f1_5_vali_serostatus <- yardstick::f_meas_vec(estimate = factor(bootIn_5_vali$sero_pred[idx], levels = c("negative", "positive")),
                                                truth = factor(bootIn_5_vali$sero_real[idx], levels = c("negative", "positive")),
                                                estimator = "binary",
                                                event_level = "second")
  
  f1_model_5_vali_multi[i] <-   f1_5_vali_multi
  f1_model_5_vali_infected[i] <- f1_5_vali_infected
  f1_model_5_vali_serostatus[i] <- f1_5_vali_serostatus
}
# END Predictions based on GBC ensemble with baseline modle and single antigens
####


# Generate plots
violin_plot_diff <-
  data_frame(ATI_N = f1_model_ATI, GBC = f1_model_GBC, 
             E8_pos = f1_model_pos_E8, GBC_pos = f1_model_pos_GBC,
             ATI_E8_combined = f1_model_ATI_E8_combined, GBC_combined = f1_model_GBC_combined,
             ATI_E8_combined_val = f1_model_ATI_E8_combined_val, GBC_combined_val = f1_model_GBC_val) %>% 
  pivot_longer(cols = everything(), names_to = "model", 
               values_to = "F1") %>% 
  mutate(model = factor(model, levels = c("ATI_N", "GBC", "E8_pos", "GBC_pos",
                                          "ATI_E8_combined", "GBC_combined",
                                          "ATI_E8_combined_val", "GBC_combined_val"),
                        labels = c("ATI-N", "GBC", "E8 Serostatus", "GBC Serostatus",
                                   "ATI-N/E8 All Classes", "GBC All Classes",
                                   "ATI-N/E8 All Classes Validation", "GBC All Classes Validation"), ordered = TRUE),
         source = case_when(grepl("ATI-N", model) |grepl("E8", model)~ "Single Antigen",
                            grepl("GBC", model) ~ "ML"),
         source = factor(source, levels = c("Single Antigen",
                                            "ML"), ordered = TRUE)) %>% 
  ggplot(mapping = aes(y = F1, x = model, fill = source)) +
  geom_violin(trim = FALSE) +
  coord_flip() +
  scale_fill_manual(name = "Algorithm", values = colorblind_pal()(8)[2:8]) +
  scale_y_continuous(name = "F1 score") +
  scale_x_discrete(name = "Antigen/Model") +
  theme_minimal()

ggsave("output/SFig11_bootstrap_F1.png", width = 7, height = 3, dpi = 600)

####
# Generate table with F1-value output
# Alle Objekte aus dem Environment, die mit "f1_model" beginnen:
f1_names <- ls(pattern = "^f1_model")

# Eine Funktion, die Mittelwert und Konfidenzintervall eines Vektors ausgibt:
get_f1_stats <- function(x) {
  v <- get(x)
  c(
    mean = mean(v),
    ci_lower = quantile(v, 0.025),
    ci_upper = quantile(v, 0.975)
  )
}

# Ergebnisse in einer Tabelle zusammenfassen:
f1_stats <- sapply(f1_names, get_f1_stats)

# In ein Datenframe umwandeln und transponieren für lesbare Form:
f1_stats_df <- as.data.frame(t(f1_stats))
f1_stats_df$Model <- rownames(f1_stats_df)
rownames(f1_stats_df) <- NULL
# f1_stats_df <- f1_stats_df[, c("Model", "mean", "ci_lower", "ci_upper")]

# Ergebnis anzeigen:
print(f1_stats_df)

f1_stats_export <-
  f1_stats_df %>% 
  mutate(Model = str_remove(Model, "f1_model_")) %>% 
  filter(Model %in% c("ATI", "A35A33", "GBC", "pos_E8", "pos_B6", "pos_GBC", "ATI_E8_combined", "B6_A35A33_combined", "GBC_combined",
                      "ATI_E8_combined_val", "B6_A35A33_combined_val",  "GBC_val",
                      "E8_val", "B6_val", "GBC_sero_val", 
                      "ATI_val", "A35A33_val",  "GBC_inf_val")) %>% 
  mutate(Model = factor(Model, levels = c("ATI", "A35A33",  "GBC", "pos_E8", "pos_B6", "pos_GBC",
                                          "ATI_E8_combined", "B6_A35A33_combined", "GBC_combined",
                                          "ATI_E8_combined_val", "B6_A35A33_combined_val", "GBC_val", 
                                          "E8_val", "B6_val", "GBC_sero_val", 
                                          "ATI_val", "A35A33_val", "GBC_inf_val"),
                        labels = c("ATI-N", "A35/A33 Ratio", "GBC Infected", "E8 Serostatus", "B6 Serostatus", "GBC Serostatus",
                                   "ATI-N/E8 All Classes", "B6 A35/A33 All Classes", "GBC All Classes",
                                   "ATI-N/E8 All Classes Validation", "B6 A35/A33 All Classes Validation", "GBC All Classes Validation",
                                   "E8 Validation", "B6 Validation", "GBC Serostatus Validation",
                                   "ATI-N Validation",  "A35/A33 Validation", "GBC Infected Validation"), ordered = TRUE),
         Antigens = case_when(grepl("ATI-N/E8", Model) ~ "ATI-N/E8",
                              grepl("ATI-N", Model) ~ "ATI-N",
                              
                              grepl("B6 A35/A33", Model) ~ "B6 A35/A33",
                              grepl("A33", Model) ~ "A35/A33",
                              grepl("B6", Model) ~ "B6",
                              grepl("E8", Model) ~ "E8",
                              grepl("GBC", Model) ~ "ML"),
         Classes = case_when(grepl("ATI-N/E8", Model) ~ "Multi",
                             grepl("B6 A35/A33", Model) ~ "Multi",
                             grepl("Serostatus", Model) ~ "Binary Serostatus",
                             grepl("Infected", Model) ~ "Binary Infected",
                             grepl("ATI-N", Model) ~ "Binary Infected",
                             grepl("A35/A33", Model) ~ "Binary Infected",
                             grepl("E8", Model) ~ "Binary Serostatus",
                             grepl("B6", Model) ~ "Binary Serostatus",
                             grepl("GBC", Model) ~ "Multi",
                             TRUE ~ "Binary"),
         Method = if_else(Antigens == "ML", "ML", "Classical"), 
         Panel = if_else(grepl("Validation", Model), "Validation", "Etablishment")) %>% 
  select(Antigens, Method, Classes, Panel, Mean = mean, Lower = `ci_lower.2.5%`, Upper = `ci_upper.97.5%`) %>% 
  mutate(Mean = round(Mean, 3),
         Lower = round(Lower, 3),
         Upper = round(Upper, 3)) %>% 
  arrange(Panel, Classes, Method)

export(file = "output/f1_stats_export.xlsx", f1_stats_export)
save(file = "output/f1_stats_export.Rdata", f1_stats_export)

f1_stats_export_removal <-
  f1_stats_df %>% 
  mutate(Model = str_remove(Model, "f1_model_")) %>% 
  filter(Model %in% c("0_comb_infected", "1_comb_infected", "4_comb_infected", "5_comb_infected", "2_comb_infected", "3_comb_infected",
                      "0_comb_serostatus", "1_comb_serostatus", "4_comb_serostatus", "5_comb_serostatus", "2_comb_serostatus", "3_comb_serostatus",
                      "0_comb_multi", "1_comb_multi", "4_comb_multi", "5_comb_multi", "2_comb_multi", "3_comb_multi",
                      "0_vali_infected", "1_vali_infected", "4_vali_infected", "5_vali_infected", "2_vali_infected", "3_vali_infected",
                      "0_vali_serostatus", "1_vali_serostatus", "4_vali_serostatus", "5_vali_serostatus", "2_vali_serostatus", "3_vali_serostatus",
                      "0_vali_multi", "1_vali_multi", "4_vali_multi", "5_vali_multi", "2_vali_multi", "3_vali_multi")) %>% 
  mutate(Model = factor(Model, levels = c("0_comb_infected", "1_comb_infected", "4_comb_infected", "5_comb_infected", "2_comb_infected", "3_comb_infected",
                                          "0_comb_serostatus", "1_comb_serostatus", "4_comb_serostatus", "5_comb_serostatus", "2_comb_serostatus", "3_comb_serostatus",
                                          "0_comb_multi", "1_comb_multi", "4_comb_multi", "5_comb_multi", "2_comb_multi", "3_comb_multi",
                                          "0_vali_infected", "1_vali_infected", "4_vali_infected", "5_vali_infected", "2_vali_infected", "3_vali_infected",
                                          "0_vali_serostatus", "1_vali_serostatus", "4_vali_serostatus", "5_vali_serostatus", "2_vali_serostatus", "3_vali_serostatus",
                                          "0_vali_multi", "1_vali_multi", "4_vali_multi", "5_vali_multi", "2_vali_multi", "3_vali_multi"), ordered = TRUE),
         Antigens = case_when(grepl("^0", Model) ~ "Baseline",
                              grepl("^1", Model) ~ "ATI-N",
                              grepl("^2", Model) ~ "ATI-N E8",
                              grepl("^3", Model) ~ "ATI-N E8 D8",
                              grepl("^4", Model) ~ "D8",
                              grepl("^5", Model) ~ "E8"),
         Classes = case_when(grepl("infected", Model) ~ "Binary Infected",
                             grepl("serostatus", Model) ~ "Binary Serostatus",
                             grepl("multi", Model) ~ "Multi"),
         Method = "Classical", 
         Panel = if_else(grepl("vali", Model), "Validation", "Combined")) %>% 
  select(Antigens, Method, Classes, Panel, Mean = mean, Lower = `ci_lower.2.5%`, Upper = `ci_upper.97.5%`) %>% 
  mutate(Mean = round(Mean, 3),
         Lower = round(Lower, 3),
         Upper = round(Upper, 3)) %>% 
  arrange(Panel, Classes, Method)


# Table S14
f1_stats_export_removal$F1_summary <- sprintf("%.2f (%.2f–%.2f)", f1_stats_export_removal$Mean, 
                                              f1_stats_export_removal$Lower, 
                                              f1_stats_export_removal$Upper)
supporting_table_s14 <-
  f1_stats_export_removal %>% 
  select(Antigens, Panel, Classes, Method, F1_summary) #%>% 
# pivot_wider(names_from = Method, values_from = F1_summary)

export(supporting_table_s14, file = "output/supporting_table_s14.xlsx")
