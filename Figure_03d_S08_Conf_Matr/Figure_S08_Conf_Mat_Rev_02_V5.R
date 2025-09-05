####
# Analyse stratified ensemble learning on test dataset
# Daniel Stern
# 2025-07-09
# Version 4.0
# Changes to version 1.0
# - Include predictions for RF
# Changes to version 2.0
# - Include calculation of mean F1 values and 95% CI using bootstrap analysis
# Changes to version 3.0
# - Include calculation of mean sensitivites and specificities using bootstrap analysis
# Changes to version 4.0
# - Harmonize filtering with regards to the ensemble prediction
####

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(ggpubr)
library(ggthemes)
library(gplots)
library(caret)
library(cvms)
library(ggimage)
library(rsvg)
library(yardstick)

set.seed(123) # Set seed for bootstrap analysis

# Load datainput
load("input/ensembleCombined_Rev_02.Rdata")

# Filter LDA und ensemble > 0.5
ensembleCombinedFiltered <-
  ensembleCombined %>% 
  filter(mean_mean_conf > 0.5) %>% 
  filter(pred_ensemble > 0.5)

####
# Calculate confusion matrix ensemble
confusionMatrixEnsemble <- 
  confusionMatrix(factor(ensembleCombined$real), factor(ensembleCombined$max_col))
confusionMatrixEnsemble[["byClass"]][ , "F1"]

confusionMatrixEnsembleFiltered <- 
  confusionMatrix(factor(ensembleCombinedFiltered$real), factor(ensembleCombinedFiltered$max_col))
confusionMatrixEnsembleFiltered[["byClass"]][ , "F1"]

misClassPre <-
  ensembleCombined %>% 
  filter(real == "Pre" & max_col == "MPXV")

misClassPreFiltered <-
  ensembleCombinedFiltered %>% 
  filter(real == "Pre" & max_col == "MPXV")

####
# Calculate confusion matrix mean of LDA and GBC unstratified
confusionMatrixMean <- 
  confusionMatrix(factor(ensembleCombined$real), factor(ensembleCombined$mean_max_col))
confusionMatrixMean[["byClass"]][ , "F1"]

# Filter LDA und ensemble > 0.5 for mean
ensembleCombinedFilteredMean <-
  ensembleCombined %>% 
  filter(mean_mean_conf > 0.5) %>% 
  filter(mean_pred > 0.5)


confusionMatrixMeanFiltered <- 
  confusionMatrix(factor(ensembleCombinedFilteredMean$real), factor(ensembleCombinedFilteredMean$mean_max_col))


misClassPre <-
  ensembleCombined %>% 
  filter(real == "Pre" & max_col == "MPXV")

misClassPreFiltered <-
  ensembleCombinedFiltered %>% 
  filter(real == "Pre" & max_col == "MPXV")


####
# Filter GBC
ensembleCombinedFilteredGBC <-
  ensembleCombined %>% 
  filter(GBC_pred_ensemble > 0.5) %>% 
  filter(mean_mean_conf > 0.5)

confusionMatrixGBC <- 
  confusionMatrix(factor(ensembleCombined$real), factor(ensembleCombined$GBC_max_col))
confusionMatrixGBC[["byClass"]][ , "F1"]

confusionMatrixGBCFiltered <- 
  confusionMatrix(factor(ensembleCombinedFilteredGBC$real), factor(ensembleCombinedFilteredGBC$GBC_max_col))
confusionMatrixGBCFiltered[["byClass"]][ , "F1"]

####
# Filter LDA
ensembleCombinedFilteredLDA <-
  ensembleCombined %>% 
  filter(mean_mean_conf > 0.5) %>% 
  filter(LDA_pred_ensemble > 0.5)

confusionMatrixLDA <- 
  confusionMatrix(factor(ensembleCombined$real), factor(ensembleCombined$LDA_max_col))
confusionMatrixLDA[["byClass"]][ , "F1"]

confusionMatrixLDAFiltered <- 
  confusionMatrix(factor(ensembleCombinedFilteredLDA$real), factor(ensembleCombinedFilteredLDA$LDA_max_col))
confusionMatrixLDAFiltered[["byClass"]][ , "F1"]


## Rev 02: Include results of RF predictions
# Filter RF
ensembleCombinedFilteredRF <-
  ensembleCombined %>% 
   filter(mean_mean_conf > 0.5) %>% 
  filter(RF_pred_ensemble > 0.5) 

confusionMatrixRF <- 
  confusionMatrix(factor(ensembleCombined$real), factor(ensembleCombined$RF_max_col))
confusionMatrixRF[["byClass"]][ , "F1"]

confusionMatrixRFFiltered <- 
  confusionMatrix(factor(ensembleCombinedFilteredRF$real), factor(ensembleCombinedFilteredRF$RF_max_col))
confusionMatrixRFFiltered[["byClass"]][ , "F1"]


ensembleCombined <-
  ensembleCombined %>% 
  filter(!is.na(max_col))

save(ensembleCombined, file = "output/ensembleCombined.Rdata")

####
# PLot confusion matrices using package cvsm on filtered data
conf_mat <- confusion_matrix(targets = factor(ensembleCombinedFiltered$real),
                             predictions = factor(ensembleCombinedFiltered$max_col))

conf_mat_LDA <- confusion_matrix(targets = factor(ensembleCombinedFilteredLDA$real),
                                 predictions = factor(ensembleCombinedFilteredLDA$LDA_max_col))

conf_mat_GBC <- confusion_matrix(targets = factor(ensembleCombinedFilteredGBC$real),
                                 predictions = factor(ensembleCombinedFilteredGBC$GBC_max_col))

conf_mat_mean <- confusion_matrix(targets = factor(ensembleCombinedFilteredMean$real),
                                  predictions = factor(ensembleCombinedFilteredMean$mean_max_col))

conf_mat_RF <- confusion_matrix(targets = factor(ensembleCombinedFilteredRF$real),
                                predictions = factor(ensembleCombinedFilteredRF$RF_max_col))

####
# Plot confusion matrices using package cvsm on unfiltered data
conf_mat_all <- confusion_matrix(targets = factor(ensembleCombined$real),
                                 predictions = factor(ensembleCombined$max_col))

conf_mat_LDA_all <- confusion_matrix(targets = factor(ensembleCombined$real),
                                     predictions = factor(ensembleCombined$LDA_max_col))

conf_mat_GBC_all <- confusion_matrix(targets = factor(ensembleCombined$real),
                                     predictions = factor(ensembleCombined$GBC_max_col))

conf_mat_mean_all <- confusion_matrix(targets = factor(ensembleCombined$real),
                                      predictions = factor(ensembleCombined$mean_max_col))

conf_mat_RF_all <- confusion_matrix(targets = factor(ensembleCombined$real),
                                    predictions = factor(ensembleCombined$RF_max_col))

sum(conf_mat_all$Table[[1]])
sum(conf_mat$Table[[1]])
sum(conf_mat_LDA$Table[[1]])
sum(conf_mat_GBC$Table[[1]])
sum(conf_mat_mean$Table[[1]])
sum(conf_mat_RF$Table[[1]])


####
# Bootstrap analysis to calculate mean and 95C% CI for the validation data
# Datainput and different cases
# Ground truth filtered: $real in all cases, use different dataframes as input
# 1) Ensemble serostatus filtered: ensembleCombinedFiltered$max_col
# 2) Ensemble mean filtered: ensembleCombinedFilteredMean$mean_max_col
# 3) LDA filtered: ensembleCombinedFilteredLDA$LDA_max_col
# 4) GBC filtered: ensembleCombinedFilteredGBC$GBC_max_col
# 5) RF filtered: ensembleCombinedFilteredRF$RF_max_col
# 6) Ensemble serostatus: ensembleCombined$max_col
# 7) Ensemble mean: ensembleCombined$mean_max_col
# 8) LDA: ensembleCombined$LDA_max_col
# 9) GBC: ensembleCombined$GBC_max_col
# 10) RF: ensembleCombined$RF_max_col

# Boolean on validation panel

# Define number of draws
n_val_serostatus_filtered <- length(ensembleCombinedFiltered$max_col) #1: 1155
n_val_mean_filtered <- length(ensembleCombinedFilteredMean$mean_max_col) #2: 1129
n_val_LDA_filtered <- length(ensembleCombinedFilteredLDA$LDA_max_col) #3: 1140
n_val_GBC_filtered <- length(ensembleCombinedFilteredGBC$GBC_max_col) #4: 1206 1156
n_val_RF_filtered <- length(ensembleCombinedFilteredRF$RF_max_col) #5: 1224 1173
n_val <-length(ensembleCombined$max_col) #6: 1260

# Define number of iterations
n_boot <- 2000 

# Define output of model as numeric vectors
f1_model_serostatus_filtered <- numeric(n_boot)
f1_model_mean_filtered <- numeric(n_boot)
f1_model_LDA_filtered <- numeric(n_boot)
f1_model_GBC_filtered <- numeric(n_boot)
f1_model_RF_filtered <- numeric(n_boot)

f1_model_serostatus <- numeric(n_boot)
f1_model_mean <- numeric(n_boot)
f1_model_LDA <- numeric(n_boot)
f1_model_GBC <- numeric(n_boot)
f1_model_RF <- numeric(n_boot)

# Define output models for recall
recall_model_serostatus_filtered <- numeric(n_boot)
recall_model_mean_filtered <- numeric(n_boot)
recall_model_LDA_filtered <- numeric(n_boot)
recall_model_GBC_filtered <- numeric(n_boot)
recall_model_RF_filtered <- numeric(n_boot)

recall_model_serostatus <- numeric(n_boot)
recall_model_mean <- numeric(n_boot)
recall_model_LDA <- numeric(n_boot)
recall_model_GBC <- numeric(n_boot)
recall_model_RF <- numeric(n_boot)

# Define output models for precision
precision_model_serostatus_filtered <- numeric(n_boot)
precision_model_mean_filtered <- numeric(n_boot)
precision_model_LDA_filtered <- numeric(n_boot)
precision_model_GBC_filtered <- numeric(n_boot)
precision_model_RF_filtered <- numeric(n_boot)

precision_model_serostatus <- numeric(n_boot)
precision_model_mean <- numeric(n_boot)
precision_model_LDA <- numeric(n_boot)
precision_model_GBC <- numeric(n_boot)
precision_model_RF <- numeric(n_boot)

# Define output models for sensitivity
sens_model_serostatus_filtered <- numeric(n_boot)
sens_model_mean_filtered <- numeric(n_boot)
sens_model_LDA_filtered <- numeric(n_boot)
sens_model_GBC_filtered <- numeric(n_boot)
sens_model_RF_filtered <- numeric(n_boot)

sens_model_serostatus <- numeric(n_boot)
sens_model_mean <- numeric(n_boot)
sens_model_LDA <- numeric(n_boot)
sens_model_GBC <- numeric(n_boot)
sens_model_RF <- numeric(n_boot)

# Define output models for specificity
spec_model_serostatus_filtered <- numeric(n_boot)
spec_model_mean_filtered <- numeric(n_boot)
spec_model_LDA_filtered <- numeric(n_boot)
spec_model_GBC_filtered <- numeric(n_boot)
spec_model_RF_filtered <- numeric(n_boot)

spec_model_serostatus <- numeric(n_boot)
spec_model_mean <- numeric(n_boot)
spec_model_LDA <- numeric(n_boot)
spec_model_GBC <- numeric(n_boot)
spec_model_RF <- numeric(n_boot)


# 1: Bootstrap for serostatus filtered
for(i in 1:n_boot) {
  idx <- sample(1:n_val_serostatus_filtered, size = n_val_serostatus_filtered, replace = TRUE)
  f1_serostatus_filtered <- yardstick::f_meas_vec(estimate = factor(ensembleCombinedFiltered$max_col[idx],
                                                                    levels = c("Pre", "MVA", "MPXV")),
                                                  truth = factor(ensembleCombinedFiltered$real[idx],
                                                                 levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  recall_serostatus_filtered <- yardstick::recall_vec(estimate = factor(ensembleCombinedFiltered$max_col[idx],
                                                                        levels = c("Pre", "MVA", "MPXV")),
                                                      truth = factor(ensembleCombinedFiltered$real[idx],
                                                                     levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  precision_serostatus_filtered <- yardstick::precision_vec(estimate = factor(ensembleCombinedFiltered$max_col[idx],
                                                                              levels = c("Pre", "MVA", "MPXV")),
                                                            truth = factor(ensembleCombinedFiltered$real[idx],
                                                                           levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  sens_serostatus_filtered <- yardstick::sens_vec(estimate = factor(ensembleCombinedFiltered$max_col[idx],
                                                                    levels = c("Pre", "MVA", "MPXV")),
                                                  truth = factor(ensembleCombinedFiltered$real[idx],
                                                                 levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  spec_serostatus_filtered <- yardstick::spec_vec(estimate = factor(ensembleCombinedFiltered$max_col[idx],
                                                                    levels = c("Pre", "MVA", "MPXV")),
                                                  truth = factor(ensembleCombinedFiltered$real[idx],
                                                                 levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  
  
  f1_model_serostatus_filtered[i] <-  f1_serostatus_filtered
  recall_model_serostatus_filtered[i] <- recall_serostatus_filtered
  precision_model_serostatus_filtered[i] <- precision_serostatus_filtered
  sens_model_serostatus_filtered[i] <- sens_serostatus_filtered
  spec_model_serostatus_filtered[i] <- spec_serostatus_filtered
}

# 2: Bootstrap for serostatus mean filtered
for(i in 1:n_boot) {
  idx <- sample(1:n_val_mean_filtered, size = n_val_mean_filtered, replace = TRUE)
  f1_mean_filtered <- yardstick::f_meas_vec(estimate = factor(ensembleCombinedFilteredMean$mean_max_col[idx],
                                                              levels = c("Pre", "MVA", "MPXV")),
                                            truth = factor(ensembleCombinedFilteredMean$real[idx],
                                                           levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  recall_mean_filtered <- yardstick::recall_vec(estimate = factor(ensembleCombinedFilteredMean$mean_max_col[idx],
                                                                  levels = c("Pre", "MVA", "MPXV")),
                                                truth = factor(ensembleCombinedFilteredMean$real[idx],
                                                               levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  precision_mean_filtered <- yardstick::precision_vec(estimate = factor(ensembleCombinedFilteredMean$mean_max_col[idx],
                                                                        levels = c("Pre", "MVA", "MPXV")),
                                                      truth = factor(ensembleCombinedFilteredMean$real[idx],
                                                                     levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  sens_mean_filtered <- yardstick::sens_vec(estimate = factor(ensembleCombinedFilteredMean$mean_max_col[idx],
                                                              levels = c("Pre", "MVA", "MPXV")),
                                            truth = factor(ensembleCombinedFilteredMean$real[idx],
                                                           levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  spec_mean_filtered <- yardstick::spec_vec(estimate = factor(ensembleCombinedFilteredMean$mean_max_col[idx],
                                                              levels = c("Pre", "MVA", "MPXV")),
                                            truth = factor(ensembleCombinedFilteredMean$real[idx],
                                                           levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  f1_model_mean_filtered[i] <-  f1_mean_filtered
  recall_model_mean_filtered[i] <- recall_mean_filtered
  precision_model_mean_filtered[i] <- precision_mean_filtered
  sens_model_mean_filtered[i] <- sens_mean_filtered
  spec_model_mean_filtered[i] <- spec_mean_filtered
}

# 3: Bootstrap for LDA filtered
for(i in 1:n_boot) {
  idx <- sample(1:n_val_LDA_filtered, size = n_val_LDA_filtered, replace = TRUE)
  f1_LDA_filtered <- yardstick::f_meas_vec(estimate = factor(ensembleCombinedFilteredLDA$LDA_max_col[idx],
                                                             levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(ensembleCombinedFilteredLDA$real[idx],
                                                          levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  recall_LDA_filtered <- yardstick::recall_vec(estimate = factor(ensembleCombinedFilteredLDA$LDA_max_col[idx],
                                                                 levels = c("Pre", "MVA", "MPXV")),
                                               truth = factor(ensembleCombinedFilteredLDA$real[idx],
                                                              levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  precision_LDA_filtered <- yardstick::precision_vec(estimate = factor(ensembleCombinedFilteredLDA$LDA_max_col[idx],
                                                                       levels = c("Pre", "MVA", "MPXV")),
                                                     truth = factor(ensembleCombinedFilteredLDA$real[idx],
                                                                    levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  sens_LDA_filtered <- yardstick::sens_vec(estimate = factor(ensembleCombinedFilteredLDA$LDA_max_col[idx],
                                                             levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(ensembleCombinedFilteredLDA$real[idx],
                                                          levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  spec_LDA_filtered <- yardstick::spec_vec(estimate = factor(ensembleCombinedFilteredLDA$LDA_max_col[idx],
                                                             levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(ensembleCombinedFilteredLDA$real[idx],
                                                          levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  f1_model_LDA_filtered[i] <-  f1_LDA_filtered
  recall_model_LDA_filtered[i] <- recall_LDA_filtered
  precision_model_LDA_filtered[i] <- precision_LDA_filtered
  sens_model_LDA_filtered[i] <- sens_LDA_filtered
  spec_model_LDA_filtered[i] <- spec_LDA_filtered
}

# 4: Bootstrap for GBC filtered
for(i in 1:n_boot) {
  idx <- sample(1:n_val_GBC_filtered, size = n_val_GBC_filtered, replace = TRUE)
  f1_GBC_filtered <- yardstick::f_meas_vec(estimate = factor(ensembleCombinedFilteredGBC$GBC_max_col[idx],
                                                             levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(ensembleCombinedFilteredGBC$real[idx],
                                                          levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  recall_GBC_filtered <- yardstick::recall_vec(estimate = factor(ensembleCombinedFilteredGBC$GBC_max_col[idx],
                                                                 levels = c("Pre", "MVA", "MPXV")),
                                               truth = factor(ensembleCombinedFilteredGBC$real[idx],
                                                              levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  precision_GBC_filtered <- yardstick::precision_vec(estimate = factor(ensembleCombinedFilteredGBC$GBC_max_col[idx],
                                                                       levels = c("Pre", "MVA", "MPXV")),
                                                     truth = factor(ensembleCombinedFilteredGBC$real[idx],
                                                                    levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  sens_GBC_filtered <- yardstick::sens_vec(estimate = factor(ensembleCombinedFilteredGBC$GBC_max_col[idx],
                                                             levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(ensembleCombinedFilteredGBC$real[idx],
                                                          levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  spec_GBC_filtered <- yardstick::spec_vec(estimate = factor(ensembleCombinedFilteredGBC$GBC_max_col[idx],
                                                             levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(ensembleCombinedFilteredGBC$real[idx],
                                                          levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  f1_model_GBC_filtered[i] <-  f1_GBC_filtered
  recall_model_GBC_filtered[i] <- recall_GBC_filtered
  precision_model_GBC_filtered[i] <- precision_GBC_filtered
  sens_model_GBC_filtered[i] <- sens_GBC_filtered
  spec_model_GBC_filtered[i] <- spec_GBC_filtered
}

# 5: Bootstrap for RF filtered
for(i in 1:n_boot) {
  idx <- sample(1:n_val_RF_filtered, size = n_val_RF_filtered, replace = TRUE)
  f1_RF_filtered <- yardstick::f_meas_vec(estimate = factor(ensembleCombinedFilteredRF$RF_max_col[idx],
                                                            levels = c("Pre", "MVA", "MPXV")),
                                          truth = factor(ensembleCombinedFilteredRF$real[idx],
                                                         levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  recall_RF_filtered <- yardstick::recall_vec(estimate = factor(ensembleCombinedFilteredRF$RF_max_col[idx],
                                                                levels = c("Pre", "MVA", "MPXV")),
                                              truth = factor(ensembleCombinedFilteredRF$real[idx],
                                                             levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  precision_RF_filtered <- yardstick::precision_vec(estimate = factor(ensembleCombinedFilteredRF$RF_max_col[idx],
                                                                      levels = c("Pre", "MVA", "MPXV")),
                                                    truth = factor(ensembleCombinedFilteredRF$real[idx],
                                                                   levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  sens_RF_filtered <- yardstick::sens_vec(estimate = factor(ensembleCombinedFilteredRF$RF_max_col[idx],
                                                            levels = c("Pre", "MVA", "MPXV")),
                                          truth = factor(ensembleCombinedFilteredRF$real[idx],
                                                         levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  spec_RF_filtered <- yardstick::spec_vec(estimate = factor(ensembleCombinedFilteredRF$RF_max_col[idx],
                                                            levels = c("Pre", "MVA", "MPXV")),
                                          truth = factor(ensembleCombinedFilteredRF$real[idx],
                                                         levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  f1_model_RF_filtered[i] <-  f1_RF_filtered
  recall_model_RF_filtered[i] <- recall_RF_filtered
  precision_model_RF_filtered[i] <- precision_RF_filtered
  sens_model_RF_filtered[i] <- sens_RF_filtered
  spec_model_RF_filtered[i] <- spec_RF_filtered
}

# 6: Bootstrap for unfiltered data combined
for(i in 1:n_boot) {
  idx <- sample(1:n_val, size = n_val, replace = TRUE)
  # 6: Serostatus
  f1_serostatus <- yardstick::f_meas_vec(estimate = factor(ensembleCombined$max_col[idx],
                                                           levels = c("Pre", "MVA", "MPXV")),
                                         truth = factor(ensembleCombined$real[idx],
                                                        levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  recall_serostatus <- yardstick::recall_vec(estimate = factor(ensembleCombined$max_col[idx],
                                                               levels = c("Pre", "MVA", "MPXV")),
                                             truth = factor(ensembleCombined$real[idx],
                                                            levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  precision_serostatus <- yardstick::precision_vec(estimate = factor(ensembleCombined$max_col[idx],
                                                                     levels = c("Pre", "MVA", "MPXV")),
                                                   truth = factor(ensembleCombined$real[idx],
                                                                  levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  sens_serostatus <- yardstick::sens_vec(estimate = factor(ensembleCombined$max_col[idx],
                                                           levels = c("Pre", "MVA", "MPXV")),
                                         truth = factor(ensembleCombined$real[idx],
                                                        levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  spec_serostatus <- yardstick::spec_vec(estimate = factor(ensembleCombined$max_col[idx],
                                                           levels = c("Pre", "MVA", "MPXV")),
                                         truth = factor(ensembleCombined$real[idx],
                                                        levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  f1_model_serostatus[i] <-  f1_serostatus
  recall_model_serostatus[i] <- recall_serostatus
  precision_model_serostatus[i] <- precision_serostatus
  sens_model_serostatus[i] <- sens_serostatus
  spec_model_serostatus[i] <- spec_serostatus
  
  # 7: Mean
  f1_mean <- yardstick::f_meas_vec(estimate = factor(ensembleCombined$mean_max_col[idx],
                                                     levels = c("Pre", "MVA", "MPXV")),
                                   truth = factor(ensembleCombined$real[idx],
                                                  levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  recall_mean <- yardstick::recall_vec(estimate = factor(ensembleCombined$mean_max_col[idx],
                                                         levels = c("Pre", "MVA", "MPXV")),
                                       truth = factor(ensembleCombined$real[idx],
                                                      levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  precision_mean <- yardstick::precision_vec(estimate = factor(ensembleCombined$mean_max_col[idx],
                                                               levels = c("Pre", "MVA", "MPXV")),
                                             truth = factor(ensembleCombined$real[idx],
                                                            levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  sens_mean <- yardstick::sens_vec(estimate = factor(ensembleCombined$mean_max_col[idx],
                                                     levels = c("Pre", "MVA", "MPXV")),
                                   truth = factor(ensembleCombined$real[idx],
                                                  levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  spec_mean <- yardstick::spec_vec(estimate = factor(ensembleCombined$mean_max_col[idx],
                                                     levels = c("Pre", "MVA", "MPXV")),
                                   truth = factor(ensembleCombined$real[idx],
                                                  levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  f1_model_mean[i] <-  f1_mean
  recall_model_mean[i] <- recall_mean
  precision_model_mean[i] <- precision_mean
  sens_model_mean[i] <- sens_mean
  spec_model_mean[i] <- spec_mean
  
  # 8: LDA
  f1_LDA <- yardstick::f_meas_vec(estimate = factor(ensembleCombined$LDA_max_col[idx],
                                                    levels = c("Pre", "MVA", "MPXV")),
                                  truth = factor(ensembleCombined$real[idx],
                                                 levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  recall_LDA <- yardstick::recall_vec(estimate = factor(ensembleCombined$LDA_max_col[idx],
                                                        levels = c("Pre", "MVA", "MPXV")),
                                      truth = factor(ensembleCombined$real[idx],
                                                     levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  precision_LDA <- yardstick::precision_vec(estimate = factor(ensembleCombined$LDA_max_col[idx],
                                                              levels = c("Pre", "MVA", "MPXV")),
                                            truth = factor(ensembleCombined$real[idx],
                                                           levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  sens_LDA <- yardstick::sens_vec(estimate = factor(ensembleCombined$LDA_max_col[idx],
                                                    levels = c("Pre", "MVA", "MPXV")),
                                  truth = factor(ensembleCombined$real[idx],
                                                 levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  spec_LDA <- yardstick::spec_vec(estimate = factor(ensembleCombined$LDA_max_col[idx],
                                                    levels = c("Pre", "MVA", "MPXV")),
                                  truth = factor(ensembleCombined$real[idx],
                                                 levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  f1_model_LDA[i] <-  f1_LDA
  recall_model_LDA[i] <- recall_LDA
  precision_model_LDA[i] <- precision_LDA
  sens_model_LDA[i] <- sens_LDA
  spec_model_LDA[i] <- spec_LDA
  
  #9: GBC
  f1_GBC <- yardstick::f_meas_vec(estimate = factor(ensembleCombined$GBC_max_col[idx],
                                                    levels = c("Pre", "MVA", "MPXV")),
                                  truth = factor(ensembleCombined$real[idx],
                                                 levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  recall_GBC <- yardstick::recall_vec(estimate = factor(ensembleCombined$GBC_max_col[idx],
                                                        levels = c("Pre", "MVA", "MPXV")),
                                      truth = factor(ensembleCombined$real[idx],
                                                     levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  precision_GBC <- yardstick::precision_vec(estimate = factor(ensembleCombined$GBC_max_col[idx],
                                                              levels = c("Pre", "MVA", "MPXV")),
                                            truth = factor(ensembleCombined$real[idx],
                                                           levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  sens_GBC <- yardstick::sens_vec(estimate = factor(ensembleCombined$GBC_max_col[idx],
                                                    levels = c("Pre", "MVA", "MPXV")),
                                  truth = factor(ensembleCombined$real[idx],
                                                 levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  spec_GBC <- yardstick::spec_vec(estimate = factor(ensembleCombined$GBC_max_col[idx],
                                                    levels = c("Pre", "MVA", "MPXV")),
                                  truth = factor(ensembleCombined$real[idx],
                                                 levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  f1_model_GBC[i] <-  f1_GBC
  recall_model_GBC[i] <- recall_GBC
  precision_model_GBC[i] <- precision_GBC
  sens_model_GBC[i] <- sens_GBC
  spec_model_GBC[i] <- spec_GBC
  
  # 10: RF
  f1_RF <- yardstick::f_meas_vec(estimate = factor(ensembleCombined$RF_max_col[idx],
                                                   levels = c("Pre", "MVA", "MPXV")),
                                 truth = factor(ensembleCombined$real[idx],
                                                levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  recall_RF <- yardstick::recall_vec(estimate = factor(ensembleCombined$RF_max_col[idx],
                                                       levels = c("Pre", "MVA", "MPXV")),
                                     truth = factor(ensembleCombined$real[idx],
                                                    levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  precision_RF <- yardstick::precision_vec(estimate = factor(ensembleCombined$RF_max_col[idx],
                                                             levels = c("Pre", "MVA", "MPXV")),
                                           truth = factor(ensembleCombined$real[idx],
                                                          levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  sens_RF <- yardstick::sens_vec(estimate = factor(ensembleCombined$RF_max_col[idx],
                                                   levels = c("Pre", "MVA", "MPXV")),
                                 truth = factor(ensembleCombined$real[idx],
                                                levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  spec_RF <- yardstick::spec_vec(estimate = factor(ensembleCombined$RF_max_col[idx],
                                                   levels = c("Pre", "MVA", "MPXV")),
                                 truth = factor(ensembleCombined$real[idx],
                                                levels = c("Pre", "MVA", "MPXV")), estimator = "macro")
  f1_model_RF[i] <-  f1_RF
  recall_model_RF[i] <- recall_RF
  precision_model_RF[i] <- precision_RF
  sens_model_RF[i] <- sens_RF
  spec_model_RF[i] <- spec_RF
}

# End implementation boostrap analysis
####

####
# Implement per class sensitivity and specificity for GBC only

# Define lists
sens_model_MPXV_GBC_Filtered <- numeric(n_boot)
sens_model_MVA_GBC_Filtered <- numeric(n_boot)
sens_model_Pre_GBC_Filtered <-  numeric(n_boot)

spec_model_MPXV_GBC_Filtered <- numeric(n_boot)
spec_model_MVA_GBC_Filtered <-  numeric(n_boot)
spec_model_Pre_GBC_Filtered <-  numeric(n_boot)

sens_model_MPXV_GBC <- numeric(n_boot)
sens_model_MVA_GBC <- numeric(n_boot)
sens_model_Pre_GBC <- numeric(n_boot)

spec_model_MPXV_GBC <- numeric(n_boot)
spec_model_MVA_GBC <- numeric(n_boot)
spec_model_Pre_GBC <- numeric(n_boot)

for(i in 1:n_boot) {
  idx <- sample(1:n_val_GBC_filtered, size = n_val_GBC_filtered, replace = TRUE)
  confusionMatrixGBCFiltered <- 
    confusionMatrix(factor(ensembleCombinedFilteredGBC$real)[idx], factor(ensembleCombinedFilteredGBC$GBC_max_col)[idx])
  
  sens_MPXV_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Sensitivity"][1]
  sens_MVA_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Sensitivity"][2]
  sens_Pre_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Sensitivity"][3]
  
  spec_MPXV_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Specificity"][1]
  spec_MVA_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Specificity"][2]
  spec_Pre_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Specificity"][3]
  
  sens_model_MPXV_GBC_Filtered[i] <-  sens_MPXV_GBC_Filtered
  sens_model_MVA_GBC_Filtered[i] <- sens_MVA_GBC_Filtered
  sens_model_Pre_GBC_Filtered[i] <- sens_Pre_GBC_Filtered 
  
  spec_model_MPXV_GBC_Filtered[i] <-  spec_MPXV_GBC_Filtered
  spec_model_MVA_GBC_Filtered[i] <- spec_MVA_GBC_Filtered
  spec_model_Pre_GBC_Filtered[i] <- spec_Pre_GBC_Filtered
  
}

for(i in 1:n_boot) {
  idx <- sample(1:n_val, size = n_val, replace = TRUE)
  confusionMatrixGBC <- 
    confusionMatrix(factor(ensembleCombined$real)[idx], factor(ensembleCombined$GBC_max_col)[idx])
  
  sens_MPXV_GBC <- confusionMatrixGBC[["byClass"]][, "Sensitivity"][1]
  sens_MVA_GBC <- confusionMatrixGBC[["byClass"]][, "Sensitivity"][2]
  sens_Pre_GBC <- confusionMatrixGBC[["byClass"]][, "Sensitivity"][3]
  
  spec_MPXV_GBC <- confusionMatrixGBC[["byClass"]][, "Specificity"][1]
  spec_MVA_GBC <- confusionMatrixGBC[["byClass"]][, "Specificity"][2]
  spec_Pre_GBC <- confusionMatrixGBC[["byClass"]][, "Specificity"][3]
  
  sens_model_MPXV_GBC[i] <-  sens_MPXV_GBC
  sens_model_MVA_GBC[i] <- sens_MVA_GBC
  sens_model_Pre_GBC[i] <- sens_Pre_GBC 
  
  spec_model_MPXV_GBC[i] <-  spec_MPXV_GBC
  spec_model_MVA_GBC[i] <- spec_MVA_GBC
  spec_model_Pre_GBC[i] <- spec_Pre_GBC
}
# End implementation of per class sensitivity and specificity for GBC only
####


####
# Generate table with F1-value output
f1_names <- ls(pattern = "^f1_model")
precision_names <- ls(pattern = "^precision_model")
recall_names <- ls(pattern = "^recall_model")
sens_names <- ls(pattern = "^sens_model")
spec_names <- ls(pattern = "^spec_model")

# Calculate mean and confidence intervalls
get_f1_stats <- function(x) {
  v <- get(x)
  c(
    mean = mean(v),
    ci_lower = quantile(v, 0.025),
    ci_upper = quantile(v, 0.975)
  )
}

# Write data into one table
f1_stats <- sapply(f1_names, get_f1_stats)
f1_stats_df <- as.data.frame(t(f1_stats))
f1_stats_df$Model <- rownames(f1_stats_df)
rownames(f1_stats_df) <- NULL

precision_stats <- sapply(precision_names, get_f1_stats)
precision_stats_df <- as.data.frame(t(precision_stats))
precision_stats_df$Model <- rownames(precision_stats_df)
rownames(precision_stats_df) <- NULL

recall_stats <- sapply(recall_names, get_f1_stats)
recall_stats_df <- as.data.frame(t(recall_stats))
recall_stats_df$Model <- rownames(recall_stats_df)
rownames(recall_stats_df) <- NULL

sens_stats <- sapply(sens_names, get_f1_stats)
sens_stats_df <- as.data.frame(t(sens_stats))
sens_stats_df$Model <- rownames(sens_stats_df)
rownames(sens_stats_df) <- NULL

spec_stats <- sapply(spec_names, get_f1_stats)
spec_stats_df <- as.data.frame(t(spec_stats))
spec_stats_df$Model <- rownames(spec_stats_df)
rownames(spec_stats_df) <- NULL





f1_summary <- f1_stats_df %>% 
  mutate(Model = str_remove(Model, "f1_model_")) %>% 
  select(Model, mean, ci_lower = `ci_lower.2.5%`, ci_upper = `ci_upper.97.5%`)
f1_summary$Model <- factor(f1_summary$Model, levels = f1_summary$Model[order(f1_summary$mean, decreasing = TRUE)])
f1_summary$F1_formatted <- sprintf(
  "%.2f (%.2f–%.2f)",
  f1_summary$mean,
  f1_summary$ci_lower,
  f1_summary$ci_upper
)

precision_summary <- precision_stats_df %>% 
  mutate(Model = str_remove(Model, "precision_model_")) %>% 
  select(Model, mean, ci_lower = `ci_lower.2.5%`, ci_upper = `ci_upper.97.5%`)
precision_summary$Model <- factor(precision_summary$Model, levels = precision_summary$Model[order(precision_summary$mean, decreasing = TRUE)])
precision_summary$precision_formatted <- sprintf(
  "%.2f (%.2f–%.2f)",
  precision_summary$mean,
  precision_summary$ci_lower,
  precision_summary$ci_upper
)


recall_summary <- recall_stats_df %>% 
  mutate(Model = str_remove(Model, "recall_model_")) %>% 
  select(Model, mean, ci_lower = `ci_lower.2.5%`, ci_upper = `ci_upper.97.5%`)
recall_summary$Model <- factor(recall_summary$Model, levels = recall_summary$Model[order(recall_summary$mean, decreasing = TRUE)])
recall_summary$recall_formatted <- sprintf(
  "%.2f (%.2f–%.2f)",
  recall_summary$mean,
  recall_summary$ci_lower,
  recall_summary$ci_upper
)


sens_summary <- sens_stats_df %>% 
  mutate(Model = str_remove(Model, "sens_model_")) %>% 
  select(Model, mean, ci_lower = `ci_lower.2.5%`, ci_upper = `ci_upper.97.5%`)
sens_summary$Model <- factor(sens_summary$Model, levels = sens_summary$Model[order(sens_summary$mean, decreasing = TRUE)])
sens_summary$sens_formatted <- sprintf(
  "%.2f (%.2f–%.2f)",
  sens_summary$mean,
  sens_summary$ci_lower,
  sens_summary$ci_upper
)


spec_summary <- spec_stats_df %>% 
  mutate(Model = str_remove(Model, "spec_model_")) %>% 
  select(Model, mean, ci_lower = `ci_lower.2.5%`, ci_upper = `ci_upper.97.5%`)
spec_summary$Model <- factor(spec_summary$Model, levels = spec_summary$Model[order(spec_summary$mean, decreasing = TRUE)])
spec_summary$spec_formatted <- sprintf(
  "%.2f (%.2f–%.2f)",
  spec_summary$mean,
  spec_summary$ci_lower,
  spec_summary$ci_upper
)


table_performance_establishment <- 
  f1_summary %>% 
  left_join(precision_summary, by = "Model", suffix = c("", "_prec")) %>% 
  left_join(recall_summary, by = "Model", suffix = c("", "_rec")) %>% 
  left_join(sens_summary, by = "Model", suffix = c("", "_rec")) %>% 
  left_join(spec_summary, by = "Model", suffix = c("", "_rec"))

save(table_performance_establishment, file = "output/table_performance_establishment_Rev_04.Rdata")


# Generate table with sensitivites and specificities by class for GBC
table_GBC_class_establishment <- 
sens_summary %>% 
  left_join(spec_summary, by = "Model", suffix = c("", "_rec")) %>% 
  filter(grepl("^MVA_|^MPXV|^Pre", Model)) %>% 
  select(Model, sensitivity = sens_formatted, specificity = spec_formatted)
save(table_GBC_class_establishment, file = "output/table_GBC_class_establishment_Rev_04.Rdata")

# Plot: Means and CIs as points and error bars
ggplot(f1_summary, aes(x = Model, y = mean)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "grey40") +
  ylim(0.70, 0.90) +
  labs(
    title = "Model Comparison: F1 Scores with 95% Confidence Intervals",
    y = "F1 Score (mean ± 95% CI)",
    x = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(precision_summary, aes(x = Model, y = mean)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "grey40") +
  ylim(0.70, 0.90) +
  labs(
    title = "Model Comparison: Precision with 95% Confidence Intervals",
    y = "Precision (mean ± 95% CI)",
    x = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(recall_summary, aes(x = Model, y = mean)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "grey40") +
  ylim(0.70, 0.90) +
  labs(
    title = "Model Comparison: Recall with 95% Confidence Intervals",
    y = "Recall (mean ± 95% CI)",
    x = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Plot confusion matrices for filtered data
plotTest <-
  plot_confusion_matrix(conf_mat,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = TRUE,
                        add_zero_shading = TRUE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges") +
  ggtitle("Ensemble serostatus") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotTestLDA <-
  plot_confusion_matrix(conf_mat_LDA,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = FALSE,
                        add_zero_shading = FALSE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges") +
  ggtitle("LDA") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotTestGBC <-
  plot_confusion_matrix(conf_mat_GBC,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = FALSE,
                        add_zero_shading = FALSE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges") +
  ggtitle("GBC") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotTestMean <-
  plot_confusion_matrix(conf_mat_mean,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = FALSE,
                        add_zero_shading = FALSE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges") +
  ggtitle("Ensemble mean") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

## Rev 02: Plot RF with filtered data
plotTestRF <-
  plot_confusion_matrix(conf_mat_RF,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = FALSE,
                        add_zero_shading = FALSE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges") +
  ggtitle("RF") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



plotConfMatrixTest <-
  ggarrange(plotTestLDA, plotTestGBC, plotTestRF, ncol = 3, 
            align = "hv")

ggsave("output/plotConfusionTest_Rev_02.pdf", plotConfMatrixTest,
       width = 12, height = 4, dpi = 600)

ggsave("output/plotConfusionTest_Rev_02.png", plotConfMatrixTest,
       width = 12, height = 4, dpi = 600)


# Plot confusion matrices for all data
plotTestAll <-
  plot_confusion_matrix(conf_mat_all,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = TRUE,
                        add_zero_shading = TRUE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges") +
  ggtitle("Ensemble serostatus") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotTestLDAAll <-
  plot_confusion_matrix(conf_mat_LDA_all,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = FALSE,
                        add_zero_shading = FALSE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges") +
  ggtitle("LDA") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotTestGBCAll <-
  plot_confusion_matrix(conf_mat_GBC_all,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = FALSE,
                        add_zero_shading = FALSE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges") +
  ggtitle("GBC") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotTestMeanAll <-
  plot_confusion_matrix(conf_mat_mean_all,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = FALSE,
                        add_zero_shading = FALSE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges") +
  ggtitle("Ensemble mean") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

## Rev 02: Plot RF for all data
plotTestRFAll <-
  plot_confusion_matrix(conf_mat_RF_all,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = FALSE,
                        add_zero_shading = FALSE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges") +
  ggtitle("RF") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



plotConfMatrixTestAll <-
  ggarrange(plotTestLDAAll, plotTestGBCAll, plotTestRFAll, ncol = 3,
            align = "hv")

ggsave("output/plotConfusionTestAll_Rev_02.pdf", plotConfMatrixTestAll,
       width = 12, height = 4, dpi = 600)

ggsave("output/plotConfusionTestAll_Rev_02.png", plotConfMatrixTestAll,
       width = 12, height = 4, dpi = 600)
