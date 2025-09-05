####
# Analyse revision measurements
# Daniel Stern
# 2025-09-05
# Version 5.0
# Changes to version 2.0
# - Inclusion of RF alogrithm as baseline model
# Changes to version 3.0
# - Bootstrap analysis to calculate mean and 95% CI 
# Changes to version 4.0
# - Read classwise F1, precision and recall for GBC
# Changes to version 5.0
# - Show confusion matrix for all data
####

rm(list = ls(all.names = TRUE))

library(rio)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(gplots)
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(scales)
library(caret)
library(pROC)
library(rstatix)
library(corrplot)
library(cvms)
library(ggimage)
library(rsvg)
library(circlize)
library(yardstick)

# Set seed for bootstrap analysis
set.seed(123)

# Load datainput
load("input/ensemblePrediction_Rev_02.Rdata")
load("input/heatmap_input.Rdata")
load("input/dataInputComparison.Rdata")

heatmap_input_IgG <- 
  heatmap_input_IgG %>% 
  mutate(panel = factor(panel, levels = c("MPXV", "Pre", "MVA"), 
                        ordered = TRUE))

plotPreGroupedLDA <-
  ensemblePrediction %>%
  mutate(serostatus.delta = factor(serostatus.delta, levels = c(0,1),
                                   labels = c("negative", "positive"))) %>% 
  pivot_longer(cols = c(LDA_Pre, LDA_MVA, LDA_MPXV), names_to = "pred",
               values_to = "ensembl") %>% 
  mutate(pred = factor(pred, levels = c("LDA_Pre", "LDA_MVA", "LDA_MPXV"),
                       labels = c("Pre", "MVA", "MPXV"),
                       ordered = TRUE), 
         year_group = ifelse(year_birth < 1972, "<1972", ">=1972"),
         year_group = factor(year_group, levels = c("<1972", ">=1972"),
                             ordered = TRUE)) %>% 
  filter(real == "Pre") %>%
  filter(!is.na(ensembl)) %>% 
  group_by(year_group, pred) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(year_group) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ggplot(aes(x = year_group, y = percent , fill = pred)) +
  geom_bar(stat = "identity", position = "fill") +  # relative Häufigkeit (100%)
  geom_text(aes(label = count),
            position = position_fill(vjust = 0.5),  # Zahlen mittig im Balken
            size = 5, color = "white", fontface = "bold") +
  scale_fill_manual(values = colorblind_pal()(8)[2:8]) + # Farben für die Kategorien
  labs(
    x = "Year of birth",
    y = "Relative frequency",
    fill = "Prediction"
  ) +
  scale_y_continuous() +
  theme_pubr() +
  ggtitle("LDA") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        strip.background = element_blank(),
        legend.position = "none")


plotPreGroupedGBC <-
  ensemblePrediction %>%
  mutate(serostatus.delta = factor(serostatus.delta, levels = c(0,1),
                                   labels = c("negative", "positive"))) %>% 
  pivot_longer(cols = c(GBC_Pre, GBC_MVA, GBC_MPXV), names_to = "pred",
               values_to = "ensembl") %>% 
  mutate(pred = factor(pred, levels = c("GBC_Pre", "GBC_MVA", "GBC_MPXV"),
                       labels = c("Pre", "MVA", "MPXV"),
                       ordered = TRUE), 
         year_group = ifelse(year_birth < 1972, "<1972", ">=1972"),
         year_group = factor(year_group, levels = c("<1972", ">=1972"),
                             ordered = TRUE)) %>% 
  filter(real == "Pre") %>%
  filter(!is.na(ensembl)) %>% 
  group_by(year_group, pred) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(year_group) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ggplot(aes(x = year_group, y = percent , fill = pred)) +
  geom_bar(stat = "identity", position = "fill") +  # relative Häufigkeit (100%)
  geom_text(aes(label = count),
            position = position_fill(vjust = 0.5),  # Zahlen mittig im Balken
            size = 5, color = "white", fontface = "bold") +
  scale_fill_manual(values = colorblind_pal()(8)[2:8]) + # Farben für die Kategorien
  labs(
    x = "Year of birth",
    y = "Relative frequency",
    fill = "Prediction"
  ) +
  scale_y_continuous() +
  theme_pubr() +
  ggtitle("GBC") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        strip.background = element_blank(),
        legend.position = "none")

## Rev 02_ plot pre grouped RF
plotPreGroupedRF <-
  ensemblePrediction %>%
  mutate(serostatus.delta = factor(serostatus.delta, levels = c(0,1),
                                   labels = c("negative", "positive"))) %>% 
  pivot_longer(cols = c(RF_Pre, RF_MVA, RF_MPXV), names_to = "pred",
               values_to = "ensembl") %>% 
  mutate(pred = factor(pred, levels = c("RF_Pre", "RF_MVA", "RF_MPXV"),
                       labels = c("Pre", "MVA", "MPXV"),
                       ordered = TRUE), 
         year_group = ifelse(year_birth < 1972, "<1972", ">=1972"),
         year_group = factor(year_group, levels = c("<1972", ">=1972"),
                             ordered = TRUE)) %>% 
  filter(real == "Pre") %>%
  filter(!is.na(ensembl)) %>% 
  group_by(year_group, pred) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(year_group) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ggplot(aes(x = year_group, y = percent , fill = pred)) +
  geom_bar(stat = "identity", position = "fill") +  # relative Häufigkeit (100%)
  geom_text(aes(label = count),
            position = position_fill(vjust = 0.5),  # Zahlen mittig im Balken
            size = 5, color = "white", fontface = "bold") +
  scale_fill_manual(values = colorblind_pal()(8)[2:8]) + # Farben für die Kategorien
  labs(
    x = "Year of birth",
    y = "Relative frequency",
    fill = "Prediction"
  ) +
  scale_y_continuous() +
  theme_pubr() +
  ggtitle("RF") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        strip.background = element_blank(),
        legend.position = "none")

####
# Plot heatmap of raw data
col_fun_MVA <- colorRamp2(c(0.5, 1), c("white", colorblind_pal()(8)[3]))
col_fun_Pre <- colorRamp2(c(0.5, 1), c("white", colorblind_pal()(8)[2]))
col_fun_MPXV <- colorRamp2(c(0.5, 1), c("white", colorblind_pal()(8)[4]))

matIgG <- as.matrix(scale(heatmap_input_IgG[2:16]))
haIgG <- HeatmapAnnotation(df = data.frame(panel = heatmap_input_IgG$panel),
                           annotation_height = unit(4, "mm"))

pdf(file = "output/heatmap_IgG_white_Rev_02.pdf", width = 14, height = 10)
Heatmap(matIgG, split = factor(heatmap_input_IgG$panel, levels = c("Pre", "MVA", "MPXV" ), ordered = TRUE),
        name = "Multiplex quantfied", clustering_method_columns = "ward.D2", 
        clustering_method_rows = "ward.D2") +
  Heatmap(factor(heatmap_input_IgG$serostatus.delta, levels = c(0,1),
                 labels = c("negative", "positive"),
                 ordered = TRUE), name = "Delta-VACV", col = c(colorblind_pal()(8)[c(6,7)])) +
  
  Heatmap(heatmap_input_IgG$n_MVA_vacc, name = "MVA vaccinations", na_col = "white", col = c("grey", viridis_pal()(5)[2:5])) +
  Heatmap(heatmap_input_IgG$LDA_MVA, name = "MVA LDA", col = col_fun_MVA, na_col = "white") +
  Heatmap(heatmap_input_IgG$GBC_MVA, name = "MVA GBC", col = col_fun_MVA, na_col = "white") +
  Heatmap(heatmap_input_IgG$RF_MVA, name = "MVA RF", col = col_fun_MVA, na_col = "white") +
  #  Heatmap(heatmap_input_IgG$MVA, name = "MVA ensemble mean", col = col_fun_MVA, na_col = "white") +
  #  Heatmap(heatmap_input_IgG$Strat_MVA, name = "MVA ensemble serostatus", col = col_fun_MVA, na_col = "white") +
  Heatmap(heatmap_input_IgG$LDA_Pre, name = "Pre LDA", col = col_fun_Pre, na_col = "white") +
  Heatmap(heatmap_input_IgG$GBC_Pre, name = "Pre GBC", col = col_fun_Pre, na_col = "white") +
  Heatmap(heatmap_input_IgG$RF_Pre, name = "Pre RF", col = col_fun_Pre, na_col = "white") +
  #  Heatmap(heatmap_input_IgG$Pre, name = "Pre ensemble mean", col = col_fun_Pre, na_col = "white") +
  #  Heatmap(heatmap_input_IgG$Strat_Pre, name = "Pre ensemble serostatus", col = col_fun_Pre, na_col = "white") +
  Heatmap(heatmap_input_IgG$LDA_MPXV, name = "MPXV LDA", col = col_fun_MPXV, na_col = "white") +
  Heatmap(heatmap_input_IgG$GBC_MPXV, name = "MPXV GBC", col = col_fun_MPXV, na_col = "white") +
  Heatmap(heatmap_input_IgG$RF_MPXV, name = "MPXV RF", col = col_fun_MPXV, na_col = "white") +
  #  Heatmap(heatmap_input_IgG$MPXV, name = "MPXV ensemble mean", col = col_fun_MPXV, na_col = "white") +
  #  Heatmap(heatmap_input_IgG$Strat_MPXV, name = "MPXV ensemble serostatus", col = col_fun_MPXV, na_col = "white") +
  Heatmap(factor(heatmap_input_IgG$panel, levels = c("Pre", "MVA", "MPXV" ), 
                 ordered = TRUE), name = "Panel", col = c(colorblind_pal()(8)[2:4]))
dev.off()


matIgM <- as.matrix(scale(heatmap_input_IgM[2:16]))
ha <- HeatmapAnnotation(df = data.frame(panel = heatmap_input_IgM$panel),
                        annotation_height = unit(4, "mm"))
pdf(file = "output/heatmap_IgM_white_Rev_02.pdf", width = 14, height = 8)
Heatmap(matIgM, split = factor(heatmap_input_IgM$panel, levels = c("Pre", "MVA", "MPXV" ), ordered = TRUE),
        name = "Multiplex quantfied", clustering_method_columns = "ward.D2", 
        clustering_method_rows = "ward.D2") +
  Heatmap(factor(heatmap_input_IgM$serostatus.delta, levels = c(0,1),
                 labels = c("negative", "positive"),
                 ordered = TRUE), name = "Delta-VACV", col = c(colorblind_pal()(8)[c(6,7)])) +
  Heatmap(heatmap_input_IgM$n_MVA_vacc, name = "MVA vaccinations", na_col = "white", col = c("grey", viridis_pal()(5)[2:5])) +
  
  Heatmap(heatmap_input_IgM$LDA_MVA, name = "MVA LDA", col = c("white", colorblind_pal()(8)[3]), na_col = "white") +
  Heatmap(heatmap_input_IgM$GBC_MVA, name = "MVA GBC", col = c("white", colorblind_pal()(8)[3]), na_col = "white") +
  Heatmap(heatmap_input_IgM$RF_MVA, name = "MVA RF", col = c("white", colorblind_pal()(8)[3]), na_col = "white") +
  # Heatmap(heatmap_input_IgM$MVA, name = "MVA ensemble mean", col = c("white", colorblind_pal()(8)[3]), na_col = "white") +
  #  Heatmap(heatmap_input_IgM$Strat_MVA, name = "MVA ensemble serostatus", col = c("white", colorblind_pal()(8)[3]), na_col = "white") +
  Heatmap(heatmap_input_IgM$LDA_Pre, name = "Pre LDA", col = c("white", colorblind_pal()(8)[2]), na_col = "white") +
  Heatmap(heatmap_input_IgM$GBC_Pre, name = "Pre GBC", col = c("white", colorblind_pal()(8)[2]), na_col = "white") +
  Heatmap(heatmap_input_IgM$RF_Pre, name = "Pre RF", col = c("white", colorblind_pal()(8)[2]), na_col = "white") +
  # Heatmap(heatmap_input_IgM$Pre, name = "Pre ensemble mean", col = c("white", colorblind_pal()(8)[2]), na_col = "white") +
  #  Heatmap(heatmap_input_IgM$Strat_Pre, name = "Pre ensemble serostatus", col = c("white", colorblind_pal()(8)[2]), na_col = "white") +
  Heatmap(heatmap_input_IgM$LDA_MPXV, name = "MPXV LDA", col = c("white", colorblind_pal()(8)[4]), na_col = "white") +
  Heatmap(heatmap_input_IgM$GBC_MPXV, name = "MPXV GBC", col = c("white", colorblind_pal()(8)[4]), na_col = "white") +
  Heatmap(heatmap_input_IgM$RF_MPXV, name = "MPXV RF", col = c("white", colorblind_pal()(8)[4]), na_col = "white") +
  # Heatmap(heatmap_input_IgM$MPXV, name = "MPXV ensemble mean", col = c("white", colorblind_pal()(8)[4]), na_col = "white") +
  #  Heatmap(heatmap_input_IgM$Strat_MPXV, name = "MPXV ensemble serostatus", col = c("white", colorblind_pal()(8)[4]), na_col = "white") +
  Heatmap(factor(heatmap_input_IgM$panel, levels = c("Pre", "MVA", "MPXV" ), 
                 ordered = TRUE), name = "Panel", col = c(colorblind_pal()(8)[2:4])) 
dev.off()

####
# ROC for discrimination between different panels based on ATI-N
rets <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv",
          "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")
rocInput <- 
  heatmap_input_IgG %>% 
  mutate(predictor = if_else(real == "MPXV", "infected", "non-infected"))
roc_ATI <-
  roc(factor(rocInput$predictor, levels = c("infected", "non-infected"), ordered = TRUE), rocInput$`ATI-N`, )
ci_roc_ATI <- ci.coords(roc_ATI, x="best", input = "threshold",  best.method = c("closest.topleft"), ret = rets,
                        best.policy = "random")

####
# Calculate performance with threshold of 0.5 on all data
# Selection based for ensemble prediction and LDA due to comparability between 
# Algorithms
ensembleSelectedThreshold <-
  ensemblePrediction %>% 
  filter(pred_ensemble > 0.5) %>% 
  filter(mean_mean_conf > 0.5)

confusionMatrixSelectedThreshold <- 
  confusionMatrix(factor(ensembleSelectedThreshold$max_col), factor(ensembleSelectedThreshold $real))
confusionMatrixSelectedThreshold[["byClass"]][ , "F1"]

misclassMVA <-
  ensembleSelectedThreshold %>% 
  filter(real == "MVA" & max_col == "MPXV")
misclassPre <-
  ensembleSelectedThreshold %>% 
  filter(real == "Pre" & max_col == "MPXV")

#### 
# Export prediction data for comparison with data used for comparison with
# binary classifiers
save(file = "output/ensemblePrediction_Val.Rdata", ensemblePrediction)

####
# Export supporting table containing the data

# Select ATI-N results for merge
dataInputJoin <- 
  dataInputComparison %>% 
  select(sampleID_meta, `ATI-N`)

ensemblePredictionExport <-
  ensemblePrediction %>% 
  left_join(dataInputJoin, by = c("sampleID_meta")) %>% 
  mutate(Sample = row_number(),
         mean_mean_conf = round(mean_mean_conf, digits = 2),
         pred_ensemble = round(pred_ensemble, digits = 2),
         mean_pred_ensemble = round(mean_pred_ensemble, digits = 2)) %>% 
  select(Sample, Real = real, Pred_LDA = LDA_max_col, 
         Pred_GBC = GBC_max_col, Pred_RF = RF_max_col, Pred_ensemble_mean = mean_max_col,
         Pred_ensemble_sero = max_col,
         `Serostatus (Delta)` = serostatus_cat.delta,
         `Serostatus ATI-N-CPXV` = `ATI-N`,
         `MVA vacc. (n)`= n_MVA_vacc, `Year Birth` = year_birth,
         `Ensemble confidence serostatus` = pred_ensemble, 
         `Ensemble confidence mean` = mean_pred_ensemble,
         `LDA confidence` = mean_mean_conf)

export(ensemblePredictionExport, file = "output/STableEnseblePrediction.xlsx")



####
# PLot confusion matrices using package cvsm and calculate performance data
# Ensemble prediction
conf_mat <- confusion_matrix(targets = factor(ensembleSelectedThreshold$real),
                             predictions = factor(ensembleSelectedThreshold$max_col))


# Ensemble prediction mean
conf_mat_mean <- confusion_matrix(targets = factor(ensembleSelectedThreshold$real),
                                  predictions = factor(ensembleSelectedThreshold$mean_max_col))



# LDA prediction
conf_mat_LDA <- confusion_matrix(targets = factor(ensembleSelectedThreshold$real),
                                 predictions = factor(ensembleSelectedThreshold$LDA_max_col))


# GBC prediction
conf_mat_GBC <- confusion_matrix(targets = factor(ensembleSelectedThreshold$real),
                                 predictions = factor(ensembleSelectedThreshold$GBC_max_col))

## Rev 02: RF prediction
conf_mat_RF <- confusion_matrix(targets = factor(ensembleSelectedThreshold$real),
                                predictions = factor(ensembleSelectedThreshold$RF_max_col))



####
# PLot confusion matrices using package cvsm and calculate performance data for all datapoint
# Ensemble prediction all
conf_mat_all <- confusion_matrix(targets = factor(ensemblePrediction$real),
                                 predictions = factor(ensemblePrediction$max_col))

# Ensemble prediction mean
conf_mat_mean_all <- confusion_matrix(targets = factor(ensemblePrediction$real),
                                      predictions = factor(ensemblePrediction$mean_max_col))

# LDA prediction
conf_mat_LDA_all <- confusion_matrix(targets = factor(ensemblePrediction$real),
                                     predictions = factor(ensemblePrediction$LDA_max_col))

# GBC prediction
conf_mat_GBC_all <- confusion_matrix(targets = factor(ensemblePrediction$real),
                                     predictions = factor(ensemblePrediction$GBC_max_col))

# Rev 02: LDA prediction
conf_mat_RF_all <- confusion_matrix(targets = factor(ensemblePrediction$real),
                                    predictions = factor(ensemblePrediction$RF_max_col))





####
# Bootstrap analysis to calculate mean and 95C% CI for the validation data
# Datainput and different cases
# save dataframe with new name 
ensembleCombinedFiltered <- 
  ensembleSelectedThreshold

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
n_val_filtered <- length(ensembleCombinedFiltered$RF_max_col) #5: 124
n_val <-length(ensemblePrediction$max_col) #6: 143

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

# Bootstrap for filtered combined
for(i in 1:n_boot) {
  idx <- sample(1:n_val_filtered, size = n_val_filtered, replace = TRUE)
  # 1: Serostatus
  f1_serostatus_filtered <- yardstick::f_meas_vec(estimate = factor(ensembleCombinedFiltered$max_col[idx],
                                                                    levels = c("MPXV", "MVA", "Pre")),
                                                  truth = factor(ensembleCombinedFiltered$real[idx],
                                                                 levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  recall_serostatus_filtered <- yardstick::recall_vec(estimate = factor(ensembleCombinedFiltered$max_col[idx],
                                                                        levels = c("MPXV", "MVA", "Pre")),
                                                      truth = factor(ensembleCombinedFiltered$real[idx],
                                                                     levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  precision_serostatus_filtered <- yardstick::precision_vec(estimate = factor(ensembleCombinedFiltered$max_col[idx],
                                                                              levels = c("MPXV", "MVA", "Pre")),
                                                            truth = factor(ensembleCombinedFiltered$real[idx],
                                                                           levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  sens_serostatus_filtered <- yardstick::sens_vec(estimate = factor(ensembleCombinedFiltered$max_col[idx],
                                                                    levels = c("MPXV", "MVA", "Pre")),
                                                  truth = factor(ensembleCombinedFiltered$real[idx],
                                                                 levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  spec_serostatus_filtered <- yardstick::spec_vec(estimate = factor(ensembleCombinedFiltered$max_col[idx],
                                                                    levels = c("MPXV", "MVA", "Pre")),
                                                  truth = factor(ensembleCombinedFiltered$real[idx],
                                                                 levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  
  
  f1_model_serostatus_filtered[i] <-  f1_serostatus_filtered
  recall_model_serostatus_filtered[i] <- recall_serostatus_filtered
  precision_model_serostatus_filtered[i] <- precision_serostatus_filtered
  sens_model_serostatus_filtered[i] <- sens_serostatus_filtered
  spec_model_serostatus_filtered[i] <- spec_serostatus_filtered
  
  # 2: Mean
  f1_mean_filtered <- yardstick::f_meas_vec(estimate = factor(ensembleCombinedFiltered$mean_max_col[idx],
                                                              levels = c("MPXV", "MVA", "Pre")),
                                            truth = factor(ensembleCombinedFiltered$real[idx],
                                                           levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  recall_mean_filtered <- yardstick::recall_vec(estimate = factor(ensembleCombinedFiltered$mean_max_col[idx],
                                                                  levels = c("MPXV", "MVA", "Pre")),
                                                truth = factor(ensembleCombinedFiltered$real[idx],
                                                               levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  precision_mean_filtered <- yardstick::precision_vec(estimate = factor(ensembleCombinedFiltered$mean_max_col[idx],
                                                                        levels = c("MPXV", "MVA", "Pre")),
                                                      truth = factor(ensembleCombinedFiltered$real[idx],
                                                                     levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  sens_mean_filtered <- yardstick::sens_vec(estimate = factor(ensembleCombinedFiltered$mean_max_col[idx],
                                                              levels = c("MPXV", "MVA", "Pre")),
                                            truth = factor(ensembleCombinedFiltered$real[idx],
                                                           levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  spec_mean_filtered <- yardstick::spec_vec(estimate = factor(ensembleCombinedFiltered$mean_max_col[idx],
                                                              levels = c("MPXV", "MVA", "Pre")),
                                            truth = factor(ensembleCombinedFiltered$real[idx],
                                                           levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  f1_model_mean_filtered[i] <-  f1_mean_filtered
  recall_model_mean_filtered[i] <- recall_mean_filtered
  precision_model_mean_filtered[i] <- precision_mean_filtered
  sens_model_mean_filtered[i] <- sens_mean_filtered
  spec_model_mean_filtered[i] <- spec_mean_filtered
  
  # 3: LDA
  f1_LDA_filtered <- yardstick::f_meas_vec(estimate = factor(ensembleCombinedFiltered$LDA_max_col[idx],
                                                             levels = c("MPXV", "MVA", "Pre")),
                                           truth = factor(ensembleCombinedFiltered$real[idx],
                                                          levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  recall_LDA_filtered <- yardstick::recall_vec(estimate = factor(ensembleCombinedFiltered$LDA_max_col[idx],
                                                                 levels = c("MPXV", "MVA", "Pre")),
                                               truth = factor(ensembleCombinedFiltered$real[idx],
                                                              levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  precision_LDA_filtered <- yardstick::precision_vec(estimate = factor(ensembleCombinedFiltered$LDA_max_col[idx],
                                                                       levels = c("MPXV", "MVA", "Pre")),
                                                     truth = factor(ensembleCombinedFiltered$real[idx],
                                                                    levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  sens_LDA_filtered <- yardstick::sens_vec(estimate = factor(ensembleCombinedFiltered$LDA_max_col[idx],
                                                             levels = c("MPXV", "MVA", "Pre")),
                                           truth = factor(ensembleCombinedFiltered$real[idx],
                                                          levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  spec_LDA_filtered <- yardstick::spec_vec(estimate = factor(ensembleCombinedFiltered$LDA_max_col[idx],
                                                             levels = c("MPXV", "MVA", "Pre")),
                                           truth = factor(ensembleCombinedFiltered$real[idx],
                                                          levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  f1_model_LDA_filtered[i] <-  f1_LDA_filtered
  recall_model_LDA_filtered[i] <- recall_LDA_filtered
  precision_model_LDA_filtered[i] <- precision_LDA_filtered
  sens_model_LDA_filtered[i] <- sens_LDA_filtered
  spec_model_LDA_filtered[i] <- spec_LDA_filtered
  
  # 4: GBC
  f1_GBC_filtered <- yardstick::f_meas_vec(estimate = factor(ensembleCombinedFiltered$GBC_max_col[idx],
                                                             levels = c("MPXV", "MVA", "Pre")),
                                           truth = factor(ensembleCombinedFiltered$real[idx],
                                                          levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  recall_GBC_filtered <- yardstick::recall_vec(estimate = factor(ensembleCombinedFiltered$GBC_max_col[idx],
                                                                 levels = c("MPXV", "MVA", "Pre")),
                                               truth = factor(ensembleCombinedFiltered$real[idx],
                                                              levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  precision_GBC_filtered <- yardstick::precision_vec(estimate = factor(ensembleCombinedFiltered$GBC_max_col[idx],
                                                                       levels = c("MPXV", "MVA", "Pre")),
                                                     truth = factor(ensembleCombinedFiltered$real[idx],
                                                                    levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  sens_GBC_filtered <- yardstick::sens_vec(estimate = factor(ensembleCombinedFiltered$GBC_max_col[idx],
                                                             levels = c("MPXV", "MVA", "Pre")),
                                           truth = factor(ensembleCombinedFiltered$real[idx],
                                                          levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  spec_GBC_filtered <- yardstick::spec_vec(estimate = factor(ensembleCombinedFiltered$GBC_max_col[idx],
                                                             levels = c("MPXV", "MVA", "Pre")),
                                           truth = factor(ensembleCombinedFiltered$real[idx],
                                                          levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  f1_model_GBC_filtered[i] <-  f1_GBC_filtered
  recall_model_GBC_filtered[i] <- recall_GBC_filtered
  precision_model_GBC_filtered[i] <- precision_GBC_filtered
  sens_model_GBC_filtered[i] <- sens_GBC_filtered
  spec_model_GBC_filtered[i] <- spec_GBC_filtered
  
  # 5: RF
  f1_RF_filtered <- yardstick::f_meas_vec(estimate = factor(ensembleCombinedFiltered$RF_max_col[idx],
                                                            levels = c("MPXV", "MVA", "Pre")),
                                          truth = factor(ensembleCombinedFiltered$real[idx],
                                                         levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  recall_RF_filtered <- yardstick::recall_vec(estimate = factor(ensembleCombinedFiltered$RF_max_col[idx],
                                                                levels = c("MPXV", "MVA", "Pre")),
                                              truth = factor(ensembleCombinedFiltered$real[idx],
                                                             levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  precision_RF_filtered <- yardstick::precision_vec(estimate = factor(ensembleCombinedFiltered$RF_max_col[idx],
                                                                      levels = c("MPXV", "MVA", "Pre")),
                                                    truth = factor(ensembleCombinedFiltered$real[idx],
                                                                   levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  sens_RF_filtered <- yardstick::sens_vec(estimate = factor(ensembleCombinedFiltered$RF_max_col[idx],
                                                            levels = c("MPXV", "MVA", "Pre")),
                                          truth = factor(ensembleCombinedFiltered$real[idx],
                                                         levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  spec_RF_filtered <- yardstick::spec_vec(estimate = factor(ensembleCombinedFiltered$RF_max_col[idx],
                                                            levels = c("MPXV", "MVA", "Pre")),
                                          truth = factor(ensembleCombinedFiltered$real[idx],
                                                         levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  f1_model_RF_filtered[i] <-  f1_RF_filtered
  recall_model_RF_filtered[i] <- recall_RF_filtered
  precision_model_RF_filtered[i] <- precision_RF_filtered
  sens_model_RF_filtered[i] <- sens_RF_filtered
  spec_model_RF_filtered[i] <- spec_RF_filtered
}

# Bootstrap for unfiltered data combined
for(i in 1:n_boot) {
  idx <- sample(1:n_val, size = n_val, replace = TRUE)
  # 6: Serostatus
  f1_serostatus <- yardstick::f_meas_vec(estimate = factor(ensemblePrediction$max_col[idx],
                                                           levels = c("MPXV", "MVA", "Pre")),
                                         truth = factor(ensemblePrediction$real[idx],
                                                        levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  recall_serostatus <- yardstick::recall_vec(estimate = factor(ensemblePrediction$max_col[idx],
                                                               levels = c("MPXV", "MVA", "Pre")),
                                             truth = factor(ensemblePrediction$real[idx],
                                                            levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  precision_serostatus <- yardstick::precision_vec(estimate = factor(ensemblePrediction$max_col[idx],
                                                                     levels = c("MPXV", "MVA", "Pre")),
                                                   truth = factor(ensemblePrediction$real[idx],
                                                                  levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  sens_serostatus <- yardstick::sens_vec(estimate = factor(ensemblePrediction$max_col[idx],
                                                           levels = c("MPXV", "MVA", "Pre")),
                                         truth = factor(ensemblePrediction$real[idx],
                                                        levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  spec_serostatus <- yardstick::spec_vec(estimate = factor(ensemblePrediction$max_col[idx],
                                                           levels = c("MPXV", "MVA", "Pre")),
                                         truth = factor(ensemblePrediction$real[idx],
                                                        levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  f1_model_serostatus[i] <-  f1_serostatus
  recall_model_serostatus[i] <- recall_serostatus
  precision_model_serostatus[i] <- precision_serostatus
  sens_model_serostatus[i] <- sens_serostatus
  spec_model_serostatus[i] <- spec_serostatus
  
  # 7: Mean
  f1_mean <- yardstick::f_meas_vec(estimate = factor(ensemblePrediction$mean_max_col[idx],
                                                     levels = c("MPXV", "MVA", "Pre")),
                                   truth = factor(ensemblePrediction$real[idx],
                                                  levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  recall_mean <- yardstick::recall_vec(estimate = factor(ensemblePrediction$mean_max_col[idx],
                                                         levels = c("MPXV", "MVA", "Pre")),
                                       truth = factor(ensemblePrediction$real[idx],
                                                      levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  precision_mean <- yardstick::precision_vec(estimate = factor(ensemblePrediction$mean_max_col[idx],
                                                               levels = c("MPXV", "MVA", "Pre")),
                                             truth = factor(ensemblePrediction$real[idx],
                                                            levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  sens_mean <- yardstick::sens_vec(estimate = factor(ensemblePrediction$mean_max_col[idx],
                                                     levels = c("MPXV", "MVA", "Pre")),
                                   truth = factor(ensemblePrediction$real[idx],
                                                  levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  spec_mean <- yardstick::spec_vec(estimate = factor(ensemblePrediction$mean_max_col[idx],
                                                     levels = c("MPXV", "MVA", "Pre")),
                                   truth = factor(ensemblePrediction$real[idx],
                                                  levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  f1_model_mean[i] <-  f1_mean
  recall_model_mean[i] <- recall_mean
  precision_model_mean[i] <- precision_mean
  sens_model_mean[i] <- sens_mean
  spec_model_mean[i] <- spec_mean
  
  # 8: LDA
  f1_LDA <- yardstick::f_meas_vec(estimate = factor(ensemblePrediction$LDA_max_col[idx],
                                                    levels = c("MPXV", "MVA", "Pre")),
                                  truth = factor(ensemblePrediction$real[idx],
                                                 levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  recall_LDA <- yardstick::recall_vec(estimate = factor(ensemblePrediction$LDA_max_col[idx],
                                                        levels = c("MPXV", "MVA", "Pre")),
                                      truth = factor(ensemblePrediction$real[idx],
                                                     levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  precision_LDA <- yardstick::precision_vec(estimate = factor(ensemblePrediction$LDA_max_col[idx],
                                                              levels = c("MPXV", "MVA", "Pre")),
                                            truth = factor(ensemblePrediction$real[idx],
                                                           levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  sens_LDA <- yardstick::sens_vec(estimate = factor(ensemblePrediction$LDA_max_col[idx],
                                                    levels = c("MPXV", "MVA", "Pre")),
                                  truth = factor(ensemblePrediction$real[idx],
                                                 levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  spec_LDA <- yardstick::spec_vec(estimate = factor(ensemblePrediction$LDA_max_col[idx],
                                                    levels = c("MPXV", "MVA", "Pre")),
                                  truth = factor(ensemblePrediction$real[idx],
                                                 levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  f1_model_LDA[i] <-  f1_LDA
  recall_model_LDA[i] <- recall_LDA
  precision_model_LDA[i] <- precision_LDA
  sens_model_LDA[i] <- sens_LDA
  spec_model_LDA[i] <- spec_LDA
  
  #9: GBC
  f1_GBC <- yardstick::f_meas_vec(estimate = factor(ensemblePrediction$GBC_max_col[idx],
                                                    levels = c("MPXV", "MVA", "Pre")),
                                  truth = factor(ensemblePrediction$real[idx],
                                                 levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  recall_GBC <- yardstick::recall_vec(estimate = factor(ensemblePrediction$GBC_max_col[idx],
                                                        levels = c("MPXV", "MVA", "Pre")),
                                      truth = factor(ensemblePrediction$real[idx],
                                                     levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  precision_GBC <- yardstick::precision_vec(estimate = factor(ensemblePrediction$GBC_max_col[idx],
                                                              levels = c("MPXV", "MVA", "Pre")),
                                            truth = factor(ensemblePrediction$real[idx],
                                                           levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  sens_GBC <- yardstick::sens_vec(estimate = factor(ensemblePrediction$GBC_max_col[idx],
                                                    levels = c("MPXV", "MVA", "Pre")),
                                  truth = factor(ensemblePrediction$real[idx],
                                                 levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  spec_GBC <- yardstick::spec_vec(estimate = factor(ensemblePrediction$GBC_max_col[idx],
                                                    levels = c("MPXV", "MVA", "Pre")),
                                  truth = factor(ensemblePrediction$real[idx],
                                                 levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  f1_model_GBC[i] <-  f1_GBC
  recall_model_GBC[i] <- recall_GBC
  precision_model_GBC[i] <- precision_GBC
  sens_model_GBC[i] <- sens_GBC
  spec_model_GBC[i] <- spec_GBC
  
  # 10: RF
  f1_RF <- yardstick::f_meas_vec(estimate = factor(ensemblePrediction$RF_max_col[idx],
                                                   levels = c("MPXV", "MVA", "Pre")),
                                 truth = factor(ensemblePrediction$real[idx],
                                                levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  recall_RF <- yardstick::recall_vec(estimate = factor(ensemblePrediction$RF_max_col[idx],
                                                       levels = c("MPXV", "MVA", "Pre")),
                                     truth = factor(ensemblePrediction$real[idx],
                                                    levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  precision_RF <- yardstick::precision_vec(estimate = factor(ensemblePrediction$RF_max_col[idx],
                                                             levels = c("MPXV", "MVA", "Pre")),
                                           truth = factor(ensemblePrediction$real[idx],
                                                          levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  sens_RF <- yardstick::sens_vec(estimate = factor(ensemblePrediction$RF_max_col[idx],
                                                   levels = c("MPXV", "MVA", "Pre")),
                                 truth = factor(ensemblePrediction$real[idx],
                                                levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
  spec_RF <- yardstick::spec_vec(estimate = factor(ensemblePrediction$RF_max_col[idx],
                                                   levels = c("MPXV", "MVA", "Pre")),
                                 truth = factor(ensemblePrediction$real[idx],
                                                levels = c("MPXV", "MVA", "Pre")), estimator = "macro")
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
f1_model_MPXV_GBC_Filtered <- numeric(n_boot)
f1_model_MVA_GBC_Filtered <- numeric(n_boot)
f1_model_Pre_GBC_Filtered <-  numeric(n_boot)

precision_model_MPXV_GBC_Filtered <- numeric(n_boot)
precision_model_MVA_GBC_Filtered <- numeric(n_boot)
precision_model_Pre_GBC_Filtered <-  numeric(n_boot)

recall_model_MPXV_GBC_Filtered <- numeric(n_boot)
recall_model_MVA_GBC_Filtered <- numeric(n_boot)
recall_model_Pre_GBC_Filtered <-  numeric(n_boot)

sens_model_MPXV_GBC_Filtered <- numeric(n_boot)
sens_model_MVA_GBC_Filtered <- numeric(n_boot)
sens_model_Pre_GBC_Filtered <-  numeric(n_boot)

spec_model_MPXV_GBC_Filtered <- numeric(n_boot)
spec_model_MVA_GBC_Filtered <-  numeric(n_boot)
spec_model_Pre_GBC_Filtered <-  numeric(n_boot)

f1_model_MPXV_GBC <- numeric(n_boot)
f1_model_MVA_GBC <- numeric(n_boot)
f1_model_Pre_GBC <-  numeric(n_boot)

precision_model_MPXV_GBC <- numeric(n_boot)
precision_model_MVA_GBC <- numeric(n_boot)
precision_model_Pre_GBC <-  numeric(n_boot)

recall_model_MPXV_GBC <- numeric(n_boot)
recall_model_MVA_GBC <- numeric(n_boot)
recall_model_Pre_GBC <-  numeric(n_boot)

sens_model_MPXV_GBC <- numeric(n_boot)
sens_model_MVA_GBC <- numeric(n_boot)
sens_model_Pre_GBC <- numeric(n_boot)

spec_model_MPXV_GBC <- numeric(n_boot)
spec_model_MVA_GBC <- numeric(n_boot)
spec_model_Pre_GBC <- numeric(n_boot)

for(i in 1:n_boot) {
  idx <- sample(1:n_val_filtered, size = n_val_filtered, replace = TRUE)
  confusionMatrixGBCFiltered <- 
    confusionMatrix(factor(ensembleCombinedFiltered$real)[idx], factor(ensembleCombinedFiltered$GBC_max_col)[idx])
  
  f1_MPXV_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "F1"][1]
  f1_MVA_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "F1"][2]
  f1_Pre_GBC_Filtered <-  confusionMatrixGBCFiltered[["byClass"]][, "F1"][3]
  
  precision_MPXV_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Precision"][1]
  precision_MVA_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Precision"][2]
  precision_Pre_GBC_Filtered <-  confusionMatrixGBCFiltered[["byClass"]][, "Precision"][3]
  
  recall_MPXV_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Recall"][1]
  recall_MVA_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Recall"][2]
  recall_Pre_GBC_Filtered <-  confusionMatrixGBCFiltered[["byClass"]][, "Recall"][3]
  
  sens_MPXV_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Sensitivity"][1]
  sens_MVA_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Sensitivity"][2]
  sens_Pre_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Sensitivity"][3]
  
  spec_MPXV_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Specificity"][1]
  spec_MVA_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Specificity"][2]
  spec_Pre_GBC_Filtered <- confusionMatrixGBCFiltered[["byClass"]][, "Specificity"][3]
  
  f1_model_MPXV_GBC_Filtered[i] <-  f1_MPXV_GBC_Filtered 
  f1_model_MVA_GBC_Filtered[i] <-   f1_MVA_GBC_Filtered
  f1_model_Pre_GBC_Filtered[i] <-  f1_Pre_GBC_Filtered
  
  precision_model_MPXV_GBC_Filtered[i] <- precision_MPXV_GBC_Filtered
  precision_model_MVA_GBC_Filtered[i] <- precision_MVA_GBC_Filtered
  precision_model_Pre_GBC_Filtered[i] <-  precision_Pre_GBC_Filtered 
  
  recall_model_MPXV_GBC_Filtered[i] <- recall_MPXV_GBC_Filtered 
  recall_model_MVA_GBC_Filtered[i] <- recall_MVA_GBC_Filtered
  recall_model_Pre_GBC_Filtered[i] <-  recall_Pre_GBC_Filtered 
  
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
    confusionMatrix(factor(ensemblePrediction$real)[idx], factor(ensemblePrediction$GBC_max_col)[idx])
  
  f1_MPXV_GBC <- confusionMatrixGBC[["byClass"]][, "F1"][1]
  f1_MVA_GBC <- confusionMatrixGBC[["byClass"]][, "F1"][2]
  f1_Pre_GBC <-  confusionMatrixGBC[["byClass"]][, "F1"][3]
  
  precision_MPXV_GBC <- confusionMatrixGBC[["byClass"]][, "Precision"][1]
  precision_MVA_GBC <- confusionMatrixGBC[["byClass"]][, "Precision"][2]
  precision_Pre_GBC <-  confusionMatrixGBC[["byClass"]][, "Precision"][3]
  
  recall_MPXV_GBC <- confusionMatrixGBC[["byClass"]][, "Recall"][1]
  recall_MVA_GBC <- confusionMatrixGBC[["byClass"]][, "Recall"][2]
  recall_Pre_GBC <-  confusionMatrixGBC[["byClass"]][, "Recall"][3]
  
  sens_MPXV_GBC <- confusionMatrixGBC[["byClass"]][, "Sensitivity"][1]
  sens_MVA_GBC <- confusionMatrixGBC[["byClass"]][, "Sensitivity"][2]
  sens_Pre_GBC <- confusionMatrixGBC[["byClass"]][, "Sensitivity"][3]
  
  spec_MPXV_GBC <- confusionMatrixGBC[["byClass"]][, "Specificity"][1]
  spec_MVA_GBC <- confusionMatrixGBC[["byClass"]][, "Specificity"][2]
  spec_Pre_GBC <- confusionMatrixGBC[["byClass"]][, "Specificity"][3]
  
  f1_model_MPXV_GBC[i] <-  f1_MPXV_GBC 
  f1_model_MVA_GBC[i] <-   f1_MVA_GBC
  f1_model_Pre_GBC[i] <-  f1_Pre_GBC
  
  precision_model_MPXV_GBC[i] <- precision_MPXV_GBC
  precision_model_MVA_GBC[i] <- precision_MVA_GBC
  precision_model_Pre_GBC[i] <-  precision_Pre_GBC 
  
  recall_model_MPXV_GBC[i] <- recall_MPXV_GBC 
  recall_model_MVA_GBC[i] <- recall_MVA_GBC
  recall_model_Pre_GBC[i] <-  recall_Pre_GBC 
  
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

table_performance_validation <- 
  f1_summary %>% 
  left_join(precision_summary, by = "Model", suffix = c("", "_prec")) %>% 
  left_join(recall_summary, by = "Model", suffix = c("", "_rec")) %>% 
  left_join(sens_summary, by = "Model", suffix = c("", "_rec")) %>% 
  left_join(spec_summary, by = "Model", suffix = c("", "_rec"))

save(table_performance_validation, file = "output/table_performance_validation_Rev_02.Rdata")

# Generate table with sensitivites and specificities by class for GBC
table_GBC_class_validation <- 
  sens_summary %>% 
  left_join(spec_summary, by = "Model", suffix = c("", "_rec")) %>% 
  filter(grepl("^MVA_|^MPXV|^Pre", Model)) %>% 
  select(Model, sensitivity = sens_formatted, specificity = spec_formatted)

save(table_GBC_class_validation, file = "output/table_GBC_class_validation_Rev_02.Rdata")




####
# Recode classification of real and prediction for MPXV differentiation for better comparison with binary classifiers
# 1) On All data fpr MPXV predictions
ensemblePredictionMPXVInf <-
  ensemblePredictionExport %>% 
  mutate(Real_Inf = if_else(Real == "MPXV", "Inf", "Non Inf"),
         Pred_GBC_Inf = if_else(Pred_GBC == "MPXV", "Inf", "Non Inf"),
         Pred_LDA_Inf = if_else(Pred_LDA == "MPXV", "Inf", "Non Inf"))

conf_mat_GBC_Inf <- confusion_matrix(targets = factor(ensemblePredictionMPXVInf$Real_Inf,
                                                      levels = c("Inf", "Non Inf"), ordered = TRUE),
                                     predictions = factor(ensemblePredictionMPXVInf$Pred_GBC_Inf,
                                                          levels = c("Inf", "Non Inf"), ordered = TRUE), metrics = list(all = TRUE), 
                                     positive = 1)
conf_mat_GBC_Inf_caret <-
  confusionMatrix(positive = "Inf",
                  factor(ensemblePredictionMPXVInf$Pred_GBC_Inf,
                         levels = c("Inf", "Non Inf"), ordered = TRUE),
                  factor(ensemblePredictionMPXVInf$Real_Inf,
                         levels = c("Inf", "Non Inf"), ordered = TRUE))


# 2) On filtered data for serostatus predictions
ensemblePredictionPos <-
  ensemblePrediction %>% 
  mutate(Real_Pos = if_else(real %in% c("MPXV", "MVA"), "positive", "negative"),
         Pred_GBC_Pos = if_else(GBC_max_col %in% c("MPXV", "MVA"), "positive", "negative"),
         Pred_LDA_Pos = if_else(LDA_max_col %in% c("MPXV", "MVA"), "positive", "negative")) %>% 
  filter(!(real == "Pre" & year_birth < 1972)) %>% 
  filter(!(real == "Pre" & highBG == TRUE))

conf_mat_GBC_Pos <- confusion_matrix(targets = factor(ensemblePredictionPos$Real_Pos),
                                     predictions = factor(ensemblePredictionPos$Pred_GBC_Pos))

conf_mat_GBC_Pos_caret <-
  confusionMatrix(factor(ensemblePredictionPos$Pred_GBC_Pos,
                         levels = c("positive", "negative"), ordered = TRUE),
                  factor(ensemblePredictionPos$Real_Pos,
                         levels = c("positive", "negative"), ordered = TRUE))


conf_mat_LDA_Pos_caret <-
  confusionMatrix(factor(ensemblePredictionPos$Pred_LDA_Pos,
                         levels = c("positive", "negative"), ordered = TRUE),
                  factor(ensemblePredictionPos$Real_Pos,
                         levels = c("positive", "negative"), ordered = TRUE))

conf_mat_LDA_Pos <- confusion_matrix(targets = factor(ensemblePredictionPos$Real_Pos),
                                     predictions = factor(ensemblePredictionPos$Pred_LDA_Pos))



####
# Generate plots for filtered data
plotValidation <-
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
  ggtitle("Ensemble serostatus confident") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotValidationLDA <-
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
  # ggtitle("LDA") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotValidationGBC <-
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
  #  ggtitle("GBC") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotValidationRF <-
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
  # ggtitle("RF") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



plotPreStratified <-
  ggarrange(plotPreGroupedLDA, plotPreGroupedGBC, plotPreGroupedRF,
            align = "hv", ncol = 3, nrow = 1)

ggsave("output/Fig_S10.png",plotPreStratified ,
       width = 10, height = 4, dpi = 600)



####
# Generate plots for all data
plotValidationAll <-
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


plotValidationLDAAll <-
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


plotValidationGBCAll <-
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

# Rev 02
plotValidationRFAll <-
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



####
# Plot ensemble all

plotEnsembleAll <-
  ggarrange(plotValidationLDAAll, plotValidationGBCAll, plotValidationRFAll, plotValidationAll, ncol = 4, 
            align = "hv", labels = c("b"))

ggsave("output/plotEnsembleAll.pdf", plotEnsembleAll,
       width = 10, height = 4, dpi = 600)

