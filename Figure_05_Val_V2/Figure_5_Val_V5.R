####
# Analyse revision measurements
# Daniel Stern
# 2025-03-31
# Version 2.0
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

# Load datainput
load("input/ensemblePrediction.Rdata")
load("input/heatmap_input.Rdata")
load("input/dataInputComparison.Rdata")

heatmap_input_IgG <- 
  heatmap_input_IgG %>% 
  mutate(panel = factor(panel, levels = c("Pre", "MVA", "MPXV"), 
                        ordered = TRUE))

plotPreGroupedLDA <-
  ensemblePrediction %>%
  filter(pred_ensemble > 0.5) %>% 
  filter(mean_mean_conf > 0.5) %>% 
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
  filter(pred_ensemble > 0.5) %>% 
  filter(mean_mean_conf > 0.5) %>% 
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


####
# Plot heatmap of raw data

col_fun_MVA <- colorRamp2(c(0.5, 1), c("white", colorblind_pal()(8)[3]))
col_fun_Pre <- colorRamp2(c(0.5, 1), c("white", colorblind_pal()(8)[2]))
col_fun_MPXV <- colorRamp2(c(0.5, 1), c("white", colorblind_pal()(8)[4]))

set.seed(123)
matIgG <- as.matrix(scale(heatmap_input_IgG[2:16]))
haIgG <- HeatmapAnnotation(df = data.frame(panel = heatmap_input_IgG$panel),
                           annotation_height = unit(4, "mm"))

pdf(file = "output/heatmap_IgG_white.pdf", width = 14, height = 10)
Heatmap(matIgG, split = heatmap_input_IgG$panel,
        name = "Multiplex quantfied", clustering_method_columns = "ward.D2", 
        clustering_method_rows = "ward.D2") +
  Heatmap(factor(heatmap_input_IgG$serostatus.delta, levels = c(0,1),
                 labels = c("negative", "positive"),
                 ordered = TRUE), name = "Delta-VACV", col = c(colorblind_pal()(8)[c(6,7)])) +
  
  Heatmap(heatmap_input_IgG$n_MVA_vacc, name = "MVA vaccinations", na_col = "white", col = c("grey", viridis_pal()(5)[2:5])) +
  Heatmap(heatmap_input_IgG$LDA_MVA, name = "MVA LDA", col = col_fun_MVA, na_col = "white") +
  Heatmap(heatmap_input_IgG$GBC_MVA, name = "MVA GBC", col = col_fun_MVA, na_col = "white") +
  #  Heatmap(heatmap_input_IgG$MVA, name = "MVA ensemble mean", col = col_fun_MVA, na_col = "white") +
  #  Heatmap(heatmap_input_IgG$Strat_MVA, name = "MVA ensemble serostatus", col = col_fun_MVA, na_col = "white") +
  Heatmap(heatmap_input_IgG$LDA_Pre, name = "Pre LDA", col = col_fun_Pre, na_col = "white") +
  Heatmap(heatmap_input_IgG$GBC_Pre, name = "Pre GBC", col = col_fun_Pre, na_col = "white") +
  #  Heatmap(heatmap_input_IgG$Pre, name = "Pre ensemble mean", col = col_fun_Pre, na_col = "white") +
  #  Heatmap(heatmap_input_IgG$Strat_Pre, name = "Pre ensemble serostatus", col = col_fun_Pre, na_col = "white") +
  Heatmap(heatmap_input_IgG$LDA_MPXV, name = "MPXV LDA", col = col_fun_MPXV, na_col = "white") +
  Heatmap(heatmap_input_IgG$GBC_MPXV, name = "MPXV GBC", col = col_fun_MPXV, na_col = "white") +
  #  Heatmap(heatmap_input_IgG$MPXV, name = "MPXV ensemble mean", col = col_fun_MPXV, na_col = "white") +
  #  Heatmap(heatmap_input_IgG$Strat_MPXV, name = "MPXV ensemble serostatus", col = col_fun_MPXV, na_col = "white") +
  Heatmap(factor(heatmap_input_IgG$panel, levels = c("Pre", "MVA", "MPXV"), 
                 ordered = TRUE), name = "Panel", col = c(colorblind_pal()(8)[2:4]))
dev.off()


matIgM <- as.matrix(scale(heatmap_input_IgM[2:16]))
ha <- HeatmapAnnotation(df = data.frame(panel = heatmap_input_IgM$panel),
                        annotation_height = unit(4, "mm"))
pdf(file = "output/heatmap_IgM_white.pdf", width = 14, height = 8)
Heatmap(matIgM, split = heatmap_input_IgM$panel,
        name = "Multiplex quantfied", clustering_method_columns = "ward.D2", 
        clustering_method_rows = "ward.D2") +
  Heatmap(factor(heatmap_input_IgM$serostatus.delta, levels = c(0,1),
                 labels = c("negative", "positive"),
                 ordered = TRUE), name = "Delta-VACV", col = c(colorblind_pal()(8)[c(6,7)])) +
  Heatmap(heatmap_input_IgM$n_MVA_vacc, name = "MVA vaccinations", na_col = "white", col = c("grey", viridis_pal()(5)[2:5])) +
  
  Heatmap(heatmap_input_IgM$LDA_MVA, name = "MVA LDA", col = c("white", colorblind_pal()(8)[3]), na_col = "white") +
  Heatmap(heatmap_input_IgM$GBC_MVA, name = "MVA GBC", col = c("white", colorblind_pal()(8)[3]), na_col = "white") +
  Heatmap(heatmap_input_IgM$MVA, name = "MVA ensemble mean", col = c("white", colorblind_pal()(8)[3]), na_col = "white") +
  Heatmap(heatmap_input_IgM$Strat_MVA, name = "MVA ensemble serostatus", col = c("white", colorblind_pal()(8)[3]), na_col = "white") +
  Heatmap(heatmap_input_IgM$LDA_Pre, name = "Pre LDA", col = c("white", colorblind_pal()(8)[2]), na_col = "white") +
  Heatmap(heatmap_input_IgM$GBC_Pre, name = "Pre GBC", col = c("white", colorblind_pal()(8)[2]), na_col = "white") +
  Heatmap(heatmap_input_IgM$Pre, name = "Pre ensemble mean", col = c("white", colorblind_pal()(8)[2]), na_col = "white") +
  Heatmap(heatmap_input_IgM$Strat_Pre, name = "Pre ensemble serostatus", col = c("white", colorblind_pal()(8)[2]), na_col = "white") +
  Heatmap(heatmap_input_IgM$LDA_MPXV, name = "MPXV LDA", col = c("white", colorblind_pal()(8)[4]), na_col = "white") +
  Heatmap(heatmap_input_IgM$GBC_MPXV, name = "MPXV GBC", col = c("white", colorblind_pal()(8)[4]), na_col = "white") +
  Heatmap(heatmap_input_IgM$MPXV, name = "MPXV ensemble mean", col = c("white", colorblind_pal()(8)[4]), na_col = "white") +
  Heatmap(heatmap_input_IgM$Strat_MPXV, name = "MPXV ensemble serostatus", col = c("white", colorblind_pal()(8)[4]), na_col = "white") +
  Heatmap(factor(heatmap_input_IgM$panel, levels = c("Pre", "MVA", "MPXV"), 
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
         Pred_GBC = GBC_max_col,
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
# Generate table with assay performance parameters
table_performance_validatation <- 
  add_column(as_tibble(conf_mat[6:8]), 
             n = sum(conf_mat$Table[[1]]), 
             panel = "Validation",
             algorithm = "Ensemble serostatus", 
             filtered = "confidence > 0.5") %>% 
  select(panel, algorithm, filtered, n, F1, Sensitivity, Specificity) %>% 
  add_row(add_column(as_tibble(conf_mat_mean[6:8]), 
                     n = sum(conf_mat_mean$Table[[1]]), 
                     panel = "Validation",
                     algorithm = "Ensemble mean", 
                     filtered = "confidence > 0.5")) %>% 
  add_row(add_column(as_tibble(conf_mat_LDA[6:8]), 
                     n = sum(conf_mat_LDA$Table[[1]]), 
                     panel = "Validation",
                     algorithm = "LDA", 
                     filtered = "confidence > 0.5")) %>% 
  add_row(add_column(as_tibble(conf_mat_GBC[6:8]), 
                     n = sum(conf_mat_GBC$Table[[1]]), 
                     panel = "Validation",
                     algorithm = "GBC", 
                     filtered = "confidence > 0.5")) %>% 
  add_row(add_column(as_tibble(conf_mat_all[6:8]), 
                     n = sum(conf_mat_all$Table[[1]]), 
                     panel = "Validation",
                     algorithm = "Ensemble serostatus", 
                     filtered = "none")) %>% 
  add_row(add_column(as_tibble(conf_mat_mean_all[6:8]), 
                     n = sum(conf_mat_mean_all$Table[[1]]), 
                     panel = "Validation",
                     algorithm = "Ensemble mean", 
                     filtered = "none")) %>% 
  add_row(add_column(as_tibble(conf_mat_LDA_all[6:8]), 
                     n = sum(conf_mat_LDA_all$Table[[1]]), 
                     panel = "Validation",
                     algorithm = "LDA", 
                     filtered = "none")) %>% 
  add_row(add_column(as_tibble(conf_mat_GBC_all[6:8]), 
                     n = sum(conf_mat_GBC_all$Table[[1]]), 
                     panel = "Validation",
                     algorithm = "GBC", 
                     filtered = "none")) 

save(table_performance_validatation, file = "output/table_performance_validatation.Rdata")

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

plotValidationMean <-
  plot_confusion_matrix(conf_mat_mean,
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
  ggtitle("Ensemble mean") +
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
  ggtitle("LDA") +
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
  ggtitle("GBC") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotConfmatrices <-
  ggarrange(plotPreGroupedLDA, plotPreGroupedGBC, plotValidationLDA, plotValidationGBC,
            align = "hv", ncol = 4, labels = c("b", "", "c", "", ""), legend = F)


ggsave("output/plotConfusionVal.pdf", plotConfmatrices,
       width = 12, height = 4, dpi = 600)

plotConfmatricesHigh <-
  ggarrange(ggarrange(plotPreGroupedLDA, plotPreGroupedGBC, ncol = 2, align = "hv"), plotValidationLDA, plotValidationGBC,
            align = "hv", nrow = 3, heights = c(2,1,1), labels = c("b", "c", "c"), legend = F)


ggsave("output/plotConfusionValHigh.pdf", plotConfmatricesHigh,
       width = 4, height = 8, dpi = 600)


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
  ggtitle("Ensemble serostatus all") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotValidationMeanAll <-
  plot_confusion_matrix(conf_mat_mean_all,
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
  ggtitle("Ensemble mean") +
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


plotConfmatricesValAll <-
  ggarrange(plotValidationLDAAll, plotValidationGBCAll,
            align = "hv", ncol = 2)


ggsave("output/plotConfusionValAll.pdf", plotConfmatricesValAll,
       width = 8, height = 4, dpi = 600)


ggsave("output/plotConfusionValAll.png", plotConfmatricesValAll,
       width = 8, height = 4, dpi = 600)


####
# Plot validation of ensemble approach
plotValidationEnsemble <-
  ggarrange(plotValidation, plotValidationAll, ncol = 2, align = "hv")
ggsave("output/plotConfusionValEnsemble.png", plotValidationEnsemble,
       width = 8, height = 4, dpi = 600)
