####
# Analyse stratified ensemble learning on test dataset
# Daniel Stern
# 2025-03-24
# Version 1.0
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

# Load datainput
load("input/ensembleCombined.Rdata")



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
confusionMatrixEnsembleFiltered[["byClass"]][ , "F1"]

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
  #  filter(mean_mean_conf > 0.5) %>% 
  filter(GBC_pred_ensemble > 0.5)

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


ensembleCombined <-
  ensembleCombined %>% 
  filter(!is.na(max_col))



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

sum(conf_mat_all$Table[[1]])
sum(conf_mat$Table[[1]])
sum(conf_mat_LDA$Table[[1]])
sum(conf_mat_GBC$Table[[1]])
sum(conf_mat_mean$Table[[1]])

####
# Generate output table for joining with results from validation

sum(conf_mat$Table[[1]])

table_performance_establishment <- 
  add_column(as_tibble(conf_mat[6:8]), 
             n = sum(conf_mat$Table[[1]]), 
             panel = "Establishment",
             algorithm = "Ensemble serostatus", 
             filtered = "confidence > 0.5") %>% 
  select(panel, algorithm, filtered, n, F1, Sensitivity, Specificity) %>% 
  add_row(add_column(as_tibble(conf_mat_mean[6:8]), 
                     n = sum(conf_mat_mean$Table[[1]]), 
                     panel = "Establishment",
                     algorithm = "Ensemble mean", 
                     filtered = "confidence > 0.5")) %>% 
  add_row(add_column(as_tibble(conf_mat_LDA[6:8]), 
                     n = sum(conf_mat_LDA$Table[[1]]), 
                     panel = "Establishment",
                     algorithm = "LDA", 
                     filtered = "confidence > 0.5")) %>% 
  add_row(add_column(as_tibble(conf_mat_GBC[6:8]), 
                     n = sum(conf_mat_GBC$Table[[1]]), 
                     panel = "Establishment",
                     algorithm = "GBC", 
                     filtered = "confidence > 0.5")) %>% 
  add_row(add_column(as_tibble(conf_mat_all[6:8]), 
                     n = sum(conf_mat_all$Table[[1]]), 
                     panel = "Establishment",
                     algorithm = "Ensemble serostatus", 
                     filtered = "none")) %>% 
  add_row(add_column(as_tibble(conf_mat_mean_all[6:8]), 
                     n = sum(conf_mat_mean_all$Table[[1]]), 
                     panel = "Establishment",
                     algorithm = "Ensemble mean", 
                     filtered = "none")) %>% 
  add_row(add_column(as_tibble(conf_mat_LDA_all[6:8]), 
                     n = sum(conf_mat_LDA_all$Table[[1]]), 
                     panel = "Establishment",
                     algorithm = "LDA", 
                     filtered = "none")) %>% 
  add_row(add_column(as_tibble(conf_mat_GBC_all[6:8]), 
                     n = sum(conf_mat_GBC_all$Table[[1]]), 
                     panel = "Establishment",
                     algorithm = "GBC", 
                     filtered = "none")) 

save(table_performance_establishment, file = "output/table_performance_establishment.Rdata")

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


plotConfMatrixTest <-
  ggarrange(plotTestLDA, plotTestGBC, ncol = 2, 
            align = "hv")

ggsave("output/plotConfusionTest.pdf", plotConfMatrixTest,
       width = 8, height = 4, dpi = 600)

ggsave("output/plotConfusionTest.png", plotConfMatrixTest,
       width = 8 , height = 4, dpi = 600)


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


plotConfMatrixTestAll <-
  ggarrange(plotTestLDAAll, plotTestGBCAll, ncol = 2,
            align = "hv")

ggsave("output/plotConfusionTestAll.pdf", plotConfMatrixTestAll,
       width = 8, height = 4, dpi = 600)

ggsave("output/plotConfusionTestAll.png", plotConfMatrixTestAll,
       width = 8, height = 4, dpi = 600)
