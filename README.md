# Mpox_Multiplex_Assay
# Description of Data and Analysis Files

The following data and analysis files are contained in the subfolders accompanying the manuscript 
"Differentiation between mpox infection and MVA immunization by a novel machine learning-supported 
serological multiplex assay" by Surtees et al.

## Installation and usage
Install R (version 4.3.0) and R Studio (2023-06-2). Install all necessary packages, which are listed in each R script before running the scripts. All necessary data files (input and output) are contained within the different folders.


## Figure_02_Compare_Panels
Plot IgG and IgM results of the differen serogroups (pre-immune, MVA, MPXV) and panels (acute, epi) stratified by childhood immunisation against smallpox. Analysis and plotting is performed in the R script `Figure_2_Compare_Panels_Rev_02_V2.R`.

Depends on the following input files:
- `input/dataInputComparePanels.Rdata`

The following output files are generated:
- Figure 2: Spider plots and plots of selected antigens from the acute and epi panel. Plot of ratios of selected homologue antigen pairs.
- Figure S3: Plot of antigens not contained in figure 2
- Figure S4: Plot of IgM results
- Figure S5: Plot of ratios not contained in Figure 2


## Figure_03_ML_Performance
Reads assay performance of the comparison of different ML algorithms tested on different panels. Analysis and plotting is performed in the R script `Figure_3_Rev_02.R`.

Depends on the following input files:
- `input/dataInMeta.Rdata`
- `input/statisticalDataCombined.Rdata`

The following output files are generated:
- Fig3a.pdf: Plot of F1 scores for the comparison of the different algorithms and panels. 
- Fig3_all_xgboost_lda.pdf: Circular plot of misclassifications on all data set, stratified by childhood smallpox vaccionation for GBC and LDA
- Fig3_all_rf.pdf: Circular plot of misclassifications on all data set, stratified by childhood smallpox vaccionation for RF
- freqTableRFGBCLDA.xlsx: Frequency table underlying the circular plots in figure 3
- freqTableAll.Rdata: Frequence table for misclassificaton for other algorithms
- supTablePerformanceRev.Rdata: Data underlying Table S9 (Generated using SupportingTablesPerformanceRev.Rmd -> Generates Docx file)


## Figure_03d_S08_Conf_Matr
Generates confusion matrices for ensemble predictions on the combined acute and epi panels. Plotting is performed in the R script  `Figure_S08_Conf_Mat_Rev_02_V5.R`.

Depends on the following input files:
- `input/ensembleCombined_Rev_02.Rdata`

The following output files are generated:
- plotConfusionTest_Rev_02.png: Plot of confusion matrices for only confident predictions (> 0.5)
- plotConfusionTest_Rev_02.pdf: Plot of confusion matrices for only confident predictions (> 0.5)
- plotConfusionTestAll_Rev_02.png: Plot of confusion matrices for all predictions
- plotConfusionTestAll_Rev_02.pdf: Plot of confusion matrices for all predictions
- table_GBC_class_establishment_Rev_04.Rdata: Performance of bootstrap analysis for GBC only, class-wise stratified. Used as input to generate table S14
- table_performance_establishment_Rev_04.Rdata:Performance of bootstrap analysis for algorithms applied to the combined acute and epi panel. Used as input to generate table S13


## Figure_04_Val_V2
Analysis and generation of figures and table for the validation measurements. Analysis and plotting is performed in the R script `Figure_4_Val_Rev_02_V6.R`.

Depends on the following input files:
- `input/dataInputComparison.Rdata`
- `input/ensemblePrediction_Rev_02.Rdata`
- `input/heatmap_input.Rdata`

The following output files are generated:
- heatmap_IgG_white.pdf: Heatmap of validation panel together with plots of metadata for IgG data
- heatmap_IgM_white.pdf: Heatmap of validataion panel together with plots of metadata for IgM data
- plotEnsembleAll.pdf: Confusion matrices for validation measurements depending on all data
- STableEnsemblePrediction.xlsx: Table S10 with all data for the validation panel
- table_GBC_class_validation_Rev_02.Rdata: Performance of bootstrap analysis for GBC only, class-wise stratified. Used as input to generate table S14
- table_performance_validation_Rev_02.Rdata: Performance of bootstrap analysis for algorithms applied to the combined acute and epi panel. Used as input to generate table S13



## Figure_S01_Method_Comparison
Compare results of novel multiplex assay with reference methods (ELISA, IFA, and NT) for a panel of selected sera, which have been quantified by the different methods. Method comparison is performed in the R script `Fig_1_Method_Comparison.R`.

Depends on the following input files:
- `input/dataInputELISAMultiplex.Rdata`
- `input/dataInputNT.Rdata`
- `input/dataInputQuant.Rdata`
- `input/metadata_IFA.Rdata`

The following output files are generated:
- Figure 1: Comparison between results obtained by the multiplex assay and ELISA, IFA, and NT as reference methods. Correlation plots and plot for comparison with Delta-VACV antigen are combined in Adobe Illustrator
- Figure S1: Plot of multiplex antigens in comparison to ELISA as reference method
- parampaBa.xlsx: Parameters for Passing-Bablock-Regression save as Excel-file


## Figure_S06_Bead_coupling
Generates plot for coupling control and batch-to-batch variabiliy using in the R script analyseCoupling.R.

Depends on the following input files:
- `input/dataInputBatch.Rdata`
- `input/dataInputPlotting.Rdata`

The following output files are generated:
- plotCoupling.png: Figure S4 containing the plots for the coupling control and the batch-to-batch variability


## Figure_S06_07_Exclusion_Training
Generates plots and population-based cut-off values used to exclude samples with hints on unrecognized infections from the epi panel from ML based training and testing. Not to confuse with more robust cut-off values, which have been generated to classify sera based on binary classifiers as described in Figure_S07_ROC.

Depends on the following input files: 
- `input/dataClustering_post.Rdata`
- `input/dataInWide.Rdata`

The following output files are generated (amoung others):
- SFig_plotCombinedyoung.png: Figure S5 Comparison between antigens and panels in young population
- SFig_plotROC.png: Figure S6 ROC curve


## Figure_S08_Classical_ROC
Generate plot for the performance of single antigens based on ROC analysis as well as threshold values. 

Depends on the following input files: 
- `input/dataInWide.Rdata`
- `input/plotImprovePerformance.Rdata`

The following output files are generated:
- plotFig6.png: Figure S8
- ROC_parameters_classic.xlsx: Table S11 with performance parameters for ROC analysis for single antigens
