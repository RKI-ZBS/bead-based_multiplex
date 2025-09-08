# Mpox Multiplex Assay

**Differentiation between mpox infection and MVA immunization by a novel machine learning–supported serological multiplex assay**
*Surtees et al.*

[![R Version](https://img.shields.io/badge/R-v4.3.0-blue.svg)](https://cran.r-project.org/)
[![RStudio](https://img.shields.io/badge/RStudio-2025.05.1-blue.svg)](https://posit.co/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/RKI-ZBS/bead-based_multiplex)](https://github.com/RKI-ZBS/bead-based_multiplex/commits/main)



---

## **Table of Contents**

* [Overview](#overview)
* [Installation & Usage](#installation--usage)
* [Analysis & Figures](#analysis--figures)
* [Tables](#tables)
* [Changelog](#changelog)
* [License & Citation](#license--citation)

---

## **Overview**

This repository contains all data, R scripts, and analysis files to recreate all figures and tables from the manuscript:

> **“Differentiation between mpox infection and MVA immunization by a novel machine learning–supported serological multiplex assay”**
> by Surtees et al.

It provides:

* R scripts for analysis, figure, and table generation
* Input and output data files
* All supplementary tables and figures

---

## **Installation & Usage**

### **1. Prerequisites**

* **R** ≥ 4.3.0
* **RStudio** ≥ 2023.06.2
* Operating systems: Linux, macOS, Windows

### **2. Install Required R Packages**

Each R script specifies required packages at the top.
Install them using:

```R
install.packages(c("tidyverse", "ggplot2", "ComplexHeatmap", "caret", "yardstick"))
```

*(Replace with packages listed in the scripts.)*

### **3. Run Analyses**

* Navigate to the corresponding figure/table folder.
* Open the respective R script in RStudio.
* Run the script to generate figures/tables.
* All necessary  input files are contained within the input folder

---

## **Analysis & Figures**

| Folder | Figure(s) / Table(s) | Purpose                                      | R Script                               | Output                                                  |
| ------ | -------------- | -------------------------------------------- | ------------------------------------- | ------------------------------------------------------- |
| Figure_2_Compare_Panels | **Figure 2**, **Figure S3**, **Figure S4**, **Figure S5**   | Compare IgG & IgM across serogroups & panels | `Figure_2_Compare_Panels_Rev_02_V2.R` | Spider plots, antigen ratios, IgG and IgM plots                 |
| Figure_3ab_ML_Performance | **Figure 3**, **Table S9**   | ML performance comparison                    | `Figure_3_Rev_02.R`                   | F1 plots, circular misclassification plots, supporting table performance comparison |
| Figure_3c_Conf_Matr | **Figure 3c**  | Confusion matrices for LDA, GBC, RF                    | `Figure_03c_Rev_02.R`                   | Plot confuction matrices, tables of ensemble predicition on combined cohort |
| Figure_4_Val_V2 | **Figure 4**, **Fig_S10**, **Table_S10**,   | Validation panel analysis                    | `Figure_4_Val_Rev_02_V6.R`            | IgG/IgM heatmaps, ensemble confusion matrices validation, tables of ensemble prediction on validation cohort           |
| Figure_S1_Method_Comparison | **Figure S1**, **Figure S2**  | Comparison with ELISA/IFA/NT                 | `Fig_1_Method_Comparison.R`           | Correlation plots, Passing-Bablok regression            |
|  | **Figure S6**  | Bead coupling quality                        | `analyseCoupling.R`                   | Coupling control and variability plots                  |
| Figure_S6_Bead_coupling | **Figure S6**  | Antigen specificity and coupling efficiency                          | `Figure_S06.R`     | Coupling control plots                          |
| Figure_S7_8 | **Figure S7**, **Figure S8**  | Determine population based cutoffs in defined populations, exclude young positive samples with pre-immune status in epi cohort                         | `Figures_S7_S8_Cutoff_Exclusion.R`     | Plot with reactivity in different cohorts and ROC output with threshold                          |
| Figure_S9_Rec_Feat_Select | **Figure S9**  | Feature elimination                          | `Figure_S9.R`              | F1 impact plots after antigen removal                                        |
| Figure_S11_Classical_ROC | **Figure S11**, **Table S12** | ROC performance                              | `Figure_S11_Classical_ROC_Rev_02.R`   | ROC curves, threshold parameters                        |
| Figure_S12_Density_Plots | **Figure S12** | Density plots                                | `Fig_S12.R`                       | Density plots                                           |
| Figure_S13_Repro | **Figure S13** | Reproducibility                              | `plotComparison_V3.R`                 | Reproducibility figure                                  |

*(Full details are in each figure’s folder.)*

---

## **Tables**

| Folder | Table                | Purpose                                 | Output                                                   |
| ------ | -------------------- | --------------------------------------- |-------------------------------------------------------- |
| Table_3_Ensemble_Performance | **Table 3**, **Table S11**, **Table S12**          | Ensemble performance bootstrap analysis | `table_3.xlsx`, `table_s11.xlsx`, `table_s12.xlsx`       |                  |
| Table_S13_S14 | **Tables S13 & S14** | Single vs ML classifier comparison      | `supporting_table_s13.xlsx`, `supporting_table_s14.xlsx` |

---

## **Changelog**

* Unified naming conventions for figures & tables.
* Added structured tables summarizing inputs & outputs.
* Improved Markdown formatting for GitHub readability.
* Harmonized names and content of folders 

---

## **License & Citation**

This repository is distributed under the [MIT License](LICENSE).
If you use this code or data, please cite:

> Surtees et al. *Differentiation between mpox infection and MVA immunization by a novel machine learning–supported serological multiplex assay*, **\[Under Revision]**
