# Mpox Multiplex Assay

**Differentiation between mpox infection and MVA immunization by a novel machine learning–supported serological multiplex assay**
*Surtees et al.*

[![R Version](https://img.shields.io/badge/R-v4.3.0-blue.svg)](https://cran.r-project.org/)
[![RStudio](https://img.shields.io/badge/RStudio-2023.06.2-blue.svg)](https://posit.co/)
[![Last Commit](https://img.shields.io/github/last-commit/RKI-ZBS/bead-based_multiplex/)](https://github.com/RKI-ZBS/bead-based_multiplex/commits/main)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

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

This repository contains all data, R scripts, and analysis files supporting the manuscript:

> **“Differentiation between mpox infection and MVA immunization by a novel machine learning–supported serological multiplex assay”**
> by Surtees et al.

It provides:

* R scripts for statistical analysis and figure generation
* Raw and processed data files
* Machine learning model performance comparisons
* Supplementary tables and reproducibility assessments

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
* Ensure the required input `.Rdata` files are present.
* Run the script to generate figures/tables.

---

## **Analysis & Figures**

| Figure / Table | Purpose                                      | R Script                              | Input Files                                                                           | Output                                                  |
| -------------- | -------------------------------------------- | ------------------------------------- | ------------------------------------------------------------------------------------- | ------------------------------------------------------- |
| **Figure 2**   | Compare IgG & IgM across serogroups & panels | `Figure_2_Compare_Panels_Rev_02_V2.R` | `dataInputComparePanels.Rdata`                                                        | Spider plots, antigen ratios, IgM plots                 |
| **Figure 3**   | ML performance comparison                    | `Figure_3_Rev_02.R`                   | `dataInMeta.Rdata`, `statisticalDataCombined.Rdata`                                   | F1 plots, circular misclassification plots, freq tables |
| **Figure S8**  | Confusion matrices                           | `Figure_S08_Conf_Mat_Rev_02_V5.R`     | `ensembleCombined_Rev_02.Rdata`                                                       | Confusion matrix PDFs & PNGs                            |
| **Figure 4**   | Validation panel analysis                    | `Figure_4_Val_Rev_02_V6.R`            | `dataInputComparison.Rdata`, `ensemblePrediction_Rev_02.Rdata`, `heatmap_input.Rdata` | IgG/IgM heatmaps, ensemble confusion matrices           |
| **Figure S1**  | Comparison with ELISA/IFA/NT                 | `Fig_1_Method_Comparison.R`           | Multiple `.Rdata` inputs                                                              | Correlation plots, Passing-Bablok regression            |
| **Figure S6**  | Bead coupling quality                        | `analyseCoupling.R`                   | `dataInputBatch.Rdata`, `dataInputPlotting.Rdata`                                     | Coupling control and variability plots                  |
| **Figure S9**  | Feature elimination                          | See script                            | `dataInputFeatElim.Rdata`                                                             | F1 impact plots                                         |
| **Figure S11** | ROC performance                              | See script                            | `dataInWide.Rdata`                                                                    | ROC curves, threshold parameters                        |
| **Figure S12** | Density plots                                | See script                            | `dataInputComparePanels.Rdata`, `heatmap_input.Rdata`                                 | Density plots                                           |
| **Figure S13** | Reproducibility                              | See script                            | `dataInRep.Rdata`, `dataInRepSpoxFiltered.Rdata`                                      | Reproducibility figure                                  |

*(Full details are in each figure’s folder.)*

---

## **Tables**

| Table                | Purpose                                 | Input Files                       | Output                                                   |
| -------------------- | --------------------------------------- | --------------------------------- | -------------------------------------------------------- |
| **Table 3**          | Ensemble performance bootstrap analysis | Multiple `.Rdata`                 | `table_3.xlsx`, `table_s11.xlsx`, `table_s12.xlsx`       |
| **Table S10**        | Validation panel results                | `ensemblePrediction_Rev_02.Rdata` | `STableEnsemblePrediction.xlsx`                          |
| **Table S12**        | ROC performance parameters              | `dataInWide.Rdata`                | `ROC_parameters_classic.xlsx`                            |
| **Tables S13 & S14** | Single vs ML classifier comparison      | Multiple `.Rdata`                 | `supporting_table_s13.xlsx`, `supporting_table_s14.xlsx` |

---

## **Changelog**

* Unified naming conventions for figures & tables.
* Added structured tables summarizing inputs & outputs.
* Improved Markdown formatting for GitHub readability.
* Added badges for R version, RStudio, Zenodo DOI, and repo status.

---

## **License & Citation**

This repository is distributed under the [MIT License](LICENSE).
If you use this code or data, please cite:

> Surtees et al. *Differentiation between mpox infection and MVA immunization by a novel machine learning–supported serological multiplex assay*, **\[Journal Name, Year, DOI]**
