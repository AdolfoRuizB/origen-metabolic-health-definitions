# Metabolic health definitions in Mexican adults: oriGen analysis code

R code for all analyses in:

> Ruiz-Ballesteros AI, et al. Metabolic health definitions and cardiometabolic risk stratification in Mexican adults: a cross-sectional analysis of the oriGen project. *The Lancet Regional Health – Americas* [year; volume: pages]. DOI: [pending]

## Overview

This repository contains the scripts used to compare two operational definitions of metabolic health in adults from the oriGen cohort (Mexico):

- **Definition A**: harmonised metabolic syndrome criteria with Latin American waist-circumference cut-offs; metabolically unhealthy (MUH) = 2 or more of 5 components.
- **Definition B**: empirically derived definition of Zembic and colleagues; MUH = 1 or more of 3 components (self-reported diabetes, systolic blood pressure ≥130 mm Hg or antihypertensive treatment, and waist-to-hip ratio ≥0.95 in women or ≥1.03 in men).

The analyses cover the characteristics of the study population, agreement between definitions, the contribution of each component, discrimination of prevalent cardiovascular disease (CVD), fatty liver disease (FLD), and chronic kidney disease (CKD) by both definitions and by 27 cardiometabolic indices, and risk gradients of these indices within metabolic phenotypes.

## Data

Individual participant data are **not included** in this repository. The data analysed were obtained from the oriGen Project (data release of July, 2026) and are not publicly available. De-identified data may be made available upon specific request to the oriGen Project, subject to approval of a research proposal in accordance with the oriGen data sharing policy:

- Requests: origen.research@servicios.tecsalud.mx
- Data access forms: https://tec.mx/en/research/origen-project/Researchers

## Scripts

Each script is self-contained: it reads the oriGen questionnaire file, applies the same data cleaning and analytical population, and produces one or more outputs of the manuscript. Scripts can be run independently and in any order.

| Script | Output |
|---|---|
| `01_table1_characteristics.R` | Table 1. Characteristics of the study population, overall and by metabolic health status under Definition A and Definition B |
| `02_figure1_mh_classification.R` | Figure 1. Prevalence of MH and MUH and reclassification between definitions, overall and by BMI category |
| `03_figure2_component_contribution.R` | Figure 2. Relative contribution (general dominance) and adjusted ORs of the components of each definition for CVD, FLD, and CKD |
| `04_figure3_auc_definitions.R` | Figure 3. Discrimination (AUC) of Definition A and Definition B for CVD, FLD, and CKD, total and by BMI category |
| `05_figure4_auc_indices.R` | Figure 4. Discrimination (AUC) of 27 cardiometabolic indices for CVD, FLD, and CKD, total and by BMI category |
| `06_figure5_supp_table4_index_gradients.R` | Figure 5. Adiposity-related gradients in the predicted probability of FLD within metabolic phenotypes; Supplementary Table 4. Index-by-phenotype interaction and MH versus MUH contrasts for all index–outcome pairs |
| `07_supp_table1_diagnoses.R` | Supplementary Table 1. Prevalence of the individual diagnoses included in the composite outcomes, overall, by BMI category, and by sex |
| `08_supp_table3_by_sex_age.R` | Supplementary Table 3. Characteristics of the study population by sex and age group |

Formulas and references for the cardiometabolic indices are provided in Table 2 of the manuscript.

## How to run

1. Obtain the oriGen questionnaire data file through the oriGen Project (see *Data*).
2. In each script, set `data_path` to the location of the file:
   ```r
   data_path <- "path/to/oriGen_questionnaire_data.tsv"
   ```
3. Run the script. Figures are printed to the graphics device and tables to the console; no files are written.

`06_figure5_supp_table4_index_gradients.R` fits models for 26 indices × 3 outcomes × 2 definitions and may take several minutes.

## Requirements

- R version 4.6.0
- Packages: `dplyr`, `readr`, `stringr`, `tidyr`, `tibble`, `ggplot2` (≥3.5.0), `ggalluvial`, `patchwork`, `pROC`, `ggrepel`, `emmeans`

Install the packages with:

```r
install.packages(c("dplyr", "readr", "stringr", "tidyr", "tibble", "ggplot2",
                   "ggalluvial", "patchwork", "pROC", "ggrepel", "emmeans"))
```

## Citation

If you use this code, please cite the article above and the archived version of this repository:

> [Authors]. origen-metabolic-health-definitions (Version v1.0.0) [Computer software]. Zenodo. DOI: [pending]

## License

This code is released under the MIT License (see `LICENSE`).

## Contact

For questions about the code, contact [Adolfo Isaac Ruiz Ballesteros adolfo.ruba@gmail.com].
