# RakaiAnalysisPlosMed2026

This repository contains the analytic code used for the Rakai analysis prepared for submission to *PLOS Medicine* (2026).

The purpose of this repository is to provide transparent, reproducible analysis code supporting the results presented in the manuscript. The code is written primarily in R and is organized to separate data preparation, analysis, and figure/table generation. Input data for Rakai_cohort_code_PLOS_submission.R cannot be shared publicly but may be requested from datarequests@rhsp.org.

---

## Repository structure

RakaiAnalysisPlosMed2026/  
├── Rakai_cohort_code_PLOS_submission.R                 # R script for cohort data processing and statistical model runs for incidence, prevalence, ART, VMMC trends etc. 

├── Rakai_model_scenario_comparison_PLOS_submission.R   # R script for plotting EMOD output

├── README.md  

└── .gitignore  

Folder names may evolve as the analysis is finalized.

---

## Software requirements

- R (version ≥ 4.2 recommended)
- R packages as specified in the scripts  
  (optionally managed using `renv`)

If `renv` is used, package versions can be restored by running `renv::restore()`.

---

## Running the analysis

The analysis is designed to be run in a modular fashion.

Typical workflow:
1. Run analysis scripts to generate incidence, prevalence, ART etc. estimates
2. Plot outputs from mathematical model (EMOD)
---

## Data availability

Due to ethical, legal, or data-use restrictions, the primary input data used in this analysis may not be publicly shareable.
Where data cannot be shared:
- Cohort Script assumes the presence of locally stored input files.  This datafile may be provided upon request. 
- Model comparison script uses csv files included in the data repository
- Variable names, formats, and expected structures are documented in code comments
- Synthetic or placeholder data may be added where helpful to demonstrate workflow.

---

## Reproducibility notes

- All random processes use fixed seeds where applicable
- Results presented in the manuscript correspond to the code state at submission
- This repository reflects the final analytic approach used for publication

---

## Citation

If you use or adapt this code, please cite the associated *PLOS Medicine* article (citation to be added upon publication).

---

## Contact

For questions about the analysis or code, please contact:

**Adam Akullian**  
adam.akullian@gatesfoundation.org
Gates Foundation
