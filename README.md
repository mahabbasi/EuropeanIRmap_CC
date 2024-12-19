# Analysis Workflow for studying impact of climate change on streamflow intermittence in Europe

This repository contains the code and resources used to conduct the analysis for the paper:  
**"[Increased Streamflow Intermittence in Europe due to Climate Change Projected by Combining Global Hydrological Modeling and Machine Learning]"**  
Submitted in *[Earth's Future]*, [2024].  
Authors: [Mahdi Abbasi, Mathis Loïc Messager, Petra Döll]

---

## **Overview**
This repository documents the analytical workflow and methods used in the paper, providing a step-by-step guide to reproduce the results. The repository includes scripts for data preprocessing, exploratory data analysis, statistical modeling, and visualization. In this paper, we have studied the impact of climate change on streamflow intermittence over European reaches. The paper has been submitted to the Earth's Future (AGU) journal. 
To see the final outcomes (i.e., intermittence status of 1.5 million European reaches under climate change in the future) and the shapefiles, please visit the {}.

---
## **Code Workflow**

The workflow is structured into the following steps:

1. **Data Collection**:
    - Main functions for high-resolution predictors: `extract_dsflow_hr`, `compute_hr_predictors`, `compute_hr_interannual_predictors`
    - **Purpose**: Extract the monthly downscaled streamflow of the five GCMs and prepare the seven dynamic HR predictors
    -Main functions for low-resolution predictors: `select_modify_lr_predictors`
    - **Purpose**: Prepare the two dynamic LR predictors
    - **Inclusion**: (default: `FALSE`) due to the size of the predictors, we decided to exclude the data collection (`prepare_hr_predictors`, `prepare_lr_predictors`) plans in `_targets.R` script. The reason for this is that the raw files are not be able to be shared because of their size, in total over 3 TB. 
2. **Model execution**:    
    - Main functions: `runmodels_over_period`
    - **Purpose**: Apply the trained RF models for the five GCMs for three different periods: 1)reference period (1985-2014), 2) 2050s (2041-2100), and 3) 2080s (2071-2100).
    - **Inclusion**: (default: `FALSE`) If the raw data is available and the two previous steps have been completed, we can proceed with this step.Please note that your system should have at least 18 processors available, as the Random Forest (RF) model was trained in parallel using 18 processors.
2. **Post-processing**:  
    - Main functions: Functions in `post_processing_functions.R`
    - **Purpose**: Compute the five ecologically relevant indicators and then generate the shapefiles and figures for the paper.
    - **Inclusion**: (default: `TRUE`) the link for the required data is [link name] {url}



## **Citation**
If you use this repository or its resources, please cite the following paper:

## **Contact**
For questions or issues, please contact:

    [Mahdi Abbasi] (abbasi@em.uni-frankfurt.de)

