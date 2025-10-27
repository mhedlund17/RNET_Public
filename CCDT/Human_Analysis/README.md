# Scripts Used for Analysis of Human Data

## Quick Install Guide

### Prerequisites
Before running the MATLAB scripts, make sure you have:
- **MATLAB R2022b** or later (recommended)  
- The following MATLAB toolboxes installed:
  - Signal Processing Toolbox  
  - Statistics and Machine Learning Toolbox  
  - (Optional) Parallel Computing Toolbox for faster processing

### Installation Steps
These steps are estimated to take approximately 10 minutes.
1. **Clone or Download the Repository**
   ```bash
   git clone https://github.com/mhedlund/RNET_Public/CCDT/Human_Analysis.git
   cd Human_Analysis
   ```
   Or download the .zip file from GitHub and extract it.  
2. **Add Project to MATLAB Path**  
   Open MATLAB and run:  
   ```bash
   addpath(genpath(pwd)) % pwd prints current working directory
   savepath
   ```  
   This makes all functions and scripts accessible.

---
## File Descriptions (listed in approximate order of use)

<u>**CCDTdatabase:**</u>  
Lists subjects and task sessions to include in the analysis. This file is loaded into all analysis functions/scripts. Edit this prior to any analysis, uncommenting lines corresponding to sessions you wish to include.


<u>**loadCCDTdata:**</u>  
Processes iEEG data (called from `CCDT_graph_features` and `CCDT_power_features`).  
Can be edited to adjust preprocessing steps.  


<u>**CCDT_graph_features:**</u>  
Calculates and saves trial-by-trial graph network metrics (communicability).  


<u>**CCDT_power_features:**</u>  
Calculates and saves trial-by-trial spectral power metrics.  


<u>**CCDT_trPLV:**</u>  
Calculates and saves time-resolved phase locking value.  


<u>**CCDTanalyze:**</u>  
Performs feature selection to identify which features are correlated with trial-by-trial response time.  
&nbsp;&nbsp;• 1st section performs regression analysis.  
&nbsp;&nbsp;• 2nd section summarizes information about selected features (subject, frequency band, electrode).  


<u>**CCDT_Classifier:**</u>  
Uses selected features as input for an SVM classifier to assess performance on classifying trials as fast or slow.  
This is used to determine the best feature space and to verify behavioral specificity with null comparisons.  


<u>**CCDT_anatomical_analysis:**</u>  
Analyzes anatomical representation of selected features across subjects.  


<u>**CCDT_subregion_proportion_test:**</u>  
Performs statistical tests to determine which subregions are over- or under-represented within the selected feature space.  

<u>**Helper Files:**</u>  
Will not be run independently — they are called within interactive files.  


---
## Demo & Analysis instuctions
Demo will go through full analysis using raw data from a single subject. This will take approximately XX minutes/hours to run (tested on Windows 11, Intel Core i9-12900K, 2500MHz, 14Cores, 16GB RAM). 

The demo can be expaneded to analyze the full dataset from FigShare. After downloading this dataset, simply change the raw data directories in `loadData_natus.m`, `CCDT_graph_features.m`, `CCDT_power_features.m`, and `CCDT_trPLV.m`. Then uncomment all lines in `CCDTdatabase.m` to include all subjects and task sessions in the analysis. 

# Pseudocode for Demo
```
%add matlab pseudocode for running scripts and descriptive comments about the purpose of each step as well as the expected input and output.
---

