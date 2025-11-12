# Scripts Used for Analysis of Human Data

## Quick Install Guide

### Prerequisites
Before running the MATLAB scripts, make sure you have:
- **MATLAB R2022b** or later (recommended)  
- The following MATLAB toolboxes installed:
  - Signal Processing Toolbox  
  - Statistics and Machine Learning Toolbox
  - Image Processing Toolbox
  - (Optional) Parallel Computing Toolbox for faster processing

### Installation Steps
These steps are estimated to take approximately 10 minutes.
1. **Clone or Download the Repository**
   ```bash
   git clone https://github.com/mhedlund/RNET_Public/CCDT/Human_Analysis.git
   cd Human_Analysis
   ```
   Or download the .zip file from GitHub and extract it.
2. **Download demo dataset from Google Drive**  
   Manually download the dataset with this URL: https://drive.google.com/uc?export=download&id=11e9hJd8AekxdcqgMp3ramicjCXZMltRt
   Manually download additional helper files (too large to upload to GitHub) and add to the 'helper' folder: https://drive.google.com/uc?export=download&id=1JenyyYrOU7TEMkJszQdRJo-6j-jq8ubc
4. **Add Project to MATLAB Path**  
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


<u>**CCDT_feature_select:**</u>  
Performs feature selection to identify which features are correlated with trial-by-trial response time.


<u>**CCDT_Classifier:**</u>  
Uses selected features as input for an SVM classifier to assess performance on classifying trials as fast or slow.  
This is used to determine the best feature space and to verify behavioral specificity with null comparisons.  


<u>**CCDT_feature_details:**</u>  
Summarizes information about selected features (subject, frequency band, electrode).  


<u>**CCDT_anatomical_analysis:**</u>  
Analyzes anatomical representation of selected features across subjects.  


<u>**CCDT_subregion_proportion_test:**</u>  
Performs statistical tests to determine which subregions are over- or under-represented within the selected feature space.  

<u>**Helper Files:**</u>  
Will not be run independently — they are called within interactive files.  


---
## Demo & Analysis instuctions
Demo will go through full analysis using raw data from a single subject. This will take approximately 30-60 minutes to run (tested on Windows 11, Intel Core i9-12900K, 2500MHz, 14Cores, 16GB RAM). 

The demo can be expaneded to analyze the full dataset from FigShare. After downloading this dataset, simply change the raw data directories in `loadData_natus.m`, `CCDT_graph_features.m`, `CCDT_power_features.m`, and `CCDT_trPLV.m`. Then uncomment all lines in `CCDTdatabase.m` to include all subjects and task sessions in the analysis. 

### Pseudocode for Demo
```
% matlab pseudocode for running scripts and descriptive comments about the purpose of each step as well as the expected input and output.

%% 1. Define subjects for analysis 
% Inside CCDTdatabase.m: uncomment lines only for subjects/sessions to 
% include in analysis (for Demo, only line 1 should be uncommented).
% Save this file.

%% 2. Generate features 
% Preprocess raw iEEG data and calculate trial-by-trial features

% Inside CCDT_graph_features.m, CCDT_power_features.m: edit parameters 
% (such as time period in trial, frequency bands, preprocessing options) 
% and define output directory and file name to save features.

run("CCDT_graph_features.m")
% output: file saved with trial-by-trial graph communicability features and
% regression for 

run("CCDT_power_features.m")
% output: file saved with trial-by-trial spectral power features

% repeat for each cue to calculate features in preparatory and anticipatory periods
%% 3. Feature Selection
% Robust feature selection regression: finds features significantly correlated to
% trial-by-trial reaction time for random 80% of trials, bootstrapped 1000x)

% Inside CCDT_feature_select.m: edit parameters (such as percent of trials 
% to train on & # bootstrap iterations), define input files/directories 
% (communicability & power features saved in step 2), and define output
% directory and file name to save regressions results from each bootstrap
% iteration.

run("CCDT_feature_select.m")
% output: regression results for each bootstrap iteration for both power
% and communicability features.

% Identify best feature space: Define a feature space as features with a
% significant correlation to reaction time in at least X percent of
% bootstrap iterations. Use an SVM classifier to assess how well this
% feature space performs at classifying trials as "fast" or "slow". Assess
% this performance over multiple values of the threshold X. Choose your
% final X to be the value with the best performance that avoids overfitting 
% and underfitting. (In manuscript, we chose X = 30%.)

% Inside CCDT_Classifier.m: edit parameters, define input files/directories 
% (communicability & power features saved in step 2, regression results from 
% CCDT_feature_select.m), and define output directory and file name to save 
% SVM performance results.

thresholds = 5:5:100; % percentages
for X = thresholds
    CCDT_Classifier(X)
    % output: saves SVM performance (AUC) for threshold X
end

%% 4. Get details for selected features
% For optimal selected features, save electrode name, frequency band, 
% and average feature value across all trials. 

% Inside CCDT_feature_details.m, edit parameters (percent threshold), 
% define input files/directories (communicability & power features saved in
% step 2, regression results from step 3), and define output directory and 
% file name to save feature details.
run("CCDT_feature_details.m")
% output: feature details for both power and communicability features

%% 5. Composite anatomical region analysis
% Analyze anatomical relationships of selected features and identify
% features consistently correlated with RT across all subjects. Analyze
% intrinsic communicability values of anatomical regions.

% Inside CCDT_anatomical_analysis.m, edit parameters (percent threshold, 
% plot options), define input files/directories (communicability & power 
% features saved in step 2, regression results from step 3, feature details 
% from step 4), and define output directory and file name to save results.
run("CCDT_anatomical_analysis.m")
% output: save data from heatmap figure showing features robustly
% correlated with RT across all subjects and plot intrinsic communicability
% of each anatomical subregion

%% 6. Anatomical subregion analysis
% Analyze proportional representation of anatomical subregions with
% Chi-squared test of proportions.

% Inside CCDT_subregion_proportion_analysis.m: edit parameters (feature type), 
% define input files/directories (selected featuredetails from step 4).
run("CCDT_subregion_proportion_analysis.m")
% output: statistical test results for each subregion in unique_subregions variable

%% 7. Calculate time-resolved PLV
% Preprocess raw iEEG data and calculate time-resolved PLV

% Inside CCDT_trPLV.m: edit parameters (such as time period in trial, 
% frequency bands, preprocessing options) and define output directory and 
% file name to save data.
run("CCDT_trPLV.m")
% output: file saved with time-resolved PLV averaged over trials
% corresponding to the following conditons:
%   1. Fastest 1/3 of trials
%   2. Slowest 1/3 of trials
%   3. All non-error trials
%   4. Middle 1/3 of trials

% Repeat for each cue to calculate features in preparatory and anticipatory periods

%% 8. Time-resolved PLV Analysis
% Plot and statistically compare time-resolved PLV for fast and slow trial conditions
% for electrode pairs between any two anatomical or functional regions.

% Inside CCDT_trplv_analysis.m: edit parameters (such as regions of interest)
% and define input files/directories (time-resolved PLV data from step 7).
run("CCDT_trplv_analysis.m")

% Repeat for all pairs of interest. Adjust num_comp variable in connectivity tilt
% section based on number of regions analyzed for appropriate multiple comparison
% correction.

```

---

