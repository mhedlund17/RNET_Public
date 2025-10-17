# Scripts Used for Analysis of Human Data

## File Descriptions (listed in approximate order of use)


<u>**loadCCDTdata:**</u>  
Processes iEEG data (called from `CCDT_graph_features` and `CCDT_power_features`).  
Can be edited to adjust preprocessing steps.  


<u>**CCDT_graph_features:**</u>  
Calculates and saves trial-by-trial graph network metrics (communicability).  


<u>**CCDT_power_features:**</u>  
Calculates and saves trial-by-trial spectral power metrics.  


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


##  Helper Files
Will not be run independently — they are called within interactive files.  
