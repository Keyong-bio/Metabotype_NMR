# Metabotype_NMR
This is the repository (code) for the paper titled: "Metabotype identification via fasting and postprandial metabolomics and its association with type 2 diabetes incidence",
published on doi: 10.1186/s12933-025-02821-6, Cardiovasc Diabetol. 2025 Jul 8.

For more related questions, please feel free to send email to k.deng@lumc.nl.

Content Guide:

The Code Folder contains the step to reproduce the results.
- 00_Prepration.R: Data cleaning and metabolites measurement transformation
- 01_IPVs.R: IPVs selecting principal metabolites as inputs for consensus clustering analysis
- 02_baseline_table.R: This script describe the baseline character of participants included, across each cluster
- 03_Cox_model.R: This script conducts cox regression model for clusters (related with T2DM incidence) and how to evaluate the predictive value of cox model
- 04_Comparison_Metabolites.R: Plot the key metabolites for metabotypes with gradient T2DM risk
- 05_Comparison_Diet_score.R: comparison of dietary pattern scores across different metabotypes
- 06_SHAP_Plot_code: using python to run the random forest multi-classification model

In addition, you also need to source two functions: (using source function to load them)
- ConsensusCluster_function.R: to determine the optimal number of clusters
- KNN_Functions.R: to perform the KNN imputation for metabolomics data (fasting, postprandial and delta state)

Other .Rdata data are the inputs you need to use. 
- Due to the data restriction, we can not share the .Rdata publicly. However, the pipeline of code can be applied successfully if you use your own dataset.
