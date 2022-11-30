# biomaker-recommendation-system

This repository contains codes used in the paper 'A knowledge integration strategy for the selection of a robust multi-stress biomarkers panel for Bacillus subtilis'. This Biomarker Recommendation System takes a pool of candidate biomarker panels as input and recommend a biomarker panel based on a set of metrics that evaluate the relevance of these biomarkers to their stress sensing power. 

In */data* directory, a gene expression dataset *NicolasGeneExpression_processed.csv* which used for all models training is provided. It was also priorly used in https://github.com/neverbehym/transcriptional-biomarkers-subtilis to identify different cellular states and produce candidate biomarker panels specific to these states. You can also find *MultiClassLabel.csv* which assign samples to different cellular states and *GeneRegulatoryNetwork.csv* which includes the Gene Regulatory Network information in the same directory. 

The gene expression datasets derived from Nicolas data by stress condition are put in */data/NicolasData* and used to calulate the stress sensing index for Biomarker Recommendation System training. The external gene expression datasets covering various stress conditions are put in */data/ExternalData* and used to validating the Biomarker Recommendation System.

This Biomarker Recommendation System Computes the evaluation metrics from Biomarker Identification Model (BIM), Gene Regulatory Network (GRN), Co-Expression Network (CEN) in R 4.0.2 and writes the outputs in */results*. Please run:
```
Rscript BIM.R [Panellog_FILE] [data_FILE] [multiclass_FILE] 
Rscript GRN.R [Panellog_FILE] [GeneRegulatoryNetwork_FILE]
Rscript CEN.R [Panellog_FILE] [data_FILE] 
```
We have provided the example files, for example:
```
Rscript BIM.R ../results/CandidateBiomarkerList.csv ../data/NicolasGeneExpression_processed.csv ../data/MultiClassLabel.csv
```
The Stress Sensing Model (SSM) that tests the predictive power of a given biomarker panel to sense different stress conditions is implemented in Python 3.7. The outputs can be found in */results*. Please run:
```
Python SSM.py [Panellog_FILE] [Panel_ID] [Genelist_FILE] [data_DIR] [random/biomarker]
```

To compute the stress sensing indice of 1st candidate biomarker panel  using example training set from Nicolas data:
```
Python SSM.py  '../results/CandidateBiomarkerList.csv' 1 '../data/GeneAnnotations_full.csv' '../data/NicolasData' biomarker
```
To compute the stress sensing indice of a random set of genes of same size as 1st candidate biomarker panel:
```
Python SSM.py  '../results/CandidateBiomarkerList.csv' 1 '../data/GeneAnnotations_full.csv' '../data/NicolasData' random
```
To compute the stress sensing indice of 1st candidate biomarker panel using example test set from external datasets:
```
Python SSM.py  '../results/CandidateBiomarkerList.csv' 1 '../data/GeneAnnotations_full.csv' '../data/ExternalData' biomarker
```
A Jupyter notebook *recommendation_system.ipynb* showed the examination of these evaluation metrics and the training process of the recommendation system as well as visualisation of results.
