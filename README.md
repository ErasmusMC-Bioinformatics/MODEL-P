# MODEL-P

## Introduction
Here we give the code of the implementation of Multi-Omics DEep Learning for Prognosis-correlated subtyping (MODEL-P). MODEL-P is described in the study **Robust Deep Learning Model for Prognostic stratification of pancreatic ductal adenocarcinoma patients**, and could be used to identify prognosis-correlated PDAC subtypes and to predict prognoses of new patients.

## Data
MODEL-P was developed based on three types of multi-omics data (mRNA-sequencing, DNA methylation array and microRNA-sequencing) of 146 PDAC patients obtained from the TCGA cohort.
All the steps for the training data downloading and preprocessing was presented in the R script **TrainingData_preprocess.R**.

Afterward the prognosis-correlated subtypes identified from MODEL-P subtypes were validated on five independent datasets. Three test sets werer dpwnloaded from the ICGC Australian cohort (one mRNA-seq, one mRNA microarray, and one DNA methylation dataset). Two of these test sets came from one cohort of the Gene Expression Omnibus database (one mRNA array GSE62452 and one microRNA dataset GSE62498). The R script used to preprocess these test sets is **TestData_preprocessing.R**.

## Model
MODEL-P consist of two parts. Firstly, the PDAC prognosis-correlated subtypes were identified based on multi-omics using Autoencoder. The Python code was implemented in **AE-subtyping.ipynb**. 

Secondly, classification models were trained based on each single omics TCGA training set with the identified subtypes as class labels and validated on five independent test sets. The subtype-specific feature were selected using Python code **subtype_signatures.ipynb**. The SVM classifiers were trained and constructed using **test_model_prediction.ipynb**.

## Code dependencies
1. tensorflow 1.15.0
2. tensorflow.keras 2.2.4-tf
3. numpy 1.19.5
4. pandas 1.1.2
5. sklearn 0.23.1
