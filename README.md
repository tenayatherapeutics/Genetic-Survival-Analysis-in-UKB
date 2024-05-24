# Genetic-Survival-Analysis-in-UKB

This repository contains an scripts for performing detailed survival analysis using Cox proportional hazards and Kaplan-Meier plots. The script includes functions for loading and preparing phenotypic data, matching cases and controls, and generating survival plots. Key functionalities include:

- **Data Loading and Preparation:** Efficiently reads and processes phenotype and genotype data from UK Biobank
- **Survival Plotting:** Generates survival plots for time from diagnosis to cardiovascular outcomes using both Cox proportional hazards and Kaplan-Meier methods.
- **Matching Function:** Implements nearest neighbor matching to balance treatment and control groups based on various covariates.
- **Power Calculation:** Includes functions to perform power calculations for survival analysis.

## Table of Contents

- [Dependencies](##Dependencies)
- **[Step 1: Data Extraction and Preparation from UK Biobank](##Step-1-Preparing-and-Extracting-Data-from-UK-Biobank)**
  - [Extracting Fields for Basic Phenotype Preparation](####Extracting-Fields-for-Basic-Phenotype-Preparation)
  - [Extracting SNP Dosages](####Extracting-SNP-Dosages)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Dependencies
Dependencies include:
- `tidyverse`
- `data.table`
- `survminer`
- `survival`
- `MatchIt`
- `expss`
- `cobalt`	

## **Step 1: Preparing and Extracting Data from UK Biobank** 

This section will focus on extracting and processing data from UK Biobank to retrieve hospital admission data, biomarker measurements, gene dosages, and much more. This section will be broken into two parts, one focusing on extracting data from DNANexus and one from Hospital Inpatient Data not included in the "main" file provided by UKB. More information on this hospital inpatient data (HES-in) could be found [here](https://biobank.ndph.ox.ac.uk/ukb/ukb/docs/HospitalEpisodeStatistics.pdf).

### UKB Research Analysis Platform/DNANexus:

The UK Biobank Research Analysis Platform is required in order to extract the fields needed to build the basic phenotype file needed to run the analysis, or replicate the results presented in our paper. If repurposing this repository for other endpoints and diseases of interest, the fields list could be ammended to grab whichever ones provide the relevant traits. This step of data preparation can also be done on a local copy of the main file, though this section would only work with the UKB DNANexus platform.
___
#### **Extracting Fields for Basic Phenotype Preparation**:

Upload the notebook [```Surv_Paper_Pheno.ipynb```](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Survival_UKB_scripts/Surv_Paper_Pheno.ipynb) into the project folder of choice on DNANexus, open a JupyterLab environment and execute the code provided. If there is a need to amend the fields extracted from the main dataset, simply amend this part of the script:
```python
field_ids = ['22006', #Genetic Ethnic Grouping
             '31', #Sex
             '21022',#Age at Recruitment
             '22009',#Genetic PCs
             '4080', #Systolic Blood Pressure
             '52', #Month of Birth
             '34', #Date of Birth
             '30710', #CReactive Protein
             '20116', #Smoking Status
             '21001',#BMI
             '30780',#LDL
             '20003',#Medication
             '137', #Medication_Groups
             '30700' #Creatinine
            ] 
field_names = field_names_for_ids(field_ids)
field_names
```
Sometimes it would be faster to copy the output from ```field_names``` and manually remove any instances which are not needed for the analysis, especially if only measurements from the initial enrollment into UKB are being considered. This would greatly increase runtime for the current implementation.  
___
#### **Extracting SNP Dosages**:
Upload the notebook [```Dosage_Extraction.ipynb```](provide_link) into the project folder of choice on DNANexus, open a JupyterLab environment, and execute the code provided. 
___

## **Step 2: Obtaining Patient Medication Use Data**

This step focuses on grouping patients by medication usage, specifically into medicine classes related to cardiovascular disease (such as Statins, ACE inhibitors, Beta-Blockers, etc). Medications are classified into bins defined in the [`MedicationGroups.txt`](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Medications/MedicationGroups.txt) file. The resulting file should contain the basic phenotypes obtained in the last step, as well as columns indicating whether each patient self-reports usage for the medicine class defined. 

### **Grouping UKB Medicine Codes**

The python script [`medication_grouping.py`](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Medications/medication_grouping.py) groups UKB medication codes into bins using the medication groups file and a data coding file obtained from UK Biobank [`Data Coding 4`](https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=4). Due to the fact that in UKB similar medications with different dosages are treated as entirely different medications, with separate numerical codes, this script aims to find keywords in the coding file which can be used to classify aforementioned medication. 

Place this script and all other required files in the same directory and run. Required libraries and files can be seen here: 
```
#!/usr/bin/python3

import sys 
import subprocess
import os 
import re

med_groups = 'MedicationGroups.txt'
codings = 'coding4.tsv'
```
