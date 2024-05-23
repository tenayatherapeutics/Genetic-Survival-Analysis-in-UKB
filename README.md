# Genetic-Survival-Analysis-in-UKB

This repository contains an scripts for performing detailed survival analysis using Cox proportional hazards and Kaplan-Meier plots. The script includes functions for loading and preparing phenotypic data, matching cases and controls, and generating survival plots. Key functionalities include:

- **Data Loading and Preparation:** Efficiently reads and processes phenotype and genotype data from UK Biobank
- **Survival Plotting:** Generates survival plots for time from diagnosis to cardiovascular outcomes using both Cox proportional hazards and Kaplan-Meier methods.
- **Matching Function:** Implements nearest neighbor matching to balance treatment and control groups based on various covariates.
- **Power Calculation:** Includes functions to perform power calculations for survival analysis.

## Table of Contents

- [Dependencies](##Dependencies)
- [Data Extraction and Preparation from UK Biobank](##Preparing-and-Extracting-Data-from-UK-Biobank)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Dependencies
Dependencies include `tidyverse`, `data.table`, `survminer`, `survival`, `MatchIt`, `expss`, and `cobalt`.

## Preparing and Extracting Data from UK Biobank 

This section will focus on extracting and processing data from UK Biobank to retrieve hospital admission data, biomarker measurements, gene dosages, and much more. This section will be broken into two parts, one focusing on extracting data from DNANexus and one from Hospital Inpatient Data not included in the "main" file provided by UKB. More information on this hospital inpatient data (HES-in) could be found [here](https://biobank.ndph.ox.ac.uk/ukb/ukb/docs/HospitalEpisodeStatistics.pdf).

### UKB Research Analysis Platform/DNANexus

The UK Biobank Research Analysis Platform is required in order to extract the fields needed to build the basic phenotype file needed to run the analysis, or replicate the results presented in our paper. If repurposing this repository for other endpoints and diseases of interest, the fields list could be ammended to grab whichever ones provide the relevant traits. This step of data preparation can also be done on a local copy of the main file, though this section would only work with the UKB DNANexus platform.

#### **[Surv_Paper_Pheno.ipynb](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Survival_UKB_scripts/Surv_Paper_Pheno.ipynb)** Usage Guide:

- Simply upload this file into the project folder of choice on DNANexus, open a JupyterLab environment and execute the code provided.
- If there is a need to amend the fields extracted from the main dataset, simply amend this part of the script:
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
field_names
```
- Sometimes it would be faster to copy the output from ```field_names``` and manually remove any instances which are not needed for the analysis, especially if only measurements from the initial enrollment into UKB are being considered. This would greatly increase runtime for the current implementation  


