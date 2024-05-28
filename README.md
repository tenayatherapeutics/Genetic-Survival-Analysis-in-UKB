# Genetic Survival Analysis in UKB

This repository contains an scripts for performing detailed survival analysis using Cox proportional hazards and Kaplan-Meier plots. The script includes functions for loading and preparing phenotypic data, matching cases and controls, and generating survival plots. Key functionalities include:

- **Data Loading and Preparation:** Efficiently reads and processes phenotype and genotype data from UK Biobank
- **Survival Plotting:** Generates survival plots for time from diagnosis to cardiovascular outcomes using both Cox proportional hazards and Kaplan-Meier methods.
- **Matching Function:** Implements nearest neighbor matching to balance treatment and control groups based on various covariates.
- **Power Calculation:** Includes functions to perform power calculations for survival analysis.

## Table of Contents

<details>
<summary>Click to expand</summary>
  
- [Dependencies](#dependencies)
- [Step 1: Data Extraction and Preparation from DNANexus](#step-1-preparing-and-extracting-data-from-dnanexus)
  - [Extracting Fields for Basic Phenotype Preparation](#extracting-fields-for-basic-phenotype-preparation)
  - [Extracting SNP Dosages](#extracting-snp-dosages)
- [Step 2: Obtaining Patient Medication Use Data](#step-2-obtaining-patient-medication-use-data)
  - [Grouping UKB Medicine Codes](#grouping-ukb-medicine-codes)
  - [Creating the Binary Medication Variables](#creating-the-binary-medication-variables)
- [Step 3: Formatting Hospital Inpatient Data](#step-3-formatting-hospital-inpatient-data)
  - [Required Files and Libraries](#required-files-and-libraries)
  - [Main Functions](#main-functions)
- [Step 4: Running Survival Analysis](#step-4-running-survival-analysis)
  - [Required Files and Libraries](#required-files-and-libraries-1)
  - [Main Functions](#main-functions-1)

</details>

## Dependencies
Dependencies include:
- [`tidyverse`](https://cran.r-project.org/package=tidyverse)
- [`data.table`](https://cran.r-project.org/package=data.table)
- [`survminer`](https://cran.r-project.org/package=survminer)
- [`survival`](https://cran.r-project.org/package=survival)
- [`MatchIt`](https://cran.r-project.org/package=MatchIt)
- [`expss`](https://cran.r-project.org/package=expss)
- [`cobalt`](https://cran.r-project.org/package=cobalt)

To Install: 
```R
# List of required packages
packages <- c("tidyverse", "data.table", "survminer", "survival", "MatchIt", "expss", "cobalt")

# Check if packages are installed and install missing ones
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
```

## Overview

Steps 1-3 will focus on extracting and processing data from UK Biobank to retrieve hospital admission data, biomarker measurements, gene dosages, and much more. This section will be broken into two parts, one focusing on extracting data from DNANexus (Step 1-2) and one from Hospital Inpatient Data (Step 3) not included in the "main" file provided by UKB. More information on this hospital inpatient data (HES-in) could be found [here](https://biobank.ndph.ox.ac.uk/ukb/ukb/docs/HospitalEpisodeStatistics.pdf). Step 4 describes how to use the files created in Steps 1-3 to perform survival analysis. 

## **Step 1: Preparing and Extracting Data from DNANexus** 

The UK Biobank Research Analysis Platform is required in order to extract the fields needed to build the basic phenotype file needed to run the analysis, or replicate the results presented in our paper. If repurposing this repository for other endpoints and diseases of interest, the fields list could be ammended to grab whichever ones provide the relevant traits. This step of data preparation can also be done on a local copy of the main file, though this section would only work with the UKB DNANexus platform.

### **Extracting Fields for Basic Phenotype Preparation**:

Upload the notebook [```Surv_Paper_Pheno.ipynb```](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Survival_UKB_scripts/Surv_Paper_Pheno.ipynb) into the project folder of choice on DNANexus, open a JupyterLab environment and execute the code provided. If there is a need to amend the fields extracted from the main dataset, simply amend this part of the script:
```python3
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

### **Extracting SNP Dosages**:

Upload the notebook [```Dosage_Extraction_UKB.ipynb```](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Survival_UKB_scripts/Dosage_Extraction_UKB.ipynb) into the project folder of choice on DNANexus, open a JupyterLab environment, and execute the code provided. One or more paths to the correct dosage files and corresponding sample files on UKB must be provided for this script to work, as seen here in this example: 

```python3
# Get TOPMED imputed dosages (per chromosome)
def getImputedPath(chrom):
    imputed_path = '/mnt/project/Bulk/Imputation/Imputation from genotype (TOPmed)/ukb21007_c' + chrom[3:] + '_b0_v1.bgen'
    return imputed_path

# Return Sample Dataframe for Imputed dosages
def getSampleDF(chrom):
    sample_path = '/mnt/project/Bulk/Imputation/Imputation from genotype (TOPmed)/ukb21007_c' + chrom[3:] + '_b0_v1.sample'
    sample_df = pd.read_csv(sample_path, sep=' ')
    sample_df = sample_df.drop(index=0)
    sample_df = sample_df.drop(['missing', 'sex'], axis=1)
    sample_df = sample_df.reset_index(drop=True)
    sample_df = sample_df.rename(columns={"ID_1": "FID", "ID_2": "IID"}, errors="raise")
    return sample_df

```
To use the main dosage extraction function, simply pass a list of SNPs formatted in the same way as the example provided below. **The format for denoting chromosome number in WES data on UKB is "chr10" while in TOPMed imputed data it is simply "10".**
```python3

# ADRB1 Example
# chr:pos:ref:alt:rsnum

# WES Format
ADRB1_snps = ['chr10:114044277:A:G:rs1801252', 
             'chr10:114045297:G:C:rs1801253']
# TOPMed Imputed Format
ADRB1_snps = ['10:114044277:A:G:rs1801252', 
             '10:114045297:G:C:rs1801253']
```
Once dosages for the SNPs of interest are successfully obtained, write the resulting dataframe to a `csv` or `txt`. This can be done in DNANexus like so:

```python3
file_path = 'dosage.txt'
multipleSNPdosages.to_csv(file_path, sep='\t', index=False)

%%bash 
dx upload --destination ./path/to/working/dir/ *.txt
```

### Expected Outputs 
- **dosage.txt**: Dosage file consisting of IIDs and RSNUMs.
- **example_survival_base_pheno.txt**: Basic phenotype file with measurements available in the main dataset, columns required to extract medication use data, and covariates such as age, sex, and principle components. 
___

## **Step 2: Obtaining Patient Medication Use Data**

This step focuses on grouping patients by medication usage, specifically into medicine classes related to cardiovascular disease (such as Statins, ACE inhibitors, Beta-Blockers, etc). Medications are classified into bins defined in the [`MedicationGroups.txt`](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Medications/MedicationGroups.txt) file. The resulting file should contain the basic phenotypes obtained in the last step, as well as columns indicating whether each patient self-reports usage for the medicine class defined. 

### **Grouping UKB Medicine Codes**

The python script [`medication_grouping.py`](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Medications/medication_grouping.py) groups UKB medication codes into bins using the medication groups file and a data coding file obtained from UK Biobank [`Data Coding 4`](https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=4). Due to the fact that in UKB similar medications with different dosages are treated as entirely different medications, with separate numerical codes, this script aims to find keywords in the coding file which can be used to classify aforementioned medication. 

Place this script and all other required files in the same directory and run. Required libraries and files can be seen here: 
```python3
#!/usr/bin/python3

import sys 
import subprocess
import os 
import re

med_groups = 'MedicationGroups.txt'
codings = 'coding4.tsv'
```
### **Creating the Binary Medication Variables** 

The R script [`Append_Med_to_Pheno`](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Survival_UKB_scripts/Append_Med_to_Pheno.R) will take the outputs from previous steps and create the binary variables needed for downstream analysis. The required files can be seen here: 
```R
#### Basic phenofile and requisite data preparations ----

MedGroup <- fread('MedicationGroups.txt')
codings <- fread('coding4.tsv')
result <- fread('merged_med_groups_stringent_anticoag.txt') #Result from running the medications.py script 
base_pheno <- fread('example_survival_base_pheno.txt')
```

### Expected Outputs
- **example_pheno_with_meds.txt**: Phenofile created from Step 1 merged with binary medication status columns created in this step.

The resulting dataframe at the end of this step should look something like this (with many more columns): 

| IID  | Trait 1 | Statins | ACE Inhibitors | ...|
| ---| ---| ---| ---| ---|
| 1 | 32 | 1 | 0 | ... |
| 2  | 51  | 0  | 1 | ... |

___
## **Step 3: Formatting Hospital Inpatient Data** 

This next step will format the hospital inpatient data available on UK Biobank, extracting diagnosis dates for ICD9/10 codes of interest so that time to event can be calculated for downstream survival analysis. The R script [`retrieve_diag_dates_UKB.R`](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Survival_UKB_scripts/retrieve_diag_dates_UKB.R) will take the raw inpatient data provided by UKB and transform it into a mixture of date/time, binary, and numerical variables used to fit Cox Proportional Hazard and Kaplan Meier survival models. Required files and libraries are listed here:

### Required Files and Libraries

```R
library(tidyverse)
library(data.table)
library(ukbtools)
library(tidyfast)
library(stringr)
library(tidyfst)

big <- fread('~/ukbiobank/Paper_Scripts/datafiles/example_pheno_with_meds.txt') #Phenofile with base traits and binary medication variables

#### Read in Dataframes for HESIN record-level data obtained from UK Biobank----

hesin_all <- fread('hesin.txt') #Main HESIN file with patient hospital stay data
hesin_diag <- fread('hesin_diag.txt') #HESIN file with ICD9/10 diagnosis codes per patient
hesin_operations <- fread('hesin_oper.txt') #HESIN file with OPCS4 diagnosis codes per patient

death_date <- fread('death.txt') #HESIN file listing death dates for applicable patients 
death_cause <- fread('death_cause.txt') #HESIN file listing ICD9/10 causes of death for corresponding patients 

```
### Main Functions

The script utilises two main functions in order to extract diagnosis dates from the HESIN data, which are described below. The second function `get_nth_diag_date_after` was implemented to look for the date of the second independant follow-up hospital admission for the same disease. Both return a dataframe with `IID/EID` as the first column and `IDates containing diagnosis dates` in the second column. 

```R
# Retrieves the first diagnosis date of a specific disease, defined by a group of ICD codes:

get_first_diag_date(all = hesin_all, # First three required inputs correspond to HESIN files 
                    diag = hesin_diag,
                    operations = hesin_operations,
                    icd_codes = afib_codes, # List of ICD codes used to diagnose disease of interest
                    opcs4_codes = afib_opcs4, # List of OPCS4 codes used for disease definition
                    pheno_name = 'Afib') # Name of disease 


# Retreives the first available diagnosis date of a specific disease after N number of days from first listed diagnosis date:

get_nth_diag_date_after(all = hesin_all, # First six fields are identical to the previous function
                        diag = hesin_diag,
                        operations = hesin_operations,
                        icd_codes = afib_codes,
                        opcs4_codes = afib_opcs4,
                        pheno_name = 'Afib',
                        n = days) # Number of days after the first diagnosis date 
```
### Expected Outputs 

Successfully running this script should yield the following files: 
- **pheno_with_diagnosis_dates.txt**: File with measurements, binary medications, and newly appended diagnosis dates for various diseases of interest 
- **thirty_day_rehos.txt**: File with dates corresponding to rehospitalization after thirty days 
- **before_AF.txt**: File with other disease onset status before first diagnosis date of Atrial Fibrillation
- **before_HF.txt**: File with other disease onset status before the first diagnosis date of all-cause Heart Failure
- **before_IHD.txt**: File with other disease onset status before the first diagnosis of Ischemic Heart Disease
- **cv_outcome_raw.txt**: File with outcome dates corresponding to date of certain cardiovascular operations, death, or study end (denoted as Sept 1, 2023)

___
## Step 4: Running Survival Analysis

After extracting and formatting patient data from UK Biobank survival analysis can be performed to determine whether specific SNP genotype carriers have increased or decreased risk for certain diseases. An example runthrough with all the required files produced in previous steps can be found in the [`Example`](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/tree/main/Example) folder in this repository. To go through the example, run `Example.R` in the aforementioned folder. Ad-hoc power calculations are also performed at the end of the script, but can be ignored if not needed. 

### Required Files and Libraries
```R
library(tidyverse)
library(data.table)
library(survminer)
library(survival)
library(MatchIt)  
library(expss)

#### Load phenotype file with diagnosis dates 
set.seed(1)
raw_master_pheno <- fread('pheno_with_diagnosis_dates.txt')

#### Outcome Dates File, created in Step 3
death_dates <- fread('cv_outcome_raw.txt')

#### Dosage file for SNP of interest
ACE_dosage <- fread('example_dosage.txt')

#### Rehospitalization dates
second_hosp <- fread('thirty_day_rehos.txt')

#### Status of other diseases before first diagnosis date of disease used as endpoint
before_HF <- fread('before_HF.txt') 

```
In the example script the files corresponding to `raw_master_pheno` and `second_hosp` are merged into one named `example_pheno_with_AF_rehospitalization.txt`, which why a separate dataframe for rehospitalization dates is not present in the example script. 

### Main Functions 

The script utilises two main functions to perform matching and survival model fitting/plotting; `plotsurv` and `matchit_all`. 

```R
# Matching function: Matches individuals in the treatment and control groups, denoted by genotype, returns dataframe with matched patients.

matchit_all(subset_df = pheno_filtered, # subset_df needs the dataframe with the participants and columns for input
            list_of_conditions = list(target_df, age_pheno_diag),   # list of conditions includes other columns which you may want to add 
            phenotype = NA, # set phenotype to NA generally, if not it'll grab columns automatically for phenotype of interest 
            SNP = rsnum, # SNP variable takes the genotype variable in binary format as the matching variable (treatment vs controls)
            omit_cols = c(example_matching_terms) # omit_cols takes a vector of columns names to omit
)  

# Plotting and model fit function for Cox Proportional Hazards.
# Prints counts of individuals experiencing the outcome of interest, grouped by genotype, as well as the CoxPH model summary. 
 
plotsurv(rsnum = 'raw', 
         gene_name = gn, 
         pheno_df = pheno_filtered, 
         phenotype = eval(pheno1), 
         endpoint = 'CV_outcome', 
         cut_time = 10, # cut_time represents max time to outcome, patients with times to outcome greater than this number will be coded as controls
         ylim = 0.0, # Ylim represents the floor of the yaxis 
         adj=TRUE) # Adj=True would add adjustments to the formula, defaulting as age and sex
```
### Expected Outputs

- Survival plots and model summaries for unmatched and matched data
___
## Example 

As previously mentioned, an example runthrough of the main survival pipelines (excluding data filtering and creation) can be found in the [`Example`](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/tree/main/Example) folder. The example workflow follows the analysis of the effects of `SNP1` on time to severe cardiovasular outcome from first diagnosis of Atrial Fibrillation (AF). This folder also contains the requisite files needed for running survival analysis, scrambled data files mimicking the files created if steps 1-3 are followed successfully. These include:
- ['example_CV_outcome.txt'](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Example/example_CV_outcome.txt): Example file containing the dates for death, study end, and cardiovascular operations.
- ['example_before.txt'](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Example/example_before.txt): Example file containing disease diagnosis status for various comorbidities, calculated as binary variables of whether the diagnosis date was before diagnosis of AF. 
- ['example_dosage.txt'](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Example/example_dosage.txt): Example dosage file.
- ['example_pheno_with_AF_rehospitalization.txt'](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/Example/example_pheno_with_AF_rehospitalization.txt): Example file containing baseline demographics information, dates for first disease diagnosis, and dates for first Atrial Fibrillation rehospitalizations (after 30 days from first diagnosis).

To begin the tutorial, place all files into a folder named `Example` (which automatically should be the case if using `git clone`), and load in the required files: 
```R
setwd('./Example/')
library(tidyverse)
library(data.table)
library(survminer)
library(survival)
library(MatchIt)  
library(expss)
library(cobalt)
```


## Contributing

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/tenayatherapeutics/Genetic-Survival-Analysis-in-UKB/blob/main/LICENSE) file for details.

