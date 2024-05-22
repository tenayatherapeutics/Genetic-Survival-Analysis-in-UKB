.libPaths('/home/lzhang@TENAYA.local/my_bioinfo/rpackages')

library(tidyverse)
library(data.table)
library(ukbtools)
library(tidyfast)
library(stringr)
library(tidyfst)

setwd('~/ukbiobank/Paper_Scripts/datafiles/HES/')

big <- fread('~/ukbiobank/Paper_Scripts/datafiles/example_pheno_with_meds.txt')
#### Read in Dataframes for HESIN record-level data ----

hesin_all <- fread('hesin.txt') 
hesin_diag <- fread('hesin_diag.txt') 
hesin_operations <- fread('hesin_oper.txt') 

#### Function to split grep list into R list ----

split_string <- function(input_string) {
  # Remove '^' and '$' characters from the input string
  cleaned_string <- gsub("[\\^\\$]", "", input_string)
  
  # Use strsplit with regular expressions to split the cleaned string into elements
  elements <- unlist(strsplit(cleaned_string, "\\|"))
  
  # Define a function to split items ending with '[0-9]'
  split_items <- function(item) {
    if (grepl("\\[0-9\\]$", item)) {
      prefix <- sub("\\[0-9\\]$", "", item)
      digits <- 0:9
      new_items <- paste0(prefix, digits)
      return(new_items)
    } else {
      return(item)
    }
  }
  
  # Use lapply to apply the split_items function to each item
  split_elements <- unlist(lapply(elements, split_items))
  
  return(split_elements)
}

#### Function to take first diag dates according to different codes ----

get_first_diag_date <- function(all, diag, operations, icd_codes, opcs4_codes, pheno_name) {
  
  diag_col <- paste0(pheno_name, '_first_diag_date')
  discharge_col <- paste0(pheno_name, '_discharge_date')
  
  first_diag_ICD_df <- all %>%
    select(eid, ins_index, dsource, admidate, disdate, epistart, epiend) %>%
    left_join(diag, by=c('eid', 'ins_index')) %>%
    mutate_if(is.character, list(~na_if(.,""))) %>%
    select(-arr_index, -diag_icd9_nb, -diag_icd10_nb) %>%
    filter(diag_icd10 %in% icd_codes | diag_icd9 %in% icd_codes) %>%
    mutate({{diag_col}} := as.IDate(coalesce(admidate, epistart), format = "%d/%m/%Y")) %>%
    mutate({{discharge_col}} := as.IDate(coalesce(disdate, epiend), format = "%d/%m/%Y")) %>% 
    #Pivot so that we don't have separate columns for icd9 and icd10
    
    dt_pivot_longer(cols=c(diag_icd9, diag_icd10),
                    names_to = 'Diag_Type',
                    values_to = 'Code', 
                    values_drop_na = T) %>%
    filter(!is.na(Code) & !is.na(get(diag_col))) %>%
    select(eid, level, eval(diag_col), eval(discharge_col), Diag_Type, Code) 
  
  print('Finished collecting ICD10 Codes ')
  print(length(unique(first_diag_ICD_df$eid)))
  
  first_diag_OPCS4_df <- all %>%
    select(eid, ins_index, dsource, admidate, disdate, epistart, epiend) %>%
    left_join(operations, by=c('eid', 'ins_index')) %>%
    mutate_if(is.character, list(~na_if(.,""))) %>%
    select(-arr_index, -oper3_nb, -oper4_nb) %>%
    filter(oper3 %in% opcs4_codes | oper4 %in% opcs4_codes) %>%
    mutate(oper3 = as.character(oper3)) %>%
    mutate({{diag_col}} := as.IDate(coalesce(opdate, epistart, admidate), format = "%d/%m/%Y")) %>%
    mutate({{discharge_col}} := as.IDate(NA)) %>% 
    dt_pivot_longer(cols=c(oper3, oper4),
                    names_to = 'Diag_Type',
                    values_to = 'Code', 
                    values_drop_na = T) %>%
    filter(!is.na(Code) & !is.na(get(diag_col))) %>%
    select(eid, level, eval(diag_col), eval(discharge_col), Diag_Type, Code) 
  
  print('Finished collecting OPCS4 Codes')  
  print(length(unique((first_diag_OPCS4_df$eid))))
  
  combined <- rbind(first_diag_ICD_df, first_diag_OPCS4_df) %>%
    group_by(eid) %>%
    filter(get(diag_col) == min(get(diag_col))) %>%
    #filter(get(discharge_col) == max(get(discharge_col))) %>%
    distinct(eid, .keep_all = T) %>%
    ungroup()
  
  print('Finished combining and collecting OPCS4/ICD10 Codes')
  
  return(combined)
}

#### Function to get first diagnosis date after n days from first diag----

get_nth_diag_date_after <- function(all, diag, operations, icd_codes, opcs4_codes, pheno_name, n) {
  
  diag_col <- paste0(pheno_name, '_after_',n,'_days_first_diag_date')
  
  first_diag_ICD_df <- all %>%
    select(eid, ins_index, dsource, admidate, disdate, epistart, epiend) %>%
    left_join(diag, by=c('eid', 'ins_index')) %>%
    mutate_if(is.character, list(~na_if(.,""))) %>%
    select(-arr_index, -diag_icd9_nb, -diag_icd10_nb) %>%
    filter(diag_icd10 %in% icd_codes | diag_icd9 %in% icd_codes) %>%
    mutate({{diag_col}} := as.IDate(coalesce(admidate, epistart), format = "%d/%m/%Y")) %>%
    select(eid, eval(diag_col)) 
  
  print('Finished collecting ICD10 Codes ')
  print(length(unique(first_diag_ICD_df$eid)))
  
  first_diag_OPCS4_df <- all %>%
    select(eid, ins_index, dsource, admidate, disdate, epistart, epiend) %>%
    left_join(operations, by=c('eid', 'ins_index')) %>%
    mutate_if(is.character, list(~na_if(.,""))) %>%
    select(-arr_index, -oper3_nb, -oper4_nb) %>%
    filter(oper3 %in% opcs4_codes | oper4 %in% opcs4_codes) %>%
    mutate(oper3 = as.character(oper3)) %>%
    mutate({{diag_col}} := as.IDate(coalesce(opdate, epistart, admidate), format = "%d/%m/%Y")) %>%
    select(eid, eval(diag_col)) 
  
  print('Finished collecting OPCS4 Codes')  
  print(length(unique((first_diag_OPCS4_df$eid))))
  
  combined <- rbind(first_diag_ICD_df, first_diag_OPCS4_df) %>%
    group_by(eid) %>%
    arrange(get(diag_col)) %>%
    mutate(first_date = min(get(diag_col))) %>%
    ungroup() %>%
    mutate(thirty = as.numeric(as.numeric(difftime(get(diag_col), first_date, units='days')) > eval(n))) %>%
    filter(thirty ==1) %>%
    group_by(eid) %>%
    filter(eval(diag_col) == min(eval(diag_col))) %>%
    select(eid, eval(diag_col)) %>%
    distinct(eid, .keep_all = T) %>%
    ungroup() 
  
  print('Finished combining and collecting OPCS4/ICD10 Codes')
  
  return(combined)
}


# + ICD and OPCS4 Codes ----

hf_codes <- c('I110', 'I130', 'I132','I255','I420','I425','I428',
              'I429','I500','I501','I509','4254','4280','4281','4289')
afib_codes <- c('I480','I481','I482','I483','I484','I489','4273', 'I48')
afib_opcs4 = c('K571','K621','K622','K623','K624','X501','X502')
t1d <- c('E10','E100','E101','E102','E103','E104','E105','E106','E107','E108','E109')
t2d <- c('E11','E110','E111','E112','E113','E114','E115','E116','E117','E118','E119')
chronic_renal_failure <- 'N18[0-9]|I120|I131|I132|585[0-9]'
crf_opcs4 <- 'M01[0-9]'
MI <- 'I21[0-9]|410[0-9]'
Mitrial_Valve_Disorder <- 'I05[0-9]|I080|I081|I34[0-9]|I390|^4240$|3940|3942|3949'
MVD_OPCS4 <- '^K25[0-9]|K301|K341|K351'
Aortic_Valve_Stenosis <- 'I06[0-9]|I080|I082|I083|I35[0-9]|I391|^4241$|^3951$|^3959$'
AVS_OPCS4 <- '^K26[0-9]|^K302'
Ischemic_Heart_Disease <- 'I20[0-9]|I21[0-9]|I22[0-9]|I23[0-9]|I24[0-9]|I25[0-9]|^410[0-9]|^411[0-9]|^412[0-9]|^413[0-9]|^414[0-9]'
COPD <- 'J43|J44|J42|J41|J43[0-9]|J44[0-9]|J42[0-9]|J41[0-9]|^491[0-9]|^492[0-9]|^496[0-9]'
HCM <- c('I421', 'I422', '4251')
DCM_codes <- c('I429', 'I420', '4254')
Cardiac_Arrhythmias <- 'I49[0-9]|I49'
Ventricular_Arrhythmias <- 'I47|I47[0-9]'
ARVC <- c('I428')
CAD <- 'I21[0-9]|I22[0-9]|I23[0-9]|I252' #numbers may not match bc cause of death also used for grouping
CAD_OPSC4 <- '^K49[0-9]|^K50[0-9]|^K75[0-9]|^K40[0-9]|^K41[0-9]|^K42[0-9]|^K43[0-9]|^K44[0-9]|^K45[0-9]|^K46[0-9]'
Unspecific_Stroke <- 'I60|I61|I62|I63|I64|G45|G45[0-9]|I60[0-9]|I61[0-9]|I62[0-9]|I63[0-9]|^430[0-9]|^431[0-9]|^432[0-9]|^433[0-9]|^434[0-9]|^435[0-9]|^436[0-9]'
Ischemic_Stroke <- 'I63|I63[0-9]|^433[0-9]|^434[0-9]'
Hemorrhagic_Stroke <- 'I60|I61|I62|I60[0-9]|I61[0-9]|I62[0-9]|^430[0-9]|^431[0-9]|^432[0-9]'


opcs4_endpoints <- c('K011', 'K012', 'K018', 'K019', 'K021', 'K022', 'K023', 'K024', 
                     'K025', 'K026', 'K028', 'K029')

OPCS4_Defib <- split_string('K59[0-9]')
OPCS4_Pacemaker <- split_string('K60[0-9]|K61[0-9]|K73[0-9]|K74[0-9]')
OPCS4_Ablation <- c('K622', 'K623', 'K621', 'K574')
Pacemaker_ICD10 <- c('Z950', 'Z450')
myomectomy <- c('K245', 'K246', 'K247')

# Create New Phenotypes ICD9+10 ----

afib_df <- get_first_diag_date(hesin_all, hesin_diag, hesin_operations, afib_codes, afib_opcs4, 'Afib')
hf_df <- get_first_diag_date(hesin_all, hesin_diag, hesin_operations, hf_codes, NULL, 'HF')
t2d_df <- get_first_diag_date(hesin_all, hesin_diag, hesin_operations, t2d, NULL, 'T2D')
t1d_df <- get_first_diag_date(hesin_all, hesin_diag, hesin_operations, t1d, NULL, 'T1D')
chronic_renal_failure_df <- get_first_diag_date(hesin_all, hesin_diag, hesin_operations, 
                                                split_string(chronic_renal_failure), split_string(crf_opcs4), 'Chronic_Renal_Failure')
aortic_valve_stenosis_df <- get_first_diag_date(hesin_all, hesin_diag, hesin_operations, split_string(Aortic_Valve_Stenosis), 
                                                split_string(AVS_OPCS4), 'Aortic_Valve_Stenosis')
COPD_df <- get_first_diag_date(hesin_all, hesin_diag, hesin_operations, split_string(COPD), 
                               NULL, 'COPD')
HCM_df <- get_first_diag_date(hesin_all, hesin_diag, hesin_operations, HCM, 
                              NULL, 'HCM')
Myocardial_Infarction_df <- get_first_diag_date(hesin_all, hesin_diag, hesin_operations, split_string(MI), 
                                                NULL, 'MI')
Mitrial_Valve_Disorder_df <- get_first_diag_date(hesin_all, hesin_diag, hesin_operations, split_string(Mitrial_Valve_Disorder), 
                                                 split_string(MVD_OPCS4), 'Mitrial_Valve_Disorder')

DCM_df <- get_first_diag_date(hesin_all, hesin_diag, hesin_operations, DCM_codes, 
                              NULL, 'DCM')

Cardiac_Arrhythmia_df <- get_first_diag_date(hesin_all, hesin_diag,hesin_operations, split_string(Cardiac_Arrhythmias), NULL, 'Cardiac_Arrhythmias')
Ventricular_Arrhythmia_df <- get_first_diag_date(hesin_all, hesin_diag,hesin_operations, split_string(Ventricular_Arrhythmias), NULL, 'Ventricular_Arrhythmias')
Ischemic_HD_df <- get_first_diag_date(hesin_all, hesin_diag,hesin_operations, split_string(Ischemic_Heart_Disease), split_string(CAD_OPSC4), 'Ischemic_Heart_Disease')
ARVC_df <- get_first_diag_date(hesin_all, hesin_diag,hesin_operations, ARVC, NULL, 'ARVC')
Unspecific_Stroke_df <- get_first_diag_date(hesin_all, hesin_diag,hesin_operations, split_string(Unspecific_Stroke), 
                                            NULL, 'Unspecific_Stroke')
Ischemic_Stroke_df <- get_first_diag_date(hesin_all, hesin_diag,hesin_operations, split_string(Ischemic_Stroke), 
                                          NULL, 'Ischemic_Stroke')
Hemorrhagic_Stroke_df <- get_first_diag_date(hesin_all, hesin_diag,hesin_operations, split_string(Hemorrhagic_Stroke), 
                                             NULL, 'Hemorrhagic_Stroke')

# Get nth hospitalization dates after n days----
days <- 30

afib_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag, hesin_operations, afib_codes, afib_opcs4, 'Afib', days)
hf_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag, hesin_operations, hf_codes, NULL, 'HF', days)
t2d_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag, hesin_operations, t2d, NULL, 'T2D', days)
t1d_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag, hesin_operations, t1d, NULL, 'T1D', days)
chronic_renal_failure_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag, hesin_operations, 
                                                    split_string(chronic_renal_failure), split_string(crf_opcs4), 'Chronic_Renal_Failure', days)
aortic_valve_stenosis_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag, hesin_operations, split_string(Aortic_Valve_Stenosis), 
                                                    split_string(AVS_OPCS4), 'Aortic_Valve_Stenosis', days)

COPD_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag, hesin_operations, split_string(COPD), 
                                   NULL, 'COPD', days)
HCM_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag, hesin_operations, HCM, 
                                  NULL, 'HCM', days)
Myocardial_Infarction_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag, hesin_operations, split_string(MI), 
                                                    NULL, 'MI', days)
Mitrial_Valve_Disorder_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag, hesin_operations, split_string(Mitrial_Valve_Disorder), 
                                                     split_string(MVD_OPCS4), 'Mitrial_Valve_Disorder', days)

DCM_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag, hesin_operations, DCM_codes, 
                                  NULL, 'DCM', days)

Cardiac_Arrhythmia_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag,hesin_operations, split_string(Cardiac_Arrhythmias), NULL, 'Cardiac_Arrhythmias', days)
Ventricular_Arrhythmia_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag,hesin_operations, split_string(Ventricular_Arrhythmias), NULL, 'Ventricular_Arrhythmias', days)
Ischemic_HD_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag,hesin_operations, split_string(Ischemic_Heart_Disease), split_string(CAD_OPSC4), 'Ischemic_Heart_Disease', days)
ARVC_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag,hesin_operations, ARVC, NULL, 'ARVC', days)
Unspecific_Stroke_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag,hesin_operations, split_string(Unspecific_Stroke), 
                                                NULL, 'Unspecific_Stroke', days)
Ischemic_Stroke_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag,hesin_operations, split_string(Ischemic_Stroke), 
                                              NULL, 'Ischemic_Stroke', days)
Hemorrhagic_Stroke_n_df <- get_nth_diag_date_after(hesin_all, hesin_diag,hesin_operations, split_string(Hemorrhagic_Stroke), 
                                                 NULL, 'Hemorrhagic_Stroke', days)

#### Merge all first diagnosis dates dataframes into one giant dataframe ----

# Grab all DF names

df_names <- ls(pattern = "_df$")

# Create a list of data frames
df_list <- mget(df_names)

# Merge the data frames by the "eid" column
merged_df <- Reduce(function(x, y) merge(x, y, by = "eid", all = TRUE), df_list) %>%
  #select(-level, -Diag_Type, -Code) %>%
  select(-contains('.'), -contains(c('discharge', 'level', 'Diag_Type', 'Code')))

#### Get new causes of death/cv outcome (SKIP IF UNNEEDED) ---- 

# + Load in new death dates, operation codes, all IDs ----

death_date <- fread('death.txt')
death_cause <- fread('death_cause.txt')

old <- big %>%
  select(IID) %>%
  mutate(study_end_date = as.IDate('2023-09-01'))

outcome_opcs4_codes <- c('K011', 'K012', 'K018', 'K019', 'K021', 'K022', 'K023', 'K024', 
                         'K025', 'K026', 'K028', 'K029')
date_of_birth <- big %>%
  select(IID, birth_date)

# + Data cleaning ----

death_combined <- left_join(death_date, death_cause, by=c('eid', 'ins_index')) %>%
  select(-source, -dsource)

cancer_death <- death_combined %>%
  filter(cause_icd10 %like% 'C') %>%
  pull(eid)

cv_death <- death_combined %>%
  filter(cause_icd10 %like% 'I') %>%
  pull(eid)

just_death_date <- death_combined %>%
  distinct(eid, .keep_all = T) %>%
  select('IID' = eid, date_of_death ) %>%
  mutate(date_of_death = as.IDate(date_of_death, format='%d/%m/%Y'))

outcome_opcs4_df <- get_first_diag_date(hesin_all, hesin_diag, hesin_operations, NULL, outcome_opcs4_codes, 'CV_Outcome' )

cv_operation_date <- outcome_opcs4_df %>%
  select('IID' = eid, 'date_of_cv_operation' = CV_Outcome_first_diag_date)

# + Merge and Create death dataframe ----

all_death_dates <-  Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), list(old, just_death_date, cv_operation_date)) %>%
  mutate(outcome_date = coalesce(date_of_cv_operation, date_of_death, study_end_date)) %>%
  mutate(CV_outcome = ifelse(IID %in% cv_death| !is.na(date_of_cv_operation), 1, 0)) %>%
  mutate(cancer_outcome = ifelse(IID %in% cancer_death, 1, 0)) %>%
  #mutate(all_cause_outcome = ifelse(!(IID %in% cancer_death) & !(IID %in% cv_death) & !is.na(date_of_death), 1, 0))
  mutate(all_cause_outcome = ifelse(!is.na(date_of_death), 1, 0))

#### Calculate time to outcome for each disease ----
combined_first_diag_dates <- merged_df %>%
  rename('IID' = eid)

cv_death_dates <- all_death_dates %>%
  select(IID, outcome_date, CV_outcome) %>%
  left_join(date_of_birth, by='IID') %>%
  left_join(combined_first_diag_dates, by='IID') %>%
  mutate(across(ends_with("first_diag_date"), 
                ~ as.numeric(difftime(outcome_date, .), units = 'days')/365,
                .names = "{strsplit(.col, '_first_diag_date')}_years_to_outcome")) %>%
  mutate(across(ends_with("first_diag_date"), 
                ~ as.numeric(difftime(. , birth_date), units = 'days')/365,
                .names = "age_{strsplit(.col, '_first_diag_date')}_diag")) %>%
  mutate(across(ends_with("first_diag_date"), 
                ~ ifelse(!is.na(.), 1,0),
                .names = "valid_{strsplit(.col, '_first_diag_date')}")) %>%
  select(-outcome_date, -birth_date) 

#### Combine into one huge pheno file 

combined_final <- big %>%
  left_join(cv_death_dates)

first_diag_dates <- combined_final %>%
  select(IID,  contains('first_diag_date')) %>%
  select(-contains('after_30_days'))

thirty_day_rehos <- combined_final %>%
  select(IID, contains('after_30_days'))

before_Isch <- first_diag_dates %>%
  mutate(across(-1, ~ as.numeric(. <= Ischemic_Heart_Disease_first_diag_date), .names = "{.col}_before_Ischemic_Heart_Disease")) %>%
  select(IID, contains('before')) %>%
  mutate_all(~replace_na(., 0))

before_AF <- first_diag_dates %>%
  mutate(across(-1, ~ as.numeric(. <= Afib_first_diag_date), .names = "{.col}_before_AF")) %>%
  select(IID, contains('before')) %>%
  mutate_all(~replace_na(., 0))

before_HF <- first_diag_dates %>%
  mutate(across(-1, ~ as.numeric(. <= HF_first_diag_date), .names = "{.col}_before_HF")) %>%
  select(IID, contains('before')) %>%
  mutate_all(~replace_na(., 0))

pheno_with_diag <- combined_final %>%
  select(-contains('after_30_days'))

#Writing phenofiles for use in survival 

try(dir.create(path = './phenofiles/'))

fwrite(pheno_with_diag, 
       file='./phenofiles/pheno_with_diagnosis_dates.txt', 
       sep = '\t', 
       quote=F, 
       na=NA,
       row.names = F)

fwrite(thirty_day_rehos, 
       file='./phenofiles/thirty_day_rehos.txt', 
       sep = '\t', 
       quote=F, 
       na=NA,
       row.names = F)

fwrite(before_AF, 
       file='./phenofiles/before_AF.txt', 
       sep = '\t', 
       quote=F, 
       na=NA,
       row.names = F)

fwrite(before_HF, 
       file='./phenofiles/before_HF.txt', 
       sep = '\t', 
       quote=F, 
       na=NA,
       row.names = F)

fwrite(before_Isch, 
       file='./phenofiles/before_IHD.txt', 
       sep = '\t', 
       quote=F, 
       na=NA,
       row.names = F)

fwrite(all_death_dates, 
       file='./phenofiles/cv_outcome_raw.txt', 
       sep = '\t', 
       quote=F, 
       na=NA,
       row.names = F)

#### QC ---- 

check <- fread('~/ukbiobank/phenotypes/All_Disease_HESIN_April25_2024.txt') %>%
  left_join(fread('~/ukbiobank/bulk/record-level-HES/rehospitalizations_after_thirty_days_UKB.txt')) %>%
  select(IID, CV_outcome, contains(c('after_30_days', 'valid_', 'years_to_outcome', '_diag'))) %>%
  select(-contains(c('CAD', 'myomectomy')))
