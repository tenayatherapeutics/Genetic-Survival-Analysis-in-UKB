.libPaths('/home/lzhang@TENAYA.local/my_bioinfo/rpackages/')

library(tidyverse)
library(data.table)

setwd('/home/lzhang@TENAYA.local/ukbiobank/Paper_Scripts/')

#### Basic phenofile and requisite data preparations ----

MedGroup <- fread('MedicationGroups.txt')
codings <- fread('coding4.tsv')
result <- fread('merged_med_groups_stringent_anticoag.txt') #Result from running the medications.py script 
base_pheno <- fread('example_survival_base_pheno.txt') 

medcols <- base_pheno %>%
  select(starts_with('Med')) %>%
  colnames()

countscols <- base_pheno %>%
  select(starts_with('RX')) %>%
  colnames()

All_Med <- unite(base_pheno, col='Medication_List', all_of(medcols), sep=',', na.rm=T) %>%
  select(IID, Medication_List)

Med_Counts <- unite(base_pheno, col='Medication_Count', all_of(countscols), sep=',', na.rm=T) %>%
  select(IID, Medication_Count) %>%
  filter(Medication_Count != '')

#### Creating Drug Dataframe ----

pivoted <- result %>%
  pivot_wider(names_from = Drug, values_from = Coding)

cols <- pivoted %>%
  select(-Group) %>%
  colnames()

final <- unite(pivoted, col='Drug_Numbers', all_of(cols), sep=',', na.rm=T)

patient_all <- All_Med %>%
  mutate(med = str_split(Medication_List, ','))

#Turns groups into vector of drug codes, check to see if each patient has medication codes within the groups  
for (grouping in final$Group) {
  print(grouping)
  assign(grouping, unlist(str_split(subset(final, Group == grouping)$Drug_Numbers, ',')))
  patient_all[[grouping]] <- sapply(patient_all$med, function(x) ifelse(any(x %in% get(grouping)), 1, 0))
}


patient_final <- patient_all %>%
  select(-med) 

patient_no_list <- patient_final %>%
  select(-Medication_List, `Anti-Coagulants`) %>%
  filter(IID %in% Med_Counts$IID) 

#### Repeat process for anti-coagulants ---- 

anticoag <- result %>%
  filter(Group == 'Anti-Coagulants') %>%
  filter(Drug %like% 'phenindione|warfarin|marevan|plavix|clopidogrel|aspirin 75mg')

pivoted <- anticoag %>%
  pivot_wider(names_from = Drug, values_from = Coding)

cols <- pivoted %>%
  select(-Group) %>%
  colnames()

final <- unite(pivoted, col='Drug_Numbers', all_of(cols), sep=',', na.rm=T)

patient_all <- All_Med %>%
  mutate(med = str_split(Medication_List, ','))

#Turns groups into vector of drug codes, check to see if each patient has medication codes within the groups  
for (grouping in final$Group) {
  print(grouping)
  assign(grouping, unlist(str_split(subset(final, Group == grouping)$Drug_Numbers, ',')))
  patient_all[[grouping]] <- sapply(patient_all$med, function(x) ifelse(any(x %in% get(grouping)), 1, 0))
}


patient_final <- patient_all %>%
  select(-med) 

anti_coag_list <- patient_final %>%
  select(-Medication_List) %>%
  filter(IID %in% Med_Counts$IID) 

#### Recombine anti-coag and add medications to big pheno file to be used in further analysis 

final_pheno_file <- base_pheno %>%
  select(-starts_with(c('Med', 'RX'))) %>%
  left_join(patient_no_list) %>%
  left_join(anti_coag_list)

just_medi <- patient_no_list %>%
  left_join(anti_coag_list)

fwrite(final_pheno_file, 'example_pheno_with_meds.txt', sep='\t', na=NA, quote=F, row.names = F)
#fwrite(just_medi, '~/ukbiobank/Medication_Groups/new_medications_May172024.txt', sep = '\t', na=NA, quote=F, row.names = F)

