.libPaths('/home/lzhang@TENAYA.local/my_bioinfo/rpackages')

setwd('~/ukbiobank/Paper_Scripts/datafiles/HES/phenofiles/')
library(tidyverse)
library(data.table)
library(survminer)
library(survival)
library(MatchIt)  
library(expss)

#### Load master pheno ----
set.seed(1)
raw_master_pheno <- fread('pheno_with_diagnosis_dates.txt') 

plotsurv <- function(rsnum, gene_name, pheno_df, phenotype, endpoint, cut_time, ylim, adj=FALSE) {
  
  death_dates <- fread('cv_outcome_raw.txt') %>%
    select(IID, date_of_cv_operation, date_of_death, contains('study')) %>%
    #mutate(end_date = coalesce(date_of_death, study_end_date)) 
    mutate(end_date = coalesce(date_of_cv_operation,date_of_death, study_end_date))
  
  pheno_df <- pheno_df %>%
    left_join(death_dates)
  
  var <- eval(rsnum) #Plot for combined heterozygous + homozygous
  #print(var)
  age_diag <- paste0('age_', phenotype, '_diag')
  #print(age_diag)
  
  if((endpoint %like% 'outcome')) {
    years_to_outcome <- paste0(phenotype, '_years_to_outcome')
    #print(years_to_outcome)
    
    pheno_df$outcome_20<- ifelse(pheno_df[[eval(years_to_outcome)]] > cut_time, 0, pheno_df[[eval(endpoint)]])
    pheno_df$years_to_outcome_20 <- ifelse(pheno_df[[eval(years_to_outcome)]] > cut_time, cut_time, pheno_df[[eval(years_to_outcome)]])
    #Exclude population of people with years to outcome < 0
    #pheno_df$years_to_outcome_20 <- ifelse(pheno_df[[eval('years_to_outcome_20')]] < 0, 0, pheno_df[[eval('years_to_outcome_20')]])
    print(summary(pheno_df$years_to_outcome_20))
    look <<- pheno_df
    if(adj == TRUE) {
      coxformula <- formula(paste0('Surv(years_to_outcome_20, outcome_20)~', var,'+sex+', eval(age_diag)))
      #coxformula <- formula(paste0('Surv(years_to_outcome_20, outcome_20)~', var,'+sex'))
    } else {
      coxformula <- formula(paste0('Surv(years_to_outcome_20, outcome_20)~', var))
    }
    
    #coxformula <- formula(paste0('Surv(years_to_outcome_20, outcome_20)~', 'sex'))
    
    
    pheno_df <- apply_labels(pheno_df, outcome_20 = 'Outcome')
    print(knitr::kable(cross_cases(data=pheno_df, get('outcome_20'), get(var)), caption=paste(gene_name, rsnum, 'Crosstabs')))
    
  } else {
    first_pheno_diag <- paste0(phenotype, '_first_diag_date')
    print(first_pheno_diag)
    first_endpoint_diag <- paste0(endpoint, '_first_diag_date')
    print(first_endpoint_diag)
    pheno_df$outcome <- as.numeric(!is.na(pheno_df[[eval(first_endpoint_diag)]]) | pheno_df[[eval('CV_outcome')]] == 1)
    #pheno_df$outcome <- as.numeric(!is.na(pheno_df[[eval(first_endpoint_diag)]]) )
    pheno_df$outcome_date <- pmin(pheno_df[[eval(first_endpoint_diag)]], pheno_df[[eval('end_date')]], na.rm=T)
    pheno_df$years_to_outcome <- as.numeric(difftime(pheno_df[[eval('outcome_date')]], pheno_df[[first_pheno_diag]]), units = 'days')/365
    pheno_df$outcome_20<- ifelse(pheno_df[['years_to_outcome']] > cut_time, 0, pheno_df[['outcome']])
    pheno_df$years_to_outcome_20 <- ifelse(pheno_df[[eval('years_to_outcome')]] > cut_time, cut_time, pheno_df[[eval('years_to_outcome')]])
    print(summary(pheno_df$years_to_outcome_20))
    
    look <<- pheno_df
    if(adj == TRUE) {
      coxformula <- formula(paste0('Surv(years_to_outcome_20, outcome_20)~', var,'+sex+', eval(age_diag)))
      #coxformula <- formula(paste0('Surv(years_to_outcome_20, outcome_20)~', var,'+sex'))
    } else {
      coxformula <- formula(paste0('Surv(years_to_outcome_20, outcome_20)~', var))
    }
    #
    pheno_df <- apply_labels(pheno_df, outcome_20 = 'Outcome')
    print(knitr::kable(cross_cases(data=pheno_df, get('outcome_20'), get(var)), caption=paste(gene_name, rsnum, 'Crosstabs')))
    
  }
  
  print(noquote('###########################################################################'))
  print(noquote('###########################################################################'))
  print(noquote('###########################################################################'))
  
  print(coxformula)
  
  fit <- coxph(coxformula,pheno_df)
  print(summary(fit))
  
  plot <- ggadjustedcurves(fit,
                           method = "average",
                           variable = eval(var),
                           #variable = 'subclass',
                           data = pheno_df,
                           ggtheme = theme_bw(), # Change ggplot2 theme
                           #title=paste("Effect of", gn, rsnum, "mutation on Afib Patients"),
                           title=paste("Effect of", gene_name, rsnum, "on time to", endpoint, "from", phenotype, "diagnosis"),
                           xlab = paste('Years to', endpoint),
                           ylim = c(ylim,1))
  
  #print(plot)
  
}

#### Create Matchit Function ----

matchit_all <- function(subset_df, list_of_conditions, phenotype, SNP, omit_cols) { 
  #Basic Variable settings
  basic <- c('IID', 'sex', 'CV_outcome', 'CReactive_protein', 'smoking0', 
             'BMI', 'LDL_1', 'first_systolic_bp', 
             'Diuretics', 'Beta-blockers', 'Nitrates', 'ACEi', 'Calcium-Channel Blockers', 
             'Vasodilator antihypertensive drugs', 'Centrally acting antihypertensive', 'Alpha-Adrenoceptor Blocking Drugs', 'Statins', 
             'Digoxin', 'ARB', 'Diabetic', 'Ezetimibe', 'Anti-Coagulants', 'age')
  
  years_to_outcome <- paste0(phenotype,'_years_to_outcome')
  age_of_diag <- paste0('age_',phenotype, '_diag')
  
  if(is.na(phenotype)) {
    combined <- basic
  } else {
    combined <- c(basic, years_to_outcome, age_of_diag)
  }
  
  
  
  #Join all df and filter 
  selected_cols <- subset_df %>%
    select(all_of(combined), contains('before')) 
  
  selected_phenos <<- reduce(list_of_conditions, ~left_join(.x, .y, by='IID'), .init = selected_cols) %>%
    na.omit()
  
  names(selected_phenos) <- gsub(" ", "_", names(selected_phenos))
  names(selected_phenos) <- gsub("-", "_", names(selected_phenos))
  
  print(colnames(selected_phenos))
  
  #Filter out columns we don't want to match with 
  filtered_cols <- colnames(selected_phenos) [! colnames(selected_phenos) %in% c('IID', 'CV_outcome', years_to_outcome, age_of_diag
                                                                                 ,'Nitrates', omit_cols)]
  
  Formula <- formula(paste(paste0(SNP, '~'), eval(paste(filtered_cols, sep='+', collapse='+'))))
  print(Formula)
  
  inbalance <<- matchit(Formula, data=selected_phenos, method = 'nearest', distance = 'glm', link = 'probit', ratio=1)
  
  print(summary(inbalance))
  
  print(plot(summary(inbalance)))
  
  matched <- match.data(inbalance)
  return(matched)
  
}

# + ACE for AFib ---- 

ACE_dosage <- fread('~/ukbiobank/dosage/ACE1_Genecards_dosage.txt') %>%
  select(-FID) %>%
  mutate_if(is.numeric, round) 

rsnum <- "rs4968782" 
rsnum <- 'rs4363'
gn <- 'ACE'

target_df <- ACE_dosage %>%
  select(IID, eval(rsnum)) %>%
  mutate(raw = get(rsnum)) %>%
  mutate_if(is.numeric, round) %>%
  mutate({{rsnum}} := as.numeric( get(rsnum) > 0)) %>%
  mutate(flipped = abs(get(rsnum) - 1))

ACE_AF_matching <- c('flipped',rsnum, 'raw', 'age', 'CV_outcome',
                     'valid_Afib', 
                     'Ezetimibe',
                     #'Diabetic',
                     'ACEi', 'Vasodilator_antihypertensive_drugs',
                     'Centrally_acting_antihypertensive',
                     'DCM_first_diag_date_before_AF',
                     'HCM_first_diag_date_before_AF',
                     'T1D_first_diag_date_before_AF',
                     'ARVC_first_diag_date_before_AF',
                     'Hemorrhagic_Stroke_first_diag_date_before_AF',
                     'Ischemic_Stroke_first_diag_date_before_AF',
                     'first_systolic_bp',
                     'Mitrial_Valve_Disorder_first_diag_date_before_AF'
                     #'Ventricular_Arrhythmias_first_diag_date_before_AF',
                     #'Cardiac_Arrhythmias_first_diag_date_before_AF'
)

pheno1 <- 'Afib'

before_Afib <- fread('before_AF.txt') %>%
  select(-contains(c('Afib_first', 'CAD', 'MI_first', 'Myomectomy')))

master_pheno <- raw_master_pheno %>%
  left_join(get(paste0('before_', eval(pheno1))))

pheno_filtered <- master_pheno %>%
  #filter(kinship_to_other_participants == 0 ) %>%
  filter(ethnic_group == 1) %>%
  filter(ACEi == 0 & `Vasodilator antihypertensive drugs` ==0 & `Centrally acting antihypertensive` == 0 ) %>%
  left_join(target_df) %>%
  filter(get(eval(paste0('valid_', pheno1))) == 1) %>%
  filter(get(eval(paste0( pheno1, '_years_to_outcome'))) > 0) 

age_pheno_diag <- pheno_filtered %>%
  select(IID, (eval(paste0('age_', pheno1, '_diag'))))

matched_pheno_filtered <- matchit_all(subset_df = pheno_filtered,
                                      list_of_conditions = list(target_df, age_pheno_diag), 
                                      phenotype = NA,
                                      SNP = 'flipped', 
                                      omit_cols = c(ACE_AF_matching))

years_to_outcome <- pheno_filtered %>%
  select(IID, eval(paste0( pheno1, '_years_to_outcome')))


matched_with_diag <- matched_pheno_filtered %>%
  left_join(years_to_outcome)

#### Plot Loveplot 

tiff('/home/lzhang@TENAYA.local/ukbiobank/Paper_Figures/ACE_4363_AF_CV_Outcome_LovePlot_WB.tiff', 
     width = 7, height = 6, units = 'in', res = 320)  #width=7, height=6)
love.plot((inbalance),abs = TRUE, stars='raw') +
  theme(legend.box.background = element_rect(), 
        legend.box.margin = margin(1, 1, 1, 1))
dev.off()

#### + Unmatched 

plotsurv('raw', gn, pheno_filtered, eval(pheno1), 'CV_outcome', 25, 0.0, adj=TRUE)

fit <- survfit(Surv(years_to_outcome_20, outcome_20)~get('raw'), data=look)
#summary(fit)

unmatched <- ggsurvplot(fit,
                        pval = F, 
                        pval.coord = c(0, .85),
                        conf.int = FALSE,
                        risk.table = TRUE,
                        #fontsize = 8,# Add risk table
                        risk.table.col = "strata", # Change risk table color by groups,
                        #risk.table.title = '',
                        linetype = "strata", # Change line type by groups
                        #surv.median.line = "hv", # Specify median survival
                        ggtheme = theme_classic2(), # Change ggplot2 theme
                        title='',
                        xlab = 'Time to Outcome',
                        legend.title=paste(rsnum, 'Genotypes'),
                        legend='none',
                        legend.labs = c('GG','GA', 'AA'),
                        ylim=c(0.6,1))

unmatched 

##save plot as a high quality pdf:
tiff('/home/lzhang@TENAYA.local/ukbiobank/Paper_Figures/ACE_4363_AF_to_CV_Outcome_UnMatched_WB.tiff', width = 7, height = 6, units = 'in', res = 320)  #width=7, height=6)
print(unmatched)
dev.off()

#### + Matched 

plotsurv('raw', gn, matched_with_diag, eval(pheno1), 'CV_outcome', 25, 0.0, adj=FALSE)

fit <- survfit(Surv(years_to_outcome_20, outcome_20)~get('raw'), data=look)
#summary(fit)

plotty <- ggsurvplot(fit,
                     pval = F, 
                     pval.coord = c(0, .85),
                     conf.int = FALSE,
                     risk.table = TRUE,
                     #fontsize = 8,# Add risk table
                     risk.table.col = "strata", # Change risk table color by groups,
                     #risk.table.title = '',
                     linetype = "strata", # Change line type by groups
                     #surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_classic2(), # Change ggplot2 theme
                     title='',
                     xlab = 'Time to Outcome',
                     legend.title=paste(rsnum, 'Genotypes'),
                     legend='none',
                     legend.labs = c('GG','GA', 'AA'),
                     ylim=c(0.6,1))

print(plotty)

##save plot as a high quality pdf:
tiff('/home/lzhang@TENAYA.local/ukbiobank/Paper_Figures/ACE_4363_AF_to_CV_Outcome_Matched_WB.tiff', width = 7, height = 6, units = 'in', res = 320)  #width=7, height=6)
print(plotty)
dev.off()


#####################################################################################################
#HF REHOS ----
#####################################################################################################

ACE_more_curated_exclusion <- c('flipped',rsnum, 'raw', 'age', 'CV_outcome',
                                'valid_Afib', 
                                'Ezetimibe',
                                #'Diabetic',
                                'ACEi', 'Vasodilator_antihypertensive_drugs',
                                'Centrally_acting_antihypertensive',
                                'DCM_first_diag_date_before_HF',
                                'HCM_first_diag_date_before_HF',
                                'T1D_first_diag_date_before_HF',
                                'ARVC_first_diag_date_before_HF',
                                'Hemorrhagic_Stroke_first_diag_date_before_HF',
                                'Ischemic_Stroke_first_diag_date_before_HF',
                                'first_systolic_bp',
                                'Beta_blockers', 
                                'Alpha_Adrenoceptor_Blocking_Drugs',
                                'Mitrial_Valve_Disorder_first_diag_date_before_HF',
                                'Ventricular_Arrhythmias_first_diag_date_before_HF',
                                'Cardiac_Arrhythmias_first_diag_date_before_HF')

pheno1 <- 'HF'
pheno2 <- 'HF_after_30_days'

second_hosp <- fread('thirty_day_rehos.txt') %>%
  select(contains(c(eval(pheno1), 'IID'))) 

before_HF<- fread('before_HF.txt') %>%
  select(-contains(c('HF_first', 'CAD', 'MI_first')))

master_pheno <- raw_master_pheno %>%
  left_join(second_hosp) %>%
  left_join(get(paste0('before_', eval(pheno1))))

pheno_filtered <- master_pheno %>%
  filter(kinship_to_other_participants == 0) %>%
  #filter(ethnic_group == 1) %>%
  filter(ACEi == 0 & `Beta-blockers` ==0 & 
           `Centrally acting antihypertensive` == 0 & `Alpha-Adrenoceptor Blocking Drugs` ==0 &
           `Vasodilator antihypertensive drugs` ==0) %>%
  left_join(target_df) %>%
  filter(get(eval(paste0('valid_', pheno1))) == 1) %>%
  filter(get(eval(paste0( pheno1, '_years_to_outcome'))) > 0) %>%
  filter(get(eval(paste0('age_', pheno1, '_diag'))) - get(eval(paste0('age_', pheno2, '_diag'))) < 0
         | is.na(get(eval(paste0('age_', pheno1, '_diag'))) - get(eval(paste0('age_', pheno2, '_diag'))))) 

age_pheno_diag <- pheno_filtered %>%
  select(IID, (eval(paste0('age_', pheno1, '_diag'))))

first_diags <- master_pheno %>%
  select(IID, (eval(paste0(pheno1, '_first_diag_date'))), (eval(paste0(pheno2, '_first_diag_date'))))

matched_pheno_filtered <- matchit_all(subset_df = pheno_filtered,
                                      list_of_conditions = list(target_df, age_pheno_diag), 
                                      phenotype = NA,
                                      SNP = 'flipped', 
                                      omit_cols = c(ACE_more_curated_exclusion))

matched_with_diag <- matched_pheno_filtered %>%
  left_join(first_diags)

#### Plot Loveplot 

tiff('/home/lzhang@TENAYA.local/ukbiobank/Paper_Figures/ACE_4363_HF_HF_LovePlot_WB.tiff', 
     width = 7, height = 6, units = 'in', res = 320)  #width=7, height=6)
love.plot((inbalance),abs = TRUE, stars='raw') +
  theme(legend.box.background = element_rect(), 
        legend.box.margin = margin(1, 1, 1, 1))
dev.off()

#### + Unmatched 

plotsurv('raw', gn, pheno_filtered, eval(pheno1), eval(pheno2), 2, 0.0, adj=TRUE)

fit <- survfit(Surv(years_to_outcome_20, outcome_20)~get(rsnum), data=look)
#summary(fit)

unmatched_rehos <- ggsurvplot(fit,
                              pval = F, 
                              pval.coord = c(0, .85),
                              conf.int = FALSE,
                              risk.table = TRUE,
                              #fontsize = 8,# Add risk table
                              risk.table.col = "strata", # Change risk table color by groups,
                              #risk.table.title = '',
                              linetype = "strata", # Change line type by groups
                              #surv.median.line = "hv", # Specify median survival
                              ggtheme = theme_classic2(), # Change ggplot2 theme
                              title='',
                              xlab = 'Time to Outcome',
                              legend.title=paste(rsnum, 'Genotypes'),
                              legend='none',
                              legend.labs = c('GG', 'GA/AA'),
                              ylim=c(0.5,1))

unmatched_rehos 

##save plot as a high quality pdf:
tiff('/home/lzhang@TENAYA.local/ukbiobank/Paper_Figures/ACE_4363_HF_HF_UnMatched_WB.tiff', width = 7, height = 6, units = 'in', res = 320)  #width=7, height=6)
print(unmatched_rehos)
dev.off()

#### + Matched 

plotsurv('raw', gn, matched_with_diag, eval(pheno1), eval(pheno2), 2, 0.0, adj=FALSE)

fit <- survfit(Surv(years_to_outcome_20, outcome_20)~get(rsnum), data=look)
#summary(fit)

plotty_rehos <- ggsurvplot(fit,
                           pval = F, 
                           pval.coord = c(0, .85),
                           conf.int = FALSE,
                           risk.table = TRUE,
                           #fontsize = 8,# Add risk table
                           risk.table.col = "strata", # Change risk table color by groups,
                           #risk.table.title = '',
                           linetype = "strata", # Change line type by groups
                           #surv.median.line = "hv", # Specify median survival
                           ggtheme = theme_classic2(), # Change ggplot2 theme
                           title='',
                           xlab = 'Time to Outcome',
                           legend.title=paste(rsnum, 'Genotypes'),
                           legend='none',
                           legend.labs = c('GG', 'GA/AA'),
                           ylim=c(0.5,1))

print(plotty_rehos)

##save plot as a high quality pdf:
tiff('/home/lzhang@TENAYA.local/ukbiobank/Paper_Figures/ACE_HF_HF_Matched_WB.tiff', width = 7, height = 6, units = 'in', res = 320)  #width=7, height=6)
print(plotty_rehos)
dev.off()


### Power Calc ----

library(survSNP)

# + Unmatched - ACE rs496 AF to CVD----
res1<-sim.snp.expsurv.power(.94, n=23823, raf=0.613, erate=0.11,
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]

# + Matched - ACE rs496 AF to CVD----
res1<-sim.snp.expsurv.power(.85, n=6400, raf=0.36, erate=0.105,
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]

# + Unmatched - ACE rs496 HF to HF----
res1<-sim.snp.expsurv.power(.93, n=7952, raf=0.61, erate=0.39,
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]

# + matched - ACE rs496 HF to HF----
res1<-sim.snp.expsurv.power(.84, n=2168, raf=0.35, erate=0.39,
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]

#####
#--------------------------------------------------------------------------
#####

# + Unmatched - ACE rs4363 AF to CVD----
res1<-sim.snp.expsurv.power(.97, n=24331, raf=0.48, erate=0.11,
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]

# + Matched - ACE rs4363 AF to CVD----
res1<-sim.snp.expsurv.power(.95, n=11904, raf=0.33, erate=0.11,
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]

# + Unmatched - ACE rs4363 HF to HF----
res1<-sim.snp.expsurv.power(1.006, n=8124, raf=0.48, erate=0.39,
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]

# + matched - ACE rs4363 HF to HF----
res1<-sim.snp.expsurv.power(1.04, n=3970, raf=0.33, erate=0.39,
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]


