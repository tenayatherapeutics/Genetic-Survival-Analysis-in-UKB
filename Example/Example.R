.libPaths('/home/lzhang@TENAYA.local/my_bioinfo/rpackages')

setwd('./Example/')
library(tidyverse)
library(data.table)
library(survminer)
library(survival)
library(MatchIt)  
library(expss)
library(cobalt)

#### Load "master pheno" file ----

set.seed(1)
raw_master_pheno <- fread('example_pheno_with_AF_rehospitalization.txt') 

plotsurv <- function(rsnum, gene_name, pheno_df, phenotype, endpoint, cut_time, ylim, adj=FALSE) {
  
  death_dates <- fread('example_CV_outcome') %>%
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

#### Example Dosage Input ----

example_dosage <- fread('example_dosage.txt') %>%
  mutate_if(is.numeric, round)

gn <- 'Example'
rsnum <- 'SNP1'


example_matching_terms <- c('flipped',rsnum, 'raw', 'age', 'CV_outcome',
                          'valid_Ischemic_Heart_Disease', 'Ezetimibe', 'valid_T2D', 
                          'Hemorrhagic_Stroke_first_diag_date_before_Ischemic_Heart_Disease',
                          'Ischemic_Stroke_first_diag_date_before_Ischemic_Heart_Disease',
                          'Diabetic', 'LDL_1', 
                          'Centrally_acting_antihypertensive',
                          'Vasodilator_antihypertensive_drugs')


target_df <- example_dosage %>%
  select(IID, eval(rsnum)) %>%
  mutate(raw = get(rsnum)) %>%
  mutate_if(is.numeric, round) %>%
  mutate({{rsnum}} := as.numeric( get(rsnum) > 0)) %>%
  mutate(flipped = abs(get(rsnum) - 1))

pheno1 <- 'Afib'

before_Afib <- fread('example_before.txt') %>%
  select(-contains(c('Afib_first', 'CAD', 'ARVC', 'Hemorrhagic', 'MI_first')))

master_pheno <- raw_master_pheno %>%
  left_join(get(paste0('before_', eval(pheno1))))

pheno_filtered <- master_pheno %>%
  left_join(target_df) %>%
  filter(get(eval(paste0('valid_', pheno1))) == 1) %>%
  filter(get(eval(paste0( pheno1, '_years_to_outcome'))) > 0) 

age_pheno_diag <- pheno_filtered %>%
  select(IID, (eval(paste0('age_', pheno1, '_diag'))))

matched_pheno_filtered <- matchit_all(subset_df = pheno_filtered,
                                      list_of_conditions = list(target_df, age_pheno_diag), 
                                      phenotype = NA,
                                      SNP = rsnum, 
                                      omit_cols = c(example_matching_terms))

years_to_outcome <- pheno_filtered %>%
  select(IID, eval(paste0( pheno1, '_years_to_outcome')))


matched_with_diag <- matched_pheno_filtered %>%
  left_join(years_to_outcome)

#### + Unmatched 

plotsurv('raw', gn, pheno_filtered, eval(pheno1), 'CV_outcome', 10, 0.0, adj=TRUE)

fit <- survfit(Surv(years_to_outcome_20, outcome_20)~get('raw'), data=look)
#summary(fit)

unmatched <- ggsurvplot(fit,
           pval = F, 
           pval.coord = c(0, .55),
           conf.int = FALSE,
           risk.table = TRUE,
          # fontsize = 8,# Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           ggtheme = theme_classic2(), # Change ggplot2 theme
           title=paste(''),
           xlab = 'Time to Outcome',
           legend.title=paste(rsnum, 'Genotypes'),
           legend='none',
           #legend.labs = c('GG', 'GT/TT'),
           legend.labs = c('XX', 'XY', 'YY'),
           ylim=c(0.2,1))

unmatched 

#### + Matched 

plotsurv('raw', gn, matched_with_diag, eval(pheno1), 'CV_outcome', 25, 0.0, adj=FALSE)

fit <- survfit(Surv(years_to_outcome_20, outcome_20)~get('raw'), data=look)
#summary(fit)

plotty <- ggsurvplot(fit,
           pval = F, 
           pval.coord = c(0, .7),
           conf.int = FALSE,
           risk.table = TRUE,
           #fontsize = 8,# Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           ggtheme = theme_classic2(), # Change ggplot2 theme
           title='',
           xlab = 'Time',
           legend.title=paste(rsnum, 'alleles'),
           legend='none',
           #legend.labs = c('GG', 'GT/TT'),
           legend.labs = c('XX', 'XY', 'YY'), #This is for rs562556
           ylim=c(0.6,1))

plotty 

#### Example file changes ---- 

new_pheno <- raw_master_pheno %>%
  #mutate(valid_Afib = 1, valid_Afib_after_30_days = 1)
  mutate(Afib_years_to_outcome = runif(30, 0.1, 10))

fwrite(new_pheno, file = 'example_pheno_with_AF_rehospitalization.txt', quote=F, na=NA)

example_dosage$SNP1 <- sample(c(rep(0, 16), sample(c(1,2), 14, replace=T)))

fwrite(example_dosage, file='example_dosage.txt', quote=F, na=NA)

### Power Calc 

library(survSNP)
GHRs<-seq(.78, .95, by=0.01)
ns<-c(2640, 48905)
rafs<-c(0.25, 0.015)
erates<-c(.12, .12)

res2<-survSNP.power.table(GHRs,ns,rafs,erates,
                          pilm=0.5,lm=1,
                          model="additive",test="additive",
                          alpha=0.000005)
res2$power=round(res2$pow0*100,2)


#Unmatched - PCSK9 rs1--
res1<-sim.snp.expsurv.power(.95, n=48905, raf=0.015, erate=0.12,
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]

#Matched - PCSK9 rs1--
res1<-sim.snp.expsurv.power(0.78, n=2640, raf=0.25, erate=0.12,
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]

#Unmatched - PCSK9 rs5--
res1<-sim.snp.expsurv.power(1.01, n=48906, raf=0.82, erate=0.12,
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]

#Matched - PCSK9 rs5--
res1<-sim.snp.expsurv.power(1.06, n=2778, raf=0.43, erate=0.11,
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]
