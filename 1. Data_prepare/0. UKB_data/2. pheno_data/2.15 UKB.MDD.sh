## Step1: disease info
library(stringr)
library(data.table)
library(dplyr)
library(magrittr)
source('/gpfs/gibbs/pi/zhao/cz354/UKB/pheno/ukbcase_v2.R')
source('/gpfs/gibbs/pi/zhao/cz354/UKB/pheno/function.R')

disease="MDD"
MDD_icd10 <- c(paste0('F32',c(0:3,8,9)),
               paste0('F33',c(0:4,8,9)))
MDD_icd9 <- NULL
MDD_ops <- NULL
MDD_s1 <- c(11) # Field 20544
MDD_s2 <- c(1286) # Field 20002

# read data
hesin=fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/ukb_hesin.tsv")
hesin_diag9=fread('/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/ukb_hesin_diag9.tsv')
hesin_diag10=fread('/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/ukb_hesin_diag10.tsv')
hesin_oper=fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/hesin_oper.txt")

death_data=fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/death_data.txt")

selfreport_20544 <- fread('/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/self_report_20544.tsv')
selfreport <- fread('/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/self_report.csv')
disease_age <- fread('/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/disease_age.csv',header=T)
cd1=fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/cd1.txt")

# from self-report, here disease_age1/2 is age corresponding to the disease in th selfreport1/2
field_selfreport1 <- c("f.eid",colnames(selfreport_20544)[!is.na(str_extract(colnames(selfreport_20544),'20544'))])
field_selfreport2 <- c('f.eid',colnames(selfreport)[!is.na(str_extract(colnames(selfreport),'20002'))])
field_diseaseage1 <- c('f.eid',colnames(disease_age)[!is.na(str_extract(colnames(disease_age),'20009'))])
field_diseaseage2 <- c('f.eid',colnames(disease_age)[!is.na(str_extract(colnames(disease_age),'20011'))])
selfreport1 <- selfreport_20544[,..field_selfreport1]
selfreport2 <- selfreport[,..field_selfreport2]
disease_age1 <- disease_age[,..field_diseaseage1]
disease_age2 <- disease_age[,..field_diseaseage2]


# Disease code in the hospital records
disease_icd10 <- get(paste0(disease,'_icd10'))
disease_icd9 <- get(paste0(disease,'_icd9'))
disease_ops <- get(paste0(disease,'_ops'))
disease_s1 <- get(paste0(disease,'_s1'))
disease_s2 <- get(paste0(disease,'_s2'))
disease_case <- ukbcase(hesin = hesin, hesin_diag10 = hesin_diag10, hesin_diag9= hesin_diag9, 
                        hesin_oper4 = hesin_oper, deathdata = death_data, icd10 = disease_icd10, 
                        icd9 = disease_icd9,oper4 = disease_ops)

### The death description in disease_case only indicate that individual's epistart time is also death time,
### It doesn't mean that there are other individuals in this case don't die.
disease_case <- age_calculation(disease_case)

# Disease code in the self-report data
Disease_res_base <- idex_case_report(selfreport1,selfreport2,disease_age1,disease_age2,disease_s1,disease_s2)
self_age=Disease_res_base[[2]]
s=selfreport1[,"f.eid"]
s$self_disease=0
s$self_disease[Disease_res_base[[1]]]=1
s$self_age=Disease_res_base[[2]]
colnames(s)[1]="eid"

# summary Disease case  
case_summary <- data.frame(eid=cd1$eid, sex=cd1$sex, age_recruit=cd1$age_recruit, 
                           year_birth=cd1$year_birth, date_recruit=cd1$date_recruit)
case_summary <- merge(case_summary, death_data[,1:2], by='eid', all.x=T)
##duplicate case is because the case of death can be multiple for one subject.
case_summary=case_summary[!duplicated(case_summary$eid),]
case_summary$Disease <- 0
case_summary$hesin <- 0
case_summary$date_recruit <- as.Date(case_summary$date_recruit)
#disease_case
case_summary <- merge(case_summary, disease_case[,c('eid','epistart')], all.x = T) %>%
  mutate(Disease_epistart=epistart) %>% select(-epistart)
case_summary$hesin[which(!is.na(case_summary$Disease_epistart))] <- 1
case_summary$Disease[which(!is.na(case_summary$Disease_epistart))] <- 1
disease_case_noepi=disease_case[is.na(disease_case$epistart),]
if(nrow(disease_case_noepi)>=1){
  case_summary[case_summary$eid %in% disease_case_noepi$eid,]$Disease <- 3
}
##### ICD cases ID
id_icd <- case_summary$eid[case_summary$Disease>0]

#Disease_res_base
## This is risk 
#case_summary$Disease[Disease_res_base[[1]]] <- 2
#--------------------- Add self-report disease -----------------------#
case_summary=merge(case_summary,s,by="eid")
case_summary$Disease=ifelse(case_summary$self_disease==1,1,case_summary$Disease)
case_summary$date_of_self=case_summary$year_birth+case_summary$self_age
### identify recurrence and incidence for Disease
case_summary$Disease[which((case_summary$Disease_epistart < case_summary$date_recruit)|
                             (case_summary$self_age<=case_summary$age_recruit))] <- 2
#--------------------- Identify self-report/ICD disease -----------------------#
case_summary$Disease_type[case_summary$eid %in% id_icd] <- 'ICD_case'
case_summary$Disease_type[case_summary$eid %in% id_icd==F &
			  case_summary$Disease>0] <- 'Self-report_case'
print(table(case_summary$Disease,case_summary$Disease_type)) 

#https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=Data_providers_and_dates: the censoring date
#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/HospitalEpisodeStatistics.pdf: date of this document 
case_summary$Disease_date_end <-  apply(cbind(as.character(case_summary$date_of_self),
                                              as.character(case_summary$date_of_death),
                                              as.character(case_summary$Disease_epistart),
                                              rep('2023-06-01',length(case_summary$eid)))
                                        ,1,my.min)
case_summary$age_end <- as.numeric(str_sub(case_summary$Disease_date_end,1,4)) - as.numeric(case_summary$year_birth)
case_summary$follow_up <- as.Date(case_summary$Disease_date_end)-as.Date(case_summary$date_recruit) 

names(case_summary) <- str_replace_all(names(case_summary),'Disease',disease)
if (nrow(case_summary) > length(unique(case_summary$eid))) {
  case_summary = case_summary[!duplicated(case_summary$eid),]
}

out_path='/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/disease_refer/'
write.table(case_summary,paste0(out_path,'ukbb_',disease,'.tsv'),
	    quote=F,sep='\t',row.names=F)


## disease clear
library(data.table)

disease="MDD"
input_path = "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/disease_refer/"
output_path = "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/"
case_summary = fread(paste0(input_path,'ukbb_',disease,'.tsv'))

ukbb_MDD = case_summary[,c("eid","MDD")]
ukbb_MDD = na.omit(ukbb_MDD)
ukbb_MDD$MDD = ifelse(ukbb_MDD$MDD == 0, 0, 1)
print(table(ukbb_MDD$MDD))

write.table(ukbb_MDD,paste0(output_path,'ukbb_',disease,'.tsv'),
	    quote=F,sep='\t',row.names=F)