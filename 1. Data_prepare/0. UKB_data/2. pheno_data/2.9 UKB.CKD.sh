## Step1: disease info
library(stringr)
library(data.table)
library(dplyr)
library(magrittr)
source('/gpfs/gibbs/pi/zhao/cz354/UKB/pheno/ukbcase_v2.R')
source('/gpfs/gibbs/pi/zhao/cz354/UKB/pheno/function.R')

disease='CKD'
CKD_icd10 <- c("D593","D638","E102","E112","E132","E142",
               paste0("E85",c(0:4,8,9)),
               paste0("I12",c(8,9)),
               paste0("I13",c(1,2)),
               "I151","M103","N001","N011","N020","N021",
               "N030","N031",paste0("N04",c(0:9)),
               "N050","N051","N060","N061","N071",
               paste0("N08",c(0:5)),
               "N10",paste0("N11",c(0,1,8,9)),"N12",
               paste0("N14",c(1:4)),"N165",
               paste0("N17",c(0:2,8,9)),
               paste0("N18",c(0:5,9)),"N19",
               "N250","N990",paste0("Q61",c(0:5,8,9)),
               "T824","T861","Y602","Y841",paste0("Z49",c(0:2)),
               "Z940","Z992")
CKD_icd9 <- NULL
CKD_ops <- NULL
CKD_s1 <- NULL
CKD_s2 <- NULL


hesin=fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/ukb_hesin.tsv")
hesin_diag9=fread('/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/ukb_hesin_diag9.tsv')
hesin_diag10=fread('/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/ukb_hesin_diag10.tsv')
hesin_oper=fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/hesin_oper.txt")

death_data=fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/death_data.txt")

selfreport <- fread('/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/self_report.csv')
disease_age <- fread('/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/disease_age.csv',header=T)
cd1=fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/cd1.txt")

# from self-report, here disease_age1/2 is age corresponding to the disease in th selfreport1/2
field_selfreport1 <- c("f.eid",colnames(selfreport)[!is.na(str_extract(colnames(selfreport),'20002'))])
field_selfreport2 <- c('f.eid',colnames(selfreport)[!is.na(str_extract(colnames(selfreport),'20004'))])
field_diseaseage1 <- c('f.eid',colnames(disease_age)[!is.na(str_extract(colnames(disease_age),'20009'))])
field_diseaseage2 <- c('f.eid',colnames(disease_age)[!is.na(str_extract(colnames(disease_age),'20011'))])
selfreport1 <- selfreport[,..field_selfreport1]
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
#Disease_res_base
## This is risk 
#case_summary$Disease[Disease_res_base[[1]]] <- 2
case_summary=merge(case_summary,s,by="eid")
case_summary$Disease=ifelse(case_summary$self_disease==1,1,case_summary$Disease)
case_summary$date_of_self=ifelse(case_summary$self_disease==0,
				 case_summary$age_recruit,
				 case_summary$year_birth+case_summary$self_age)
### identify recurrence and incidence for Disease
case_summary$Disease[which((case_summary$Disease_epistart < case_summary$date_recruit)|
                             (case_summary$self_age<=case_summary$age_recruit)&
			     case_summary$Disease==1)] <- 2
print(table(case_summary$Disease)) 

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

disease="CKD"
input_path = "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/disease_refer/"
output_path = "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/"
case_summary = fread(paste0(input_path,'ukbb_',disease,'.tsv'))

ukbb_CKD = case_summary[,c("eid","CKD")]
ukbb_CKD = na.omit(ukbb_CKD)
ukbb_CKD$CKD = ifelse(ukbb_CKD$CKD == 0, 0, 1)
print(table(ukbb_CKD$CKD))

write.table(ukbb_CKD,paste0(output_path,'ukbb_',disease,'.tsv'),
	    quote=F,sep='\t',row.names=F)