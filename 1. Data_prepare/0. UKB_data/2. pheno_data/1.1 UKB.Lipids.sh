# 1 Obtain basic phenotype from ukbb
## obtain each phenotype
library(data.table)
library(stringr)

## remove id
remove_id = fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/w29900_2023-04-25.csv")

## phenotype
## lipids: 
## HDL: 30760
## LDL: 30780
## TC: 30690
## TG: 30870
table <- fread('/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_pheno/update/ukb42584.tab',header=T,quote='\t')

namelist <- colnames(table)
field <- c()
for ( i in c("f.30760.0.0","f.30780.0.0","f.30690.0.0","f.30870.0.0")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)

lipids_data <- table[,..field]
colnames(lipids_data) = c("eid","HDL","LDL","TC","TG") 
lipids_data = lipids_data[which(!(lipids_data$eid %in% remove_id $V1)),]
lipids_data$HDL = as.numeric(lipids_data$HDL)
lipids_data$LDL = as.numeric(lipids_data$LDL)
lipids_data$TC = as.numeric(lipids_data$TC)
lipids_data$TG = as.numeric(lipids_data$TG)

field = fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/10570/fields.ukb")
field_data = fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/10570/ukb671493.tab",header=T,quote='\t')

med_male_code = "6177"
med_male_code_s = "f.6177.0.0"
med_female_code = "6153"
med_female_code_s = "f.6153.0.0"
field_name = c("f.eid",med_male_code_s,med_female_code_s)
med_data = field_data[,..field_name]
med_data$med = ifelse(is.na(med_data$f.6177.0.0),med_data$f.6153.0.0,med_data$f.6177.0.0)
med_data = med_data[,c("f.eid","med")]
colnames(med_data) = c("eid","med")
med_data = na.omit(med_data)

clm_eid = med_data$eid[which(med_data$med %in% c("1"))]
no_clm_eid = med_data$eid[which(med_data$med %in% c("2","3","4","5","-7"))]

## HDL
ukbb_HDL = lipids_data[,c("eid","HDL")]
ukbb_HDL = na.omit(ukbb_HDL)
write.table(ukbb_HDL, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_HDL.tsv", row.names=F, col.names=T, quote = F)

## LDL
ukbb_LDL = lipids_data[,c("eid","LDL")]
ukbb_LDL = na.omit(ukbb_LDL)
ukbb_LDL = ukbb_LDL[which(ukbb_LDL$eid %in% c(clm_eid,no_clm_eid)),]
ukbb_LDL$LDL[which(ukbb_LDL$eid %in% clm_eid)] = ukbb_LDL$LDL[which(ukbb_LDL$eid %in% clm_eid)] / 0.7
write.table(ukbb_LDL, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_LDL.tsv", row.names=F, col.names=T, quote = F)

## TC
ukbb_TC = lipids_data[,c("eid","TC")]
ukbb_TC = na.omit(ukbb_TC)
ukbb_TC = ukbb_TC[which(ukbb_TC$eid %in% c(clm_eid,no_clm_eid)),]
ukbb_TC$TC[which(ukbb_TC$eid %in% clm_eid)] = ukbb_TC$TC[which(ukbb_TC$eid %in% clm_eid)] / 0.8
write.table(ukbb_TC, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_TC.tsv", row.names=F, col.names=T, quote = F)

## logTG
ukbb_TG = lipids_data[,c("eid","TG")]
ukbb_TG = na.omit(ukbb_TG)
ukbb_logTG = ukbb_TG
ukbb_logTG$logTG = log(ukbb_TG$TG)
ukbb_logTG = ukbb_logTG[,c("eid","logTG")]
write.table(ukbb_logTG, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_logTG.tsv", row.names=F, col.names=T, quote = F)