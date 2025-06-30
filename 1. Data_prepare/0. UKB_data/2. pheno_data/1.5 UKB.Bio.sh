# 1 Obtain basic phenotype from ukbb
## obtain each phenotype
library(data.table)
library(stringr)

## remove id
remove_id = fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/w29900_2023-04-25.csv")

## phenotype
## Glu2h: 2-hour glucose: 30740
table <- fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/2008682/ukb42584.tab",header=T,quote='\t')
namelist <- colnames(table)

field <- c()
for ( i in c("f.30740.0.0")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)

Glu2h_data <- table[,..field]
colnames(Glu2h_data) = c("eid","Glu2h")
Glu2h_data = Glu2h_data[which(!(Glu2h_data$eid %in% remove_id $V1)),]
Glu2h_data$Glu2h = as.numeric(Glu2h_data$Glu2h)

ukbb_Glu2h = Glu2h_data[,c("eid","Glu2h")]
ukbb_Glu2h = na.omit(ukbb_Glu2h)
write.table(ukbb_Glu2h, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_Glu2h.tsv", row.names=F, col.names=T, quote = F)


## HbA1c: Glycated haemoglobin: 30750
field <- c()
for ( i in c("f.30750.0.0")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)

HbA1c_data <- table[,..field]
colnames(HbA1c_data) = c("eid","HbA1c")
HbA1c_data = HbA1c_data[which(!(HbA1c_data$eid %in% remove_id $V1)),]
HbA1c_data$HbA1c = as.numeric(HbA1c_data$HbA1c)

ukbb_HbA1c = HbA1c_data[,c("eid","HbA1c")]
ukbb_HbA1c = na.omit(ukbb_HbA1c)
write.table(ukbb_HbA1c, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_HbA1c.tsv", row.names=F, col.names=T, quote = F)


## eGFR: Glycated haemoglobin: 30700
field <- c()
for ( i in c("f.30700.0.0")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)

eGFR_data <- table[,..field]
colnames(eGFR_data) = c("eid","eGFR")
eGFR_data = eGFR_data[which(!(eGFR_data$eid %in% remove_id $V1)),]
eGFR_data$eGFR = as.numeric(eGFR_data$eGFR)

ukbb_eGFR = eGFR_data[,c("eid","eGFR")]
ukbb_eGFR = na.omit(ukbb_eGFR)
write.table(ukbb_eGFR, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_eGFR.tsv", row.names=F, col.names=T, quote = F)


## logCp: Glycated haemoglobin: 30710
field <- c()
for ( i in c("f.30710.0.0")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)

logCp_data <- table[,..field]
colnames(logCp_data) = c("eid","Cp")
logCp_data = logCp_data[which(!(logCp_data$eid %in% remove_id $V1)),]
logCp_data$Cp = as.numeric(logCp_data$Cp)
logCp_data$logCp = log(logCp_data$Cp)

ukbb_logCp = logCp_data[,c("eid","logCp")]
ukbb_logCp = na.omit(ukbb_logCp)
write.table(ukbb_logCp, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_logCp.tsv", row.names=F, col.names=T, quote = F)
