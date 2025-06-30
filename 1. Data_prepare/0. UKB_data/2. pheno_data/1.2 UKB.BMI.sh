# 1 Obtain basic phenotype from ukbb
## obtain each phenotype
library(data.table)
library(stringr)

## remove id
remove_id = fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/w29900_2023-04-25.csv")

## phenotype
## BMI: 21001
table <- fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/10570/ukb671493.tab",header=T,quote='\t')
namelist <- colnames(table)

field <- c()
for ( i in c("f.21001.0.0")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)

BMI_data <- table[,..field]
colnames(BMI_data) = c("eid","BMI")
BMI_data = BMI_data[which(!(BMI_data$eid %in% remove_id $V1)),]
BMI_data$BMI = as.numeric(BMI_data$BMI)

ukbb_BMI = BMI_data[,c("eid","BMI")]
ukbb_BMI = na.omit(ukbb_BMI)
write.table(ukbb_BMI, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_BMI.tsv", row.names=F, col.names=T, quote = F)

## WHR: Waist circumference / Hip circumference
## Waist circumference: 48
## Hip circumference: 49
field <- c()
for ( i in c("f.48.0.0","f.49.0.0")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)
field = field[-4]
print(field)

WHR_data <- table[,..field]
colnames(WHR_data) = c("eid","Waist","Hip")
WHR_data = WHR_data[which(!(WHR_data$eid %in% remove_id $V1)),]
WHR_data$Waist = as.numeric(WHR_data$Waist)
WHR_data$Hip = as.numeric(WHR_data$Hip)
WHR_data$WHR = WHR_data$Waist / WHR_data$Hip

ukbb_WHR = WHR_data[,c("eid","WHR")]
ukbb_WHR = na.omit(ukbb_WHR)
write.table(ukbb_WHR, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_WHR.tsv", row.names=F, col.names=T, quote = F)
