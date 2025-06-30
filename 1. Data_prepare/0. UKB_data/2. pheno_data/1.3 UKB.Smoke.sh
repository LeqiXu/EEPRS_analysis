# 1 Obtain basic phenotype from ukbb
## obtain each phenotype
library(data.table)
library(stringr)

## remove id
remove_id = fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/w29900_2023-04-25.csv")

## phenotype
## AgeSmk: Age started smoking in current smokers or Age started smoking in former smokers
## Age started smoking in current smokers: 3436
## Age started smoking in former smokers: 2867
table <- fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/10570/ukb671493.tab",header=T,quote='\t')
namelist <- colnames(table)

field <- c()
for ( i in c("f.3436.0.0","f.2867.0.0")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)

AgeSmk_data <- table[,..field]
colnames(AgeSmk_data) = c("eid","AgeSmkC","AgeSmkF")
AgeSmk_data = AgeSmk_data[which(!(AgeSmk_data$eid %in% remove_id $V1)),]
AgeSmk_data$AgeSmkC = as.numeric(AgeSmk_data$AgeSmkC)
AgeSmk_data$AgeSmkF = as.numeric(AgeSmk_data$AgeSmkF)
AgeSmk_data$AgeSmk = AgeSmk_data$AgeSmkC
AgeSmk_data$AgeSmk[is.na(AgeSmk_data$AgeSmk)] = AgeSmk_data$AgeSmkF[is.na(AgeSmk_data$AgeSmk)]

ukbb_AgeSmk = AgeSmk_data[,c("eid","AgeSmk")]
ukbb_AgeSmk = na.omit(ukbb_AgeSmk)
write.table(ukbb_AgeSmk, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_AgeSmk.tsv", row.names=F, col.names=T, quote = F)


## Smk: 20116
## Smoking Initiation (SmkI: current+ever vs. never): 1+2 vs. 0
## Smoking Cessation (SmkC: current vs. never): 2 vs. 0
field <- c()
for ( i in c("f.20116.0.0")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)

Smk_data <- table[,..field]
colnames(Smk_data) = c("eid","Smk")
Smk_data = Smk_data[which(!(Smk_data$eid %in% remove_id $V1)),]
Smk_data$Smk = as.integer(Smk_data$Smk)
Smk_data = na.omit(Smk_data)
Smk_data = Smk_data[which(Smk_data$Smk != -3),]

ukbb_SmkI = Smk_data[,c("eid","Smk")]
ukbb_SmkI$SmkI = ifelse(ukbb_SmkI$Smk == 0, 0, 1)
ukbb_SmkI = ukbb_SmkI[,c("eid","SmkI")]
write.table(ukbb_SmkI, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_SmkI.tsv", row.names=F, col.names=T, quote = F)

ukbb_SmkC = Smk_data[,c("eid","Smk")]
ukbb_SmkC = ukbb_SmkC[which(ukbb_SmkC$Smk != 1),]
ukbb_SmkC$SmkC = ifelse(ukbb_SmkC$Smk == 2, 1, 0)
ukbb_SmkC = ukbb_SmkC[,c("eid","SmkC")]
write.table(ukbb_SmkC, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_SmkC.tsv", row.names=F, col.names=T, quote = F)


## CigDay: Cigarettes per day
## Number of cigarettes currently smoked daily (current cigarette smokers): 3456
## Number of cigarettes perviously smoked daily: 2887
field <- c()
for ( i in c("f.3456.0.0","f.2887.0.0")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)

CigDay_data <- table[,..field]
colnames(CigDay_data) = c("eid","CigDayC","CigDayF")
CigDay_data = CigDay_data[which(!(CigDay_data$eid %in% remove_id $V1)),]
CigDay_data$CigDayC = as.numeric(CigDay_data$CigDayC)
CigDay_data$CigDayF = as.numeric(CigDay_data$CigDayF)
CigDay_data$CigDay = CigDay_data$CigDayC
CigDay_data$CigDay[is.na(CigDay_data$CigDay)] = CigDay_data$CigDayF[is.na(CigDay_data$CigDay)]

ukbb_CigDay = CigDay_data[,c("eid","CigDay")]
ukbb_CigDay = na.omit(ukbb_CigDay)
write.table(ukbb_CigDay, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_CigDay.tsv", row.names=F, col.names=T, quote = F)
