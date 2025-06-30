# 1 Obtain basic phenotype from ukbb
## obtain each phenotype
library(data.table)
library(stringr)

## remove id
remove_id = fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/w29900_2023-04-25.csv")

## phenotype
## DrnkWk: Drinks per week: 1558
table <- fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/10570/ukb671493.tab",header=T,quote='\t')
namelist <- colnames(table)

field <- c()
for ( i in c("f.1558.0.0")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)

DrnkWk_data <- table[,..field]
colnames(DrnkWk_data) = c("eid","DrnkWk")
DrnkWk_data = DrnkWk_data[which(!(DrnkWk_data$eid %in% remove_id $V1)),]
DrnkWk_data$DrnkWk = as.numeric(DrnkWk_data$DrnkWk)

ukbb_DrnkWk = DrnkWk_data[,c("eid","DrnkWk")]
ukbb_DrnkWk = na.omit(ukbb_DrnkWk)
write.table(ukbb_DrnkWk, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_DrnkWk.tsv", row.names=F, col.names=T, quote = F)
