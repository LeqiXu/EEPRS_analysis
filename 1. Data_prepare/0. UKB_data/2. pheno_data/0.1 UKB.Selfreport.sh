## self report for mental disease: 20544
library(data.table)
library(stringr)

table <- fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/2001150/ukb671492.tab",header=T,quote='\t')
namelist <- colnames(table)

field <- c()
for ( i in c("20544")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)

self_report_20544 <- table[,..field]
write.table(self_report_20544, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/self_report_20544.tsv", row.names=F, col.names=T, quote = F)


## self report for cancer: 20001
library(data.table)
library(stringr)

table <- fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/10570/ukb671493.tab",header=T,quote='\t')
namelist <- colnames(table)

field <- c()
for ( i in c("20001")){
add <- namelist[!is.na(str_extract(namelist,i))]
field <- c(field,add)
}
field <- c("f.eid",unique(field))
print(field)

self_report_20001 <- table[,..field]
write.table(self_report_20001, "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/self_report_20001.tsv", row.names=F, col.names=T, quote = F)