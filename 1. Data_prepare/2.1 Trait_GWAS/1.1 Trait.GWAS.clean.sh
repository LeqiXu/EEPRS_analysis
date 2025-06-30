# 1.  Sumstat with no UKB [phenotype existed in UKB]
## Organize data
pop = "EUR"

trait_list = c("HDL","LDL","TC","logTG","BMI","WHR","AgeSmk","SmkI","SmkC","CigDay","DrnkWk","Glu2h","HbA1c","eGFR","logCp",
                "AD","Angina","Asthma","AF","ADHD","ASD","BIP","BrC","CKD","CAD","HF","IBD","IS","LuC",
                "MDD","Osteoporosis","OvC","PAD","PBC","PSC","PrC","PAH","RA","SCZ","SLE","T2D","HOMA_BadjBMI","BUN","FI","FPadjBMI","FG","fibrinogen","")

for (trait in trait_list){

if (trait %in% c("HDL","LDL","TC","logTG","IS")) {
  jiaqi_trait = tolower(trait)
} else if (trait == "WHR") {
  jiaqi_trait = "WHRadjBMI"
} else if (trait == "SmkI") {
  jiaqi_trait = "smk2"
} else if (trait == "SmkC") {
  jiaqi_trait = "SmkCes"
} else if (trait == "Glu2h") {
  jiaqi_trait = "2hGlu"
} else if (trait == "logCp") {
  jiaqi_trait = "CRP"
} else if (trait == "Angina") {
  jiaqi_trait = "I9_ANGINA"
} else if (trait == "AF") {
  jiaqi_trait = "af2"
} else if (trait == "BrC") {
  jiaqi_trait = "breast_cancer"
} else if (trait == "CAD") {
  jiaqi_trait = "cad_cardiogramc4a"
} else if (trait == "HF") {
  jiaqi_trait = "I9_HEARTFAIL"
} else if (trait == "LuC") {
  jiaqi_trait = "lung_cancer"
} else if (trait == "Osteoporosis") {
  jiaqi_trait = "M13_OSTEOPOROSIS"
} else if (trait == "OvC") {
  jiaqi_trait = "ovarian_cancer"
} else if (trait == "PAD") {
  jiaqi_trait = "pad_mvp"
} else if (trait == "PrC") {
  jiaqi_trait = "prostate_cancer"
} else if (trait == "T2D") {
  jiaqi_trait = "t2d_mahajan2018b"
} else {
  jiaqi_trait = trait
}

sumstat <- readRDS(paste0("/gpfs/gibbs/pi/zhao/jh2875/shared/prs/sumstats/cleaned_sumstats/",jiaqi_trait,"_sum_rm_dup_overlapUKB.RDS"))
sumstat$Z = sumstat$BETA / sumstat$SE

if (trait %in% c("AF","PAH")) {
sumstat$SNP = sumstat$V2
} else if (trait == "CAD") {
sumstat$N = 184305
} else if (trait == "PAD") {
sumstat$N = 174993
} else if (trait == "T2D") {
sumstat$SNP = sumstat$V2
sumstat$N = 455313
}

## clean
if ("EAF" %in% colnames(sumstat)) {
sumstat_clean = sumstat[,c("SNP","CHR","POS","EA","NEA","N","EAF","BETA","SE","Z","P")]
colnames(sumstat_clean) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","SE","Z","P")
write.table(sumstat_clean, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/clean/",trait,"_",pop,"_all_clean.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

} else {
sumstat_clean = sumstat[,c("SNP","CHR","POS","EA","NEA","N","BETA","SE","Z","P")]
colnames(sumstat_clean) = c("SNP","CHR","POS","A1","A2","N","BETA","SE","Z","P")
write.table(sumstat_clean, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/clean/",trait,"_",pop,"_all_clean.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}

## ldsc
sumstat_ldsc = sumstat[,c("SNP","EA","NEA","N","P","Z")]
colnames(sumstat_ldsc) = c("SNP","A1","A2","N","P","Z")

write.table(sumstat_ldsc, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/ldsc/",trait,"_",pop,"_all_ldsc.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

## PRScsx
sumstat_PRScsx = sumstat[,c("SNP","EA","NEA","BETA","P")]
colnames(sumstat_PRScsx) = c("SNP","A1","A2","BETA","P")

write.table(sumstat_PRScsx, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/PRScsx/",trait,"_",pop,"_all_PRScsx.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")
}


## Clean sumstat with incorrect column
pop=EUR

for trait in Osteoporosis HF Angina; do

  # Filter lines with exactly 11 fields, and overwrite
  awk 'NF == 11' /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/clean/${trait}_${pop}_all_clean.txt > /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/clean/${trait}_${pop}_all_clean.tmp
  mv /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/clean/${trait}_${pop}_all_clean.tmp /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/clean/${trait}_${pop}_all_clean.txt

  # Filter lines with exactly 6 fields, and overwrite
  awk 'NF == 6' /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/ldsc/${trait}_${pop}_all_ldsc.txt > /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/ldsc/${trait}_${pop}_all_ldsc.tmp
  mv /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/ldsc/${trait}_${pop}_all_ldsc.tmp /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/ldsc/${trait}_${pop}_all_ldsc.txt

  # Filter lines with exactly 5 fields, and overwrite
  awk 'NF == 5' /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/PRScsx/${trait}_${pop}_all_PRScsx.txt > /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/PRScsx/${trait}_${pop}_all_PRScsx.tmp
  mv /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/PRScsx/${trait}_${pop}_all_PRScsx.tmp /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/PRScsx/${trait}_${pop}_all_PRScsx.txt

done


# 2. Sumstat with no UKB [phenotype not existed in UKB]
## Organize data
pop = "EUR"

trait_list = c("homocysteine","fibrinogen","BUN","IFCadjBMI","ISIadjBMI","FPadjBMI","FG","FI","HOMA_BadjBMI","HOMA_IRadjBMI")

for (trait in trait_list){

sumstat <- readRDS(paste0("/gpfs/gibbs/pi/zhao/jh2875/shared/prs/sumstats/cleaned_sumstats/",trait,"_sum_rm_dup_overlapUKB.RDS"))
sumstat$Z = sumstat$BETA / sumstat$SE

if (trait %in% c("homocysteine","IFCadjBMI","ISIadjBMI","FPadjBMI")) {
sumstat$SNP = sumstat$V2
} 

## clean
if ("EAF" %in% colnames(sumstat)) {
sumstat_clean = sumstat[,c("SNP","CHR","POS","EA","NEA","N","EAF","BETA","SE","Z","P")]
colnames(sumstat_clean) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","SE","Z","P")
write.table(sumstat_clean, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/clean/",trait,"_",pop,"_all_clean.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

} else {
sumstat_clean = sumstat[,c("SNP","CHR","POS","EA","NEA","N","BETA","SE","Z","P")]
colnames(sumstat_clean) = c("SNP","CHR","POS","A1","A2","N","BETA","SE","Z","P")
write.table(sumstat_clean, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/clean/",trait,"_",pop,"_all_clean.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}

## ldsc
sumstat_ldsc = sumstat[,c("SNP","EA","NEA","N","P","Z")]
colnames(sumstat_ldsc) = c("SNP","A1","A2","N","P","Z")

write.table(sumstat_ldsc, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/ldsc/",trait,"_",pop,"_all_ldsc.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

## PRScsx
sumstat_PRScsx = sumstat[,c("SNP","EA","NEA","BETA","P")]
colnames(sumstat_PRScsx) = c("SNP","A1","A2","BETA","P")

write.table(sumstat_PRScsx, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/PRScsx/",trait,"_",pop,"_all_PRScsx.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")
}

# 3. Clean data
module load miniconda
conda activate ldsc

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D homocysteine fibrinogen BUN IFCadjBMI ISIadjBMI FPadjBMI FG FI HOMA_BadjBMI HOMA_IRadjBMI"

for trait in ${trait_list}; do
python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/munge_sumstats.py \
--sumstats /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/ldsc/${trait}_EUR_all_ldsc.txt \
--out /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/munge_data/${trait}_EUR_ldsc
done