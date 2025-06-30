## Step1: Obtain clean GWAS for MIX
library(data.table)

pop = "EUR"

trait_list = c("HDL","LDL","TC","logTG","BMI","WHR","AgeSmk","SmkI","SmkC","CigDay","DrnkWk","Glu2h","HbA1c","eGFR","logCp",
                "AD","Angina","Asthma","AF","ADHD","ASD","BIP","BrC","CKD","CAD","HF","IBD","IS","LuC",
                "MDD","Osteoporosis","OvC","PAD","PBC","PSC","PrC","PAH","RA","SCZ","SLE","T2D")

for (trait in trait_list){

sumstat_clean <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/clean/",trait,"_",pop,"_all_clean.txt"))

sumstat_MIX = sumstat_clean[,c("SNP","A1","A2","BETA","SE","Z","P","N")]

write.table(sumstat_MIX,paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/MIX/",trait, "_", pop, "_all_MIX.txt"),quote=F,sep='\t',row.names=F,col.names=T)
}


## Step2: Subsample GWAS by MIX
pop=EUR
i=1

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}
do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/subsample_data/clean/${trait}_prune_${pop}_tune_GWAS_approxTRUE_ratio3.00_repeat4.txt" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_subsample_generate_${trait}_prune
#SBATCH --output=out_GWAS_subsample_generate_${trait}_prune.txt

module load miniconda
conda activate py_env

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--sst_file=/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/MIX/${trait}_${pop}_all_MIX.txt \
--pop=${pop} \
--prune_snplist=/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/${pop}_prune_pval1_r20.5_wc250_${i}.snplist \
--indep_approx=TRUE \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/subsample_data/clean \
--out_name=${trait}_prune

EOT
fi
done


## Step3: Format training GWAS and tuning GWAS
## training GWAS
library(readr)
library(data.table)

pop="EUR"

trait_list = c("HDL","LDL","TC","logTG","BMI","WHR","AgeSmk","SmkI","SmkC","CigDay","DrnkWk","Glu2h","HbA1c","eGFR","logCp",
                "AD","Angina","Asthma","AF","ADHD","ASD","BIP","BrC","CKD","CAD","HF","IBD","IS","LuC",
                "MDD","Osteoporosis","OvC","PAD","PBC","PSC","PrC","PAH","RA","SCZ","SLE","T2D")

for (trait in trait_list){
for (rpt in c(1:4)){
# clean sumstat
sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/subsample_data/clean/",trait,"_prune_",pop,"_train_GWAS_approxTRUE_ratio3.00_repeat",rpt,".txt"))

PRScsx_data = sumstat_data[,c("SNP","A1","A2","BETA","P")]
print(paste0(trait, " SNP number: ", nrow(PRScsx_data)))
write.table(PRScsx_data, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/subsample_data/PRScsx/",trait,"_prune_",pop,"_train_PRScsx_repeat",rpt,".txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")
}
}


## tuning GWAS
library(data.table)

pop="EUR"

trait_list = c("HDL","LDL","TC","logTG","BMI","WHR","AgeSmk","SmkI","SmkC","CigDay","DrnkWk","Glu2h","HbA1c","eGFR","logCp",
                "AD","Angina","Asthma","AF","ADHD","ASD","BIP","BrC","CKD","CAD","HF","IBD","IS","LuC",
                "MDD","Osteoporosis","OvC","PAD","PBC","PSC","PrC","PAH","RA","SCZ","SLE","T2D")

for (trait in trait_list){
for (rpt in c(1:4)){

sumstat_data_prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/subsample_data/clean/",trait,"_prune_",pop,"_tune_GWAS_approxTRUE_ratio3.00_repeat",rpt,".txt"))

sumstat_MIX_prune_snplist = sumstat_data_prune_snplist[,c("SNP","A1","A2","BETA","SE","Z","P","N")]
write.table(sumstat_MIX_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/subsample_data/MIX/",trait,"_prune_",pop,"_tune_MIX_approxTRUE_ratio3.00_repeat",rpt,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}
}