## Step1: Generate summary statistics [grace]
## PCA
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/word2vec100_PCA30_glm.txt"
> $job_file  # Empty the job file if it already exists

train_type=train

if [[ ${train_type} == "train" ]]; then
    echo "module load PLINK/2; plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/EUR_hm3 --double-id --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/PRScsx_snplist.txt --keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_train_2over3_doubleid.tsv --pheno /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/individual_embedding/word2vec100_PCA30_EUR_${train_type}_embed_normalize.txt --covar /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/agesex20PC.csv --covar-variance-standardize --glm hide-covar cols=chrom,pos,alt,ref,a1freq,orbeta,se,tz,p,nobs --out /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm/word2vec100_PCA30_${train_type}_EUR_UKB" >> $job_file
fi

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/word2vec100_PCA30_glm.txt --partition=scavenge,day --requeue --mem=100G --cpus-per-task=24 --ntasks=1 --nodes=1 --time=1-00:00:00 --mail-type=ALL
sbatch dsq-word2vec100_PCA30_glm-$(date +%Y-%m-%d).sh

## ICA
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/word2vec100_ICA30_glm.txt"
> $job_file  # Empty the job file if it already exists

train_type=train

if [[ ${train_type} == "train" ]]; then
    echo "module load PLINK/2; plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/EUR_hm3 --double-id --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/PRScsx_snplist.txt --keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_train_2over3_doubleid.tsv --pheno /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/individual_embedding/word2vec100_ICA30_EUR_${train_type}_embed_normalize.txt --covar /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/agesex20PC.csv --covar-variance-standardize --glm hide-covar cols=chrom,pos,alt,ref,a1freq,orbeta,se,tz,p,nobs --out /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm/word2vec100_ICA30_${train_type}_EUR_UKB" >> $job_file
fi


module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/word2vec100_ICA30_glm.txt --partition=scavenge,day --requeue --mem=100G --cpus-per-task=24 --ntasks=1 --nodes=1 --time=1-00:00:00 --mail-type=ALL
sbatch dsq-word2vec100_ICA30_glm-$(date +%Y-%m-%d).sh


## Step2 Clean summary statistics
## PCA
library(data.table)

train_type = "train"

for (i in c(1:30)){

    # sumstat obtain
    sumstat = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm/word2vec100_PCA30_",train_type,"_EUR_UKB.PC",i,".glm.linear"))
    sumstat = sumstat[, A2 := ifelse(REF == A1, ALT, REF)]

    # clean sumstat
    sumstat_clean = sumstat[,c("ID","#CHROM","POS","A1","A2","OBS_CT","A1_FREQ","BETA","SE","T_STAT","P")]
    colnames(sumstat_clean) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","SE","Z","P")
    write.table(sumstat_clean, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/clean/word2vec100_PCA30_",train_type,"_EUR_UKB_PC",i,"_clean.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")
    
    # PRScsx
    PRScsx_clean = sumstat_clean[,c("SNP","A1","A2","BETA","P")]
    write.table(PRScsx_clean, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/PRScsx/word2vec100_PCA30_",train_type,"_EUR_UKB_PC",i,"_PRScsx.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}

## ICA
library(data.table)

train_type = "train"

for (i in c(1:30)){

    # sumstat obtain
    sumstat = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm/word2vec100_ICA30_",train_type,"_EUR_UKB.V",i,".glm.linear"))
    sumstat = sumstat[, A2 := ifelse(REF == A1, ALT, REF)]

    # clean sumstat
    sumstat_clean = sumstat[,c("ID","#CHROM","POS","A1","A2","OBS_CT","A1_FREQ","BETA","SE","T_STAT","P")]
    colnames(sumstat_clean) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","SE","Z","P")
    write.table(sumstat_clean, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/clean/word2vec100_ICA30_",train_type,"_EUR_UKB_V",i,"_clean.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")
    
    # PRScsx
    PRScsx_clean = sumstat_clean[,c("SNP","A1","A2","BETA","P")]
    write.table(PRScsx_clean, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/PRScsx/word2vec100_ICA30_",train_type,"_EUR_UKB_V",i,"_PRScsx.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}


## Step3 Clean prune summary statistics
## PCA
library(data.table)

pop = "EUR"
train_type = "train"

for (i in c(1:30)){
    # PRScsx for full snplist
    PRScsx_clean = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/PRScsx/word2vec100_PCA30_",train_type,"_EUR_UKB_PC",i,"_PRScsx.txt"))

    # Prune snplist
    prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/",pop,"_prune_pval1_r20.5_wc250_1.snplist"), header = FALSE)

    # PRScsx for prune snplist
    PRScsx_clean_prune_snplist = PRScsx_clean[which(PRScsx_clean$SNP %in% prune_snplist$V1),]

    write.table(PRScsx_clean_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/PRScsx/word2vec100_PCA30_",train_type,"_EUR_UKB_PC",i,"_prune_",pop,"_PRScsx.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}

## ICA
library(data.table)

for (pop in c("EUR")){

if (pop == "EUR"){
    train_type = "train"
}

for (i in c(1:30)){
    # PRScsx for full snplist
    PRScsx_clean = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/PRScsx/word2vec100_ICA30_",train_type,"_EUR_UKB_V",i,"_PRScsx.txt"))

    # Prune snplist
    prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/",pop,"_prune_pval1_r20.5_wc250_1.snplist"), header = FALSE)

    # PRScsx for prune snplist
    PRScsx_clean_prune_snplist = PRScsx_clean[which(PRScsx_clean$SNP %in% prune_snplist$V1),]

    write.table(PRScsx_clean_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/PRScsx/word2vec100_ICA30_",train_type,"_EUR_UKB_V",i,"_prune_",pop,"_PRScsx.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}
}

## Step3 Clean unecessary data
## PCA
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm

for train_type in full train; do
rm -rf word2vec100_PCA30_${train_type}_EUR_UKB*.glm.linear
rm -rf word2vec100_PCA30_${train_type}_EUR_UKB*.log
done

## ICA
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm

for train_type in full train; do
rm -rf word2vec100_ICA30_${train_type}_EUR_UKB*.glm.linear
rm -rf word2vec100_ICA30_${train_type}_EUR_UKB*.log
done