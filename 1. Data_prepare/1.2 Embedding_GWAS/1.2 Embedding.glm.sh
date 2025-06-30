## Step0: Considers hapmap3 SNPs
library(data.table)

snp_1kg = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg/snpinfo_mult_1kg_hm3"))
snplist = data.table(SNP = snp_1kg$SNP)

write.table(snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/PRScsx_snplist.txt"), row.names=F, col.names=F, quote=F, append=F)

## Step1: Generate summary statistics [grace]
## word2vec100
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/word2vec100_glm.txt"
> $job_file  # Empty the job file if it already exists

train_type=train

if [[ ${train_type} == "train" ]]; then
    echo "module load PLINK/2; plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/EUR_hm3 --double-id --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/PRScsx_snplist.txt --keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_train_2over3_doubleid.tsv --pheno /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/individual_embedding/word2vec100_EUR_${train_type}_embed_normalize.txt --covar /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/agesex20PC.csv --covar-variance-standardize --glm hide-covar cols=chrom,pos,alt,ref,a1freq,orbeta,se,tz,p,nobs --out /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm/word2vec100_${train_type}_EUR_UKB" >> $job_file
fi

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/word2vec100_glm.txt --partition=scavenge,day --requeue --mem=100G --cpus-per-task=24 --ntasks=1 --nodes=1 --time=1-00:00:00 --mail-type=ALL
sbatch dsq-word2vec100_glm-$(date +%Y-%m-%d).sh


## gpticd3072
start_col=3
end_col=3074
step=240

job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_glm.txt"
> $job_file  # Empty the job file if it already exists

current_start="${start_col}"
while [[ "${current_start}" -le "${end_col}" ]]; do
current_end=$(( current_start + step - 1 ))
if [[ "${current_end}" -gt "${end_col}" ]]; then
current_end="${end_col}"
fi

train_type=train

if [[ ${train_type} == "train" ]]; then
    echo "module load PLINK/2; plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/EUR_hm3 --double-id --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/PRScsx_snplist.txt --keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_train_2over3_doubleid.tsv --pheno /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/individual_embedding/gpticd3072_EUR_${train_type}_embed_normalize.txt --pheno-col-nums ${current_start}-${current_end} --covar /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/agesex20PC.csv --covar-variance-standardize --glm hide-covar cols=chrom,pos,alt,ref,a1freq,orbeta,se,tz,p,nobs --out /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm/gpticd3072_${train_type}_EUR_UKB" >> $job_file
fi

current_start=$(( current_end + 1 ))

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_glm.txt --partition=week --requeue --mem=100G --cpus-per-task=24 --ntasks=1 --nodes=1 --time=3-00:00:00 --mail-type=ALL
sbatch dsq-gpticd3072_glm-$(date +%Y-%m-%d).sh

## use tmp for running
## choice1
start_col=3
end_col=3074
step=24

job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_glm.txt"

> "$job_file"  # Empty out any old contents

current_start="$start_col"

while [[ "$current_start" -le "$end_col" ]]; do

    current_end=$(( current_start + step - 1 ))

    if [[ "$current_end" -gt "$end_col" ]]; then

        current_end="$end_col"

    fi

    train_type="train"

    # Write a single-line command for dsq
    echo "module load PLINK/2; \
    mkdir -p /tmp/lx94/\${SLURM_JOBID}_${current_start}_${current_end}; \
    cd /tmp/lx94/\${SLURM_JOBID}_${current_start}_${current_end}; \
    cp /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/EUR_hm3.* .; \
    cp /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_train_2over3_doubleid.tsv .; \
    cp /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/PRScsx_snplist.txt .; \
    cp /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/individual_embedding/gpticd3072_EUR_${train_type}_embed_normalize.txt .; \
    cp /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/agesex20PC.csv .; \
    plink2 --bfile EUR_hm3 --double-id --extract PRScsx_snplist.txt \
           --keep EUR_train_2over3_doubleid.tsv \
           --pheno gpticd3072_EUR_${train_type}_embed_normalize.txt \
           --pheno-col-nums ${current_start}-${current_end} \
           --covar agesex20PC.csv --covar-variance-standardize \
           --glm hide-covar cols=chrom,pos,alt,ref,a1freq,orbeta,se,tz,p,nobs \
           --out gpticd3072_${train_type}_EUR_UKB; \
    cp gpticd3072_${train_type}_EUR_UKB.Embedding_* \
       /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm/; \
    cd /; \
    rm -rf /tmp/lx94/\${SLURM_JOBID}_${current_start}_${current_end}" \
    >> "$job_file"
    current_start=$(( current_end + 1 ))

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/
dsq --job-file "$job_file" \
    --partition=scavenge,day \
    --requeue \
    --mem=50G \
    --cpus-per-task=24 \
    --ntasks=1 \
    --nodes=1 \
    --time=1-00:00:00 \
    --mail-type=ALL

sbatch dsq-gpticd3072_glm-$(date +%Y-%m-%d).sh

## choice2
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_glm.txt"
> "$job_file"  # Empty out any old contents

for emb_i in {1..3072}; do
    train_type="train"

    # Write a single-line command for dsq
    if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm/gpticd3072_${train_type}_EUR_UKB.Embedding_${emb_i}.glm.linear" ]]; then
    echo "module load PLINK/2; \
    mkdir -p /tmp/lx94/\${SLURM_JOBID}_Embedding_${emb_i}; \
    cd /tmp/lx94/\${SLURM_JOBID}_Embedding_${emb_i}; \
    cp /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/EUR_hm3.* .; \
    cp /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_train_2over3_doubleid.tsv .; \
    cp /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/PRScsx_snplist.txt .; \
    cp /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/individual_embedding/gpticd3072_EUR_${train_type}_embed_normalize.txt .; \
    cp /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/agesex20PC.csv .; \
    plink2 --bfile EUR_hm3 --double-id --extract PRScsx_snplist.txt \
           --keep EUR_train_2over3_doubleid.tsv \
           --pheno gpticd3072_EUR_${train_type}_embed_normalize.txt \
           --pheno-name Embedding_${emb_i} \
           --covar agesex20PC.csv --covar-variance-standardize \
           --glm hide-covar cols=chrom,pos,alt,ref,a1freq,orbeta,se,tz,p,nobs \
           --out gpticd3072_${train_type}_EUR_UKB; \
    cp gpticd3072_${train_type}_EUR_UKB.Embedding_${emb_i}.* \
       /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm/; \
    cd /; \
    rm -rf /tmp/lx94/\${SLURM_JOBID}_Embedding_${emb_i}" \
    >> "$job_file"
    fi

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/

dsq --job-file "$job_file" \
    --partition=scavenge,day \
    --requeue \
    --mem=50G \
    --cpus-per-task=24 \
    --ntasks=1 \
    --nodes=1 \
    --time=1-00:00:00 \
    --mail-type=ALL

sbatch dsq-gpticd3072_glm-$(date +%Y-%m-%d).sh


## Step2 Clean summary statistics
## word2vec100
library(data.table)

train_type = "train"

for (i in c(1:100)){

    # sumstat obtain
    sumstat = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm/word2vec100_",train_type,"_EUR_UKB.Embedding_",i,".glm.linear"))
    sumstat = sumstat[, A2 := ifelse(REF == A1, ALT, REF)]

    # clean sumstat
    sumstat_clean = sumstat[,c("ID","#CHROM","POS","A1","A2","OBS_CT","A1_FREQ","BETA","SE","T_STAT","P")]
    colnames(sumstat_clean) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","SE","Z","P")
    write.table(sumstat_clean, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/clean/word2vec100_",train_type,"_EUR_UKB_Embedding",i,"_clean.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")
    
    # PRScsx
    PRScsx_clean = sumstat_clean[,c("SNP","A1","A2","BETA","P")]
    write.table(PRScsx_clean, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/PRScsx/word2vec100_",train_type,"_EUR_UKB_Embedding",i,"_PRScsx.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}

## gpticd3072
library(data.table)

train_type = "train"

for (i in c(1:3072)){

    # sumstat obtain
    sumstat = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm/gpticd3072_",train_type,"_EUR_UKB.Embedding_",i,".glm.linear"))
    sumstat = sumstat[, A2 := ifelse(REF == A1, ALT, REF)]

    # clean sumstat
    sumstat_clean = sumstat[,c("ID","#CHROM","POS","A1","A2","OBS_CT","A1_FREQ","BETA","SE","T_STAT","P")]
    colnames(sumstat_clean) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","SE","Z","P")
    write.table(sumstat_clean, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/clean/gpticd3072_",train_type,"_EUR_UKB_Embedding",i,"_clean.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")
    
    # PRScsx
    PRScsx_clean = sumstat_clean[,c("SNP","A1","A2","BETA","P")]
    write.table(PRScsx_clean, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/PRScsx/gpticd3072_",train_type,"_EUR_UKB_Embedding",i,"_PRScsx.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}

## Step3 Clean prune summary statistics
## word2vec100
library(data.table)

pop = "EUR"
train_type = "train"

for (i in c(1:100)){
    # PRScsx for full snplist
    PRScsx_clean = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/PRScsx/word2vec100_",train_type,"_EUR_UKB_Embedding",i,"_PRScsx.txt"))

    # Prune snplist
    prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/",pop,"_prune_pval1_r20.5_wc250_1.snplist"), header = FALSE)

    # PRScsx for prune snplist
    PRScsx_clean_prune_snplist = PRScsx_clean[which(PRScsx_clean$SNP %in% prune_snplist$V1),]

    write.table(PRScsx_clean_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/PRScsx/word2vec100_",train_type,"_EUR_UKB_Embedding",i,"_prune_",pop,"_PRScsx.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}

## gpticd3072
library(data.table)

pop = "EUR"
train_type = "train"

for (i in c(1:3072)){
    # PRScsx for full snplist
    PRScsx_clean = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/PRScsx/gpticd3072_",train_type,"_EUR_UKB_Embedding",i,"_PRScsx.txt"))

    # Prune snplist
    prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/",pop,"_prune_pval1_r20.5_wc250_1.snplist"), header = FALSE)

    # PRScsx for prune snplist
    PRScsx_clean_prune_snplist = PRScsx_clean[which(PRScsx_clean$SNP %in% prune_snplist$V1),]

    write.table(PRScsx_clean_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/PRScsx/gpticd3072_",train_type,"_EUR_UKB_Embedding",i,"_prune_",pop,"_PRScsx.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}

## Step3 Clean unecessary data
## word2vec100
/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm

train_type=train

rm -rf word2vec100_${train_type}_EUR_UKB*.glm.linear
rm -rf word2vec100_${train_type}_EUR_UKB*.log

## gpticd3072
/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/glm

train_type=train

rm -rf gpticd3072_${train_type}_EUR_UKB*.glm.linear
rm -rf gpticd3072_${train_type}_EUR_UKB*.log