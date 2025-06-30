# 1. Estimate beta from PRScsx
for chr in {1..22}; do
for pop2 in EUR; do

job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/gpt_prune_${pop2}_PRScsx_beta_chr${chr}.txt"
> $job_file  # Empty the job file if it already exists

pop=EUR
param_phi=auto

for i in {1..3072}; do
trait="Embedding${i}"

if [[ ${pop2} == "EUR" ]]; then
train_type="train"; sample_size=207734
fi

output_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/gpticd3072_${train_type}_${pop}_UKB_${trait}_prune_${pop2}_PRScsx_${pop}_pst_eff_a1_b0.5_phi${param_phi}_chr${chr}.txt"
if [[ ! -e ${output_file} ]]; then
        echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/; python /gpfs/gibbs/pi/zhao/lx94/PRSmap/method/PRScsx/PRScsx.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/embedding_data/PRScsx/gpticd3072_${train_type}_${pop}_UKB_${trait}_prune_${pop2}_PRScsx.txt --n_gwas=${sample_size} --chrom=${chr} --pop=${pop} --out_dir=result/Embedding_PRS/ --out_name=gpticd3072_${train_type}_${pop}_UKB_${trait}_prune_${pop2}_PRScsx" >> $job_file
fi

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/
dsq --job-file gpt_prune_${pop2}_PRScsx_beta_chr${chr}.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-gpt_prune_${pop2}_PRScsx_beta_chr${chr}-$(date +%Y-%m-%d).sh

done
done

# 2. Obrain beta by chr pop for each param in each trait
library(data.table)

pop="EUR"

for (i in c(1:3072)){
trait = paste0("Embedding",i)

for (pop2 in c("EUR")){

if (pop2 == "EUR"){
train_type="train"
}

PRScsx_all <- data.table()
for (chr in 1:22){

PRScsx_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/gpticd3072_",train_type,"_",pop,"_UKB_",trait,"_prune_",pop2,"_PRScsx_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))

PRScsx_pop_chr <- PRScsx_pop_chr[,c(2,4,6)]
names(PRScsx_pop_chr) = c("SNP","A1",pop)

PRScsx_all = rbind(PRScsx_all,PRScsx_pop_chr)
    
}

write.table(PRScsx_all,paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/gpticd3072_",train_type,"_",pop,"_UKB_",trait,"_prune_",pop2,"_PRScsx.txt"),quote=F,sep='\t',row.names=F,col.names=T)

}

}


# 3. Clean previous result
for i in {1..3072}; do
trait="Embedding${i}"

for pop2 in EUR; do


if [[ ${pop2} == "EUR" ]]; then
train_type="train"
fi

for chr in {1..22}; do

rm -rf /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/gpticd3072_${train_type}_${pop}_UKB_${trait}_prune_${pop2}_PRScsx_${pop}_pst_eff_a1_b0.5_phi${param_phi}_chr${chr}.txt

done
done
done