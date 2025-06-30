# Step0: Create folder
trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do
mkdir /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MTAG_meta_PRS/${trait}
done

# Step1: hm3 PRScsx beta
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/${trait}_multi_trait_beta.txt"
> $job_file  # Empty the job file if it already exists

# sample size
if [[ ${trait} == "HDL" ]]; then
sample_size=950886
elif [[ ${trait} == "LDL" ]]; then
sample_size=950886
elif [[ ${trait} == "TC" ]]; then
sample_size=950886
elif [[ ${trait} == "logTG" ]]; then
sample_size=950886
elif [[ ${trait} == "BMI" ]]; then
sample_size=339224
elif [[ ${trait} == "WHR" ]]; then
sample_size=224459
elif [[ ${trait} == "AgeSmk" ]]; then
sample_size=175835
elif [[ ${trait} == "SmkI" ]]; then
sample_size=357235
elif [[ ${trait} == "SmkC" ]]; then
sample_size=188701
elif [[ ${trait} == "CigDay" ]]; then
sample_size=183196
elif [[ ${trait} == "DrnkWk" ]]; then
sample_size=304322
elif [[ ${trait} == "Glu2h" ]]; then
sample_size=63396
elif [[ ${trait} == "HbA1c" ]]; then
sample_size=146806
elif [[ ${trait} == "eGFR" ]]; then
sample_size=567460
elif [[ ${trait} == "logCp" ]]; then
sample_size=204402
elif [[ ${trait} == "AD" ]]; then
sample_size=63926
elif [[ ${trait} == "Angina" ]]; then
sample_size=418385
elif [[ ${trait} == "Asthma" ]]; then
sample_size=127669
elif [[ ${trait} == "AF" ]]; then
sample_size=635097
elif [[ ${trait} == "ADHD" ]]; then
sample_size=225534
elif [[ ${trait} == "ASD" ]]; then
sample_size=46350
elif [[ ${trait} == "BIP" ]]; then
sample_size=353899
elif [[ ${trait} == "BrC" ]]; then
sample_size=247173
elif [[ ${trait} == "CKD" ]]; then
sample_size=625219
elif [[ ${trait} == "CAD" ]]; then
sample_size=184305
elif [[ ${trait} == "HF" ]]; then
sample_size=428008
elif [[ ${trait} == "IBD" ]]; then
sample_size=34652
elif [[ ${trait} == "IS" ]]; then
sample_size=440328
elif [[ ${trait} == "LuC" ]]; then
sample_size=112781
elif [[ ${trait} == "MDD" ]]; then
sample_size=674452
elif [[ ${trait} == "Osteoporosis" ]]; then
sample_size=438872
elif [[ ${trait} == "OvC" ]]; then
sample_size=97898
elif [[ ${trait} == "PAD" ]]; then
sample_size=174993
elif [[ ${trait} == "PBC" ]]; then
sample_size=24510
elif [[ ${trait} == "PSC" ]]; then
sample_size=24751
elif [[ ${trait} == "PrC" ]]; then
sample_size=140306
elif [[ ${trait} == "PAH" ]]; then
sample_size=11744
elif [[ ${trait} == "RA" ]]; then
sample_size=58284
elif [[ ${trait} == "SCZ" ]]; then
sample_size=130644
elif [[ ${trait} == "SLE" ]]; then
sample_size=23210
elif [[ ${trait} == "T2D" ]]; then
sample_size=455313
else
echo "Please provide the available phenotype"
fi

param_phi=auto

input_dir="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/PRScsx/MTAG_meta/${trait}/"
out_dir="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MTAG_meta_PRS/${trait}/"

for sumstat_file in "${input_dir}"/mtag_*${trait}_*_EUR_PRScsx.txt; do
  
base_name=$(basename "${sumstat_file}")
out_name="${base_name%.txt}"

for chr in {1..22}; do

out_check="${out_dir}/${out_name}_${pop}_pst_eff_a1_b0.5_phi${param_phi}_chr${chr}.txt"
    
if [[ ! -e "${out_check}" ]]; then

echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/; python /gpfs/gibbs/pi/zhao/lx94/PRSmap/method/PRScsx/PRScsx.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} --sst_file=${sumstat_file} --n_gwas=${sample_size} --chrom=${chr} --pop=${pop} --out_dir=${out_dir} --out_name=${out_name}" >> $job_file

fi

done
done


module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/

dsq --job-file ${trait}_multi_trait_beta.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-${trait}_multi_trait_beta-$(date +%Y-%m-%d).sh
done


# Step2: Obtain beta by chr pop for each param in each trait
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

pop = "EUR"
param_ph= "auto"
trait = args[1]

out_dir   <- paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MTAG_meta_PRS/",trait,"/")

sumstat_files <- list.files(
  path       = out_dir,
  pattern    = ".*_pst_eff_a1_b0.5_phiauto_chr1.txt$", 
  full.names = TRUE
)

base_names <- sub("_pst_eff_a1_b0.5_phiauto_chr1.txt$", "", basename(sumstat_files))

for (bn in base_names) {
  
  # We'll accumulate the data in one data.table
  PRScsx_all <- data.table()
  
  # For each chromosome 1..22
  for (chr in 1:22) {
    chr_file <- file.path(out_dir, paste0(bn, "_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt"))

    # Read the chromosome-specific output
    PRScsx_pop_chr <- fread(chr_file)
    PRScsx_pop_chr <- PRScsx_pop_chr[, c(2,4,6)]
    setnames(PRScsx_pop_chr, c("SNP", "A1", pop))
    
    # Append to our full dataset
    PRScsx_all <- rbind(PRScsx_all, PRScsx_pop_chr)
  }
  
  # Write the combined data
  # The final output might be something like "<bn>_chr1to22_combined.txt"
  out_file <- file.path(out_dir, paste0(bn, ".txt"))
  
  write.table(PRScsx_all,out_file,quote=F,sep='\t',row.names=F,col.names=T)

}

vim organize_MATG_PRScsx.R

job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/organize_MATG_PRScsx.txt"
> $job_file  # Empty the job file if it already exists

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

echo "module load miniconda; conda activate r_env; Rscript /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/organize_MATG_PRScsx.R ${trait}" >> $job_file

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/

dsq --job-file organize_MATG_PRScsx.txt --partition=day --requeue --mem=20G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-organize_MATG_PRScsx-$(date +%Y-%m-%d).sh


## clean unnecessary file
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

param_phi=auto

input_dir="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/PRScsx/MTAG_meta/${trait}/"
out_dir="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MTAG_meta_PRS/${trait}/"

for sumstat_file in "${input_dir}"/mtag_*${trait}_*_EUR_PRScsx.txt; do
  
base_name=$(basename "${sumstat_file}")
out_name="${base_name%.txt}"

for chr in {1..22}; do

out_check="${out_dir}/${out_name}_${pop}_pst_eff_a1_b0.5_phi${param_phi}_chr${chr}.txt"

rm -rf ${out_check}

done
done
done
