## 0. PRScs
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_PRS/UKB_all_${trait}_PRScsx_prs_${pop}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=day
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_PRScsx_prs_${pop}
#SBATCH --output=out_PRS_${trait}_PRScsx_prs_${pop}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_PRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop}_hm3 \
--double-id \
--threads 1 \
--score ${trait}_${pop}_PRScsx.txt \
--out UKB_all_${trait}_PRScsx_prs_${pop}

EOT
fi

done

## 1. Multi trait PRScs
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/multi_trait_score.txt"
> $job_file  # Empty the job file if it already exists

pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

input_dir="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MTAG_meta_PRS/${trait}/"
out_dir="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MTAG_meta_PRS/${trait}/"

for prs_beta_file in "${input_dir}"/mtag_*${trait}_*_EUR_PRScsx_EUR.txt; do
  
base_name=$(basename "${prs_beta_file}")
out_name="${base_name%.txt}"

out_check="${out_dir}/UKB_all_${out_name}.sscore"
    
if [[ ! -e "${out_check}" ]]; then

echo "module load PLINK/2; cd ${out_dir}; plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop}_hm3 --double-id --threads 1 --score ${prs_beta_file} --out UKB_all_${out_name}" >> $job_file

fi

done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/

dsq --job-file multi_trait_score.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-multi_trait_score-$(date +%Y-%m-%d).sh
