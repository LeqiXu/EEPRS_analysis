## Step1: Linear combination with different snplists
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/subsample_prune_gpt_ica_trait_PRS_linear_weights.txt"
> $job_file  # Empty the job file if it already exists

pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

trait_significant_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/trait_embedding_corr/trait_significant_gpt_ica_embeddings.txt"

for trait in ${trait_list}; do

selection_criterion="NO"
non_negative_weights="FALSE"
weight_name="linear_weights_approxTRUE"

## only considers embeddings with significant correlations
embeddings_list=$(grep -w "^${trait}" "${trait_significant_file}"  | cut -f2)
IFS=' ' read -r -a embedding_array <<< "${embeddings_list}"

embedding_files_array=()
for emb_i in "${embedding_array[@]}"; do
  embedding_files_array+=("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/gpticd3072_ICA67_train_EUR_UKB_V${emb_i}_prune_${pop}_PRScsx.txt")
done
embedding_files_str=$(IFS=,; echo "${embedding_files_array[*]}")

if [[ ! -z "$embedding_files_str" ]]; then
for rpt in {1..4}; do

## trait-specific prs
trait_prs_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_subsample_PRS/${trait}_prune_${pop}_train_PRScsx_repeat${rpt}.txt"

## final prs_files
prs_files="${trait_prs_file},${embedding_files_str}"

output_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Final_weight/${trait}_gpticd3072_ICA67_MIXPRS_weight_repeat${rpt}_${pop}_${weight_name}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/subsample_data/MIX/${trait}_prune_${pop}_tune_MIX_approxTRUE_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=${prs_files} --selection_criterion=${selection_criterion} --non_negative_weights=${non_negative_weights} --out_dir=result/Final_weight --out_name=${trait}_gpticd3072_ICA67_MIXPRS_weight_repeat${rpt}" >> $job_file
fi

done
fi

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/subsample_prune_gpt_ica_trait_PRS_linear_weights.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_gpt_ica_trait_PRS_linear_weights-$(date +%Y-%m-%d).sh



## Step2: Obtain final MIXPRS weight
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/subsample_prune_gpt_ica_trait_MIXPRS_weight.txt"
> $job_file  # Empty the job file if it already exists

pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

trait_significant_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/trait_embedding_corr/trait_significant_gpt_ica_embeddings.txt"

for trait in ${trait_list}; do

weight_name="linear_weights_approxTRUE"

sst_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/MIX/${trait}_${pop}_all_MIX.txt"

## only considers embeddings with significant correlations
embeddings_list=$(grep -w "^${trait}" "${trait_significant_file}"  | cut -f2)
IFS=' ' read -r -a embedding_array <<< "${embeddings_list}"

embedding_files_array=()
for emb_i in "${embedding_array[@]}"; do
  embedding_files_array+=("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/gpticd3072_ICA67_train_EUR_UKB_V${emb_i}_prune_${pop}_PRScsx.txt")
done
embedding_files_str=$(IFS=,; echo "${embedding_files_array[*]}")

## trait-specific prs
trait_prs_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_PRS/${trait}_${pop}_PRScsx.txt"

## final prs_files
prs_files="${trait_prs_file},${embedding_files_str}"

## weight_files
weight_file1="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Final_weight/${trait}_gpticd3072_ICA67_MIXPRS_weight_repeat1_${pop}_${weight_name}.txt"
weight_file2="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Final_weight/${trait}_gpticd3072_ICA67_MIXPRS_weight_repeat2_${pop}_${weight_name}.txt"
weight_file3="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Final_weight/${trait}_gpticd3072_ICA67_MIXPRS_weight_repeat3_${pop}_${weight_name}.txt"
weight_file4="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Final_weight/${trait}_gpticd3072_ICA67_MIXPRS_weight_repeat4_${pop}_${weight_name}.txt"

output_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/${trait}_gpticd3072_ICA67_${pop}_MIXPRS.txt"

if [[ ! -z "$embedding_files_str" ]]; then
if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_final_combine.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=${sst_file} --pop=${pop} --prs_beta_file=${prs_files} --weight_file=${weight_file1},${weight_file2},${weight_file3},${weight_file4} --out_dir=result/MIXPRS --out_name=${trait}_gpticd3072_ICA67" >> $job_file
fi
else
    cp /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_PRS/${trait}_${pop}_PRScsx.txt /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/${trait}_gpticd3072_ICA67_${pop}_MIXPRS.txt
fi

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/subsample_prune_gpt_ica_trait_MIXPRS_weight.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_gpt_ica_trait_MIXPRS_weight-$(date +%Y-%m-%d).sh
