## Step1: hm3 PRScsx beta
# Calculate PRS
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/trait_PRScs.txt"
> $job_file  # Empty the job file if it already exists

pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

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

for chr in {1..22}; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_PRS/${trait}_${pop}_PRScsx_${pop}_pst_eff_a1_b0.5_phi${param_phi}_chr${chr}.txt" ]]; then

echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/; python /gpfs/gibbs/pi/zhao/lx94/PRSmap/method/PRScsx/PRScsx.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} --sst_file=data/summary_data/PRScsx/${trait}_${pop}_all_PRScsx.txt --n_gwas=${sample_size} --chrom=${chr} --pop=${pop} --out_dir=result/Trait_PRS/ --out_name=${trait}_${pop}_PRScsx" >> $job_file

fi
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/trait_PRScs.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=1-00:00:00 --mail-type=ALL
sbatch dsq-trait_PRScs-$(date +%Y-%m-%d).sh


## Step2: Obtain beta by chr pop for each param in each trait
library(data.table)

pop="EUR"

trait_list = c("HDL","LDL","TC","logTG","BMI","WHR","AgeSmk","SmkI","SmkC","CigDay","DrnkWk","Glu2h","HbA1c","eGFR","logCp",
                "AD","Angina","Asthma","AF","ADHD","ASD","BIP","BrC","CKD","CAD","HF","IBD","IS","LuC",
                "MDD","Osteoporosis","OvC","PAD","PBC","PSC","PrC","PAH","RA","SCZ","SLE","T2D")

for (trait in trait_list){

PRScsx_all <- data.table()
for (chr in c(1:22)){

PRScsx_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_PRS/",trait,"_",pop,"_PRScsx_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))

PRScsx_pop_chr <- PRScsx_pop_chr[,c(2,4,6)]
names(PRScsx_pop_chr) = c("SNP","A1",pop)

PRScsx_all = rbind(PRScsx_all,PRScsx_pop_chr)
    
}

write.table(PRScsx_all,paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_PRS/",trait,"_",pop,"_PRScsx.txt"),quote=F,sep='\t',row.names=F,col.names=T)
}


## Step3: Clean chr file
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do
for chr in {1..22}; do

param_phi=auto

rm -rf /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_PRS/${trait}_${pop}_PRScsx_${pop}_pst_eff_a1_b0.5_phi${param_phi}_chr${chr}.txt

done
done