## Step1: hm3 PRScsx beta
# Calculate PRS
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/trait_subsample_PRScs.txt"
> $job_file  # Empty the job file if it already exists

pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

# sample size
if [[ ${trait} == "HDL" ]]; then
  # original sample_size=950886
  sample_size=713165    # 3/4 * 950886 = 713164.5 → 713165
elif [[ ${trait} == "LDL" ]]; then
  # original sample_size=950886
  sample_size=713165
elif [[ ${trait} == "TC" ]]; then
  # original sample_size=950886
  sample_size=713165 
elif [[ ${trait} == "logTG" ]]; then
  # original sample_size=950886
  sample_size=713165 
elif [[ ${trait} == "BMI" ]]; then
  # original sample_size=339224
  sample_size=254418    # 3/4 * 339224 = 254418 (exact)  
elif [[ ${trait} == "WHR" ]]; then
  # original sample_size=224459
  sample_size=168344    # 3/4 * 224459 = 168344.25 → 168344  
elif [[ ${trait} == "AgeSmk" ]]; then
  # original sample_size=175835
  sample_size=131876    # 3/4 * 175835 = 131876.25 → 131876  
elif [[ ${trait} == "SmkI" ]]; then
  # original sample_size=357235
  sample_size=267926    # 3/4 * 357235 = 267926.25 → 267926  
elif [[ ${trait} == "SmkC" ]]; then
  # original sample_size=188701
  sample_size=141526    # 3/4 * 188701 = 141525.75 → 141526
elif [[ ${trait} == "CigDay" ]]; then
  # original sample_size=183196
  sample_size=137397    # 3/4 * 183196 = 137397 (exact)
elif [[ ${trait} == "DrnkWk" ]]; then
  # original sample_size=304322
  sample_size=228242    # 3/4 * 304322 = 228241.5 → 228242
elif [[ ${trait} == "Glu2h" ]]; then
  # original sample_size=63396
  sample_size=47547     # 3/4 * 63396 = 47547 (exact)
elif [[ ${trait} == "HbA1c" ]]; then
  # original sample_size=146806
  sample_size=110105    # 3/4 * 146806 = 110104.5 → 110105
elif [[ ${trait} == "eGFR" ]]; then
  # original sample_size=567460
  sample_size=425595    # 3/4 * 567460 = 425595 (exact)
elif [[ ${trait} == "logCp" ]]; then
  # original sample_size=204402
  sample_size=153302    # 3/4 * 204402 = 153301.5 → 153302
elif [[ ${trait} == "AD" ]]; then
  # original sample_size=63926
  sample_size=47945     # 3/4 * 63926 = 47944.5 → 47945
elif [[ ${trait} == "Angina" ]]; then
  # original sample_size=418385
  sample_size=313789    # 3/4 * 418385 = 313788.75 → 313789
elif [[ ${trait} == "Asthma" ]]; then
  # original sample_size=127669
  sample_size=95752     # 3/4 * 127669 = 95751.75 → 95752
elif [[ ${trait} == "AF" ]]; then
  # original sample_size=635097
  sample_size=476323    # 3/4 * 635097 = 476322.75 → 476323
elif [[ ${trait} == "ADHD" ]]; then
  # original sample_size=225534
  sample_size=169151    # 3/4 * 225534 = 169150.5 → 169151
elif [[ ${trait} == "ASD" ]]; then
  # original sample_size=46350
  sample_size=34763     # 3/4 * 46350 = 34762.5 → 34763
elif [[ ${trait} == "BIP" ]]; then
  # original sample_size=353899
  sample_size=265424    # 3/4 * 353899 = 265424.25 → 265424
elif [[ ${trait} == "BrC" ]]; then
  # original sample_size=247173
  sample_size=185380    # 3/4 * 247173 = 185379.75 → 185380
elif [[ ${trait} == "CKD" ]]; then
  # original sample_size=625219
  sample_size=468914    # 3/4 * 625219 = 468914.25 → 468914
elif [[ ${trait} == "CAD" ]]; then
  # original sample_size=184305
  sample_size=138229    # 3/4 * 184305 = 138228.75 → 138229
elif [[ ${trait} == "HF" ]]; then
  # original sample_size=428008
  sample_size=321006    # 3/4 * 428008 = 321006 (exact)
elif [[ ${trait} == "IBD" ]]; then
  # original sample_size=34652
  sample_size=25989     # 3/4 * 34652 = 25989 (exact)
elif [[ ${trait} == "IS" ]]; then
  # original sample_size=440328
  sample_size=330246    # 3/4 * 440328 = 330246 (exact)
elif [[ ${trait} == "LuC" ]]; then
  # original sample_size=112781
  sample_size=84586     # 3/4 * 112781 = 84585.75 → 84586
elif [[ ${trait} == "MDD" ]]; then
  # original sample_size=674452
  sample_size=505839    # 3/4 * 674452 = 505839 (exact)
elif [[ ${trait} == "Osteoporosis" ]]; then
  # original sample_size=438872
  sample_size=329154    # 3/4 * 438872 = 329154 (exact)
elif [[ ${trait} == "OvC" ]]; then
  # original sample_size=97898
  sample_size=73424     # 3/4 * 97898 = 73423.5 → 73424
elif [[ ${trait} == "PAD" ]]; then
  # original sample_size=174993
  sample_size=131245    # 3/4 * 174993 = 131244.75 → 131245
elif [[ ${trait} == "PBC" ]]; then
  # original sample_size=24510
  sample_size=18383     # 3/4 * 24510 = 18382.5 → 18383
elif [[ ${trait} == "PSC" ]]; then
  # original sample_size=24751
  sample_size=18563     # 3/4 * 24751 = 18563.25 → 18563
elif [[ ${trait} == "PrC" ]]; then
  # original sample_size=140306
  sample_size=105230    # 3/4 * 140306 = 105229.5 → 105230
elif [[ ${trait} == "PAH" ]]; then
  # original sample_size=11744
  sample_size=8808      # 3/4 * 11744 = 8808 (exact)
elif [[ ${trait} == "RA" ]]; then
  # original sample_size=58284
  sample_size=43713     # 3/4 * 58284 = 43713 (exact)
elif [[ ${trait} == "SCZ" ]]; then
  # original sample_size=130644
  sample_size=97983     # 3/4 * 130644 = 97983 (exact)
elif [[ ${trait} == "SLE" ]]; then
  # original sample_size=23210
  sample_size=17408     # 3/4 * 23210 = 17407.5 → 17408
elif [[ ${trait} == "T2D" ]]; then
  # original sample_size=455313
  sample_size=341485    # 3/4 * 455313 = 341484.75 → 341485
else
  echo "Please provide an available phenotype."
fi

param_phi=auto

for rpt in {1..4}; do
for chr in {1..22}; do
    
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_subsample_PRS/${trait}_prune_${pop}_train_PRScsx_repeat${rpt}_${pop}_pst_eff_a1_b0.5_phi${param_phi}_chr${chr}.txt" ]]; then
                
echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/; python /gpfs/gibbs/pi/zhao/lx94/PRSmap/method/PRScsx/PRScsx.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} --sst_file=data/subsample_data/PRScsx/${trait}_prune_${pop}_train_PRScsx_repeat${rpt}.txt --n_gwas=${sample_size} --chrom=${chr} --pop=${pop} --out_dir=result/Trait_subsample_PRS/ --out_name=${trait}_prune_${pop}_train_PRScsx_repeat${rpt}" >> $job_file
    
fi

done
done

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/trait_subsample_PRScs.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=1-00:00:00 --mail-type=ALL
sbatch dsq-trait_subsample_PRScs-$(date +%Y-%m-%d).sh


## Step2: Obtain beta by chr pop for each param in each trait for each repeat
library(data.table)

pop="EUR"

trait_list = c("HDL","LDL","TC","logTG","BMI","WHR","AgeSmk","SmkI","SmkC","CigDay","DrnkWk","Glu2h","HbA1c","eGFR","logCp",
                "AD","Angina","Asthma","AF","ADHD","ASD","BIP","BrC","CKD","CAD","HF","IBD","IS","LuC",
                "MDD","Osteoporosis","OvC","PAD","PBC","PSC","PrC","PAH","RA","SCZ","SLE","T2D")

for (trait in trait_list){
for (rpt in c(1:4)){

PRScsx_all <- data.table()

for (chr in c(1:22)){

PRScsx_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_subsample_PRS/",trait,"_prune_",pop,"_train_PRScsx_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))

PRScsx_pop_chr <- PRScsx_pop_chr[,c(2,4,6)]
names(PRScsx_pop_chr) = c("SNP","A1",pop)

PRScsx_all = rbind(PRScsx_all,PRScsx_pop_chr)
    
}

write.table(PRScsx_all,paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_subsample_PRS/",trait,"_prune_",pop,"_train_PRScsx_repeat",rpt,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}
}


## Step3: Clean chr file
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do
for rpt in {1..4}; do
for chr in {1..22}; do

param_phi=auto

rm -rf /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_subsample_PRS/${trait}_prune_${pop}_train_PRScsx_repeat${rpt}_${pop}_pst_eff_a1_b0.5_phi${param_phi}_chr${chr}.txt

done
done
done