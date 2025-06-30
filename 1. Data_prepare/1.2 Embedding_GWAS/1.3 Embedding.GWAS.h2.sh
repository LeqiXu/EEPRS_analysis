## Word2Vec
## Step1: LDSC format data
library(data.table)

for (train_type in c("train")){

for (i in c(1:100)){
    # sumstat for clean format
    sumstat_clean = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/clean/word2vec100_",train_type,"_EUR_UKB_Embedding",i,"_clean.txt"))

    # sumstat for ldsc format
    sumstat_ldsc = sumstat_clean[,c("SNP","A1","A2","N","P","Z")]
    write.table(sumstat_ldsc, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/ldsc/word2vec100_",train_type,"_EUR_UKB_Embedding",i,"_ldsc.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}
}


## Step2: run LDSC munge
for train_type in train; do
for i in {1..100}; do

python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/munge_sumstats.py \
--sumstats /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/ldsc/word2vec100_${train_type}_EUR_UKB_Embedding${i}_ldsc.txt \
--out /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/word2vec100_${train_type}_EUR_UKB_Embedding${i}_ldsc

done
done

## Step3: calculate heritability
for train_type in train; do
for i in {1..100}; do

python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/ldsc.py \
--h2 /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/word2vec100_${train_type}_EUR_UKB_Embedding${i}_ldsc.sumstats.gz \
--ref-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ \
--w-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ \
--out /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/word2vec100_${train_type}_EUR_UKB_Embedding${i}_h2

done
done


## Step4: Obtain h2 table
setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/")

get_embedding_number <- function(filename) {
  sub("^word2vec100_train_EUR_UKB_Embedding([0-9]+)_h2\\.log$", "\\1", filename)
}

files <- list.files(
  pattern = "^word2vec100_train_EUR_UKB_Embedding[0-9]+_h2\\.log$",
  full.names = FALSE
)

results <- data.frame(
  file = character(),
  embedding = integer(),
  h2 = numeric(),
  h2_se = numeric(),
  stringsAsFactors = FALSE
)

# Regex pattern for extracting h2 and SE from the log line
h2_pattern <- "Total Observed scale h2:\\s+([-0-9.]+) \\(([-0-9.]+)\\)"

for (f in files) {
  # Read the file contents
  lines <- readLines(f, warn = FALSE)
  
  # Locate line with "Total Observed scale h2: X (Y)"
  line_of_interest <- grep(h2_pattern, lines, value = TRUE)
  
  if (length(line_of_interest) == 1) {
    match_info <- regexec(h2_pattern, line_of_interest)
    captures   <- regmatches(line_of_interest, match_info)[[1]]
    h2_val     <- as.numeric(captures[2])  # the numeric for h2
    h2_se_val  <- as.numeric(captures[3])  # the numeric for h2_se
  } else {
    h2_val     <- NA_real_
    h2_se_val  <- NA_real_
  }
  
  # Parse the embedding number from the filename
  emb_num_str <- get_embedding_number(basename(f))
  emb_num     <- as.integer(emb_num_str)  # convert "99" to 99
  
  # Add one row to results
  results <- rbind(results, data.frame(
    file = f,
    embedding = emb_num,
    h2 = h2_val,
    h2_se = h2_se_val,
    stringsAsFactors = FALSE
  ))
}

results$embedding_approach = "Word2Vec"
results

write.table(results, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/UKB_word2vec100_train_h2.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

## Step5: Plot heritability
library(ggplot2)

results$h2_lower <- results$h2 - 1.96 * results$h2_se
results <- results %>% mutate(above_zero = (h2 - 1.96 * h2_se) > 0)

results <- results %>% 
  arrange(desc(h2)) %>%
  mutate(embedding_ordered = factor(as.character(embedding), levels = unique(as.character(embedding))))

# Create a ggplot
ggplot(results, aes(x = embedding_ordered, y = h2, color=above_zero)) +
  geom_point() +
  geom_errorbar(aes(ymin = h2 - 1.96 * h2_se, ymax = h2 + 1.96 * h2_se),
                width = 0.3) +
  labs(
    x = "Embedding Number",
    y = NULL,
    title = "Heritability Estimates (Word2Vec)"
  ) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none") + 
  ylim(c(-0.01,0.1))

## Width: 1800 Height: 800



## GPT
## Step1: LDSC format data
library(data.table)

for (train_type in c("train")){

for (i in c(1:3072)){
    # sumstat for clean format
    sumstat_clean = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/clean/gpticd3072_",train_type,"_EUR_UKB_Embedding",i,"_clean.txt"))

    # sumstat for ldsc format
    sumstat_ldsc = sumstat_clean[,c("SNP","A1","A2","N","P","Z")]
    write.table(sumstat_ldsc, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/ldsc/gpticd3072_",train_type,"_EUR_UKB_Embedding",i,"_ldsc.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}
}


## Step2: run LDSC munge
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_munge.txt"
> $job_file  # Empty the job file if it already exists

train_type=train

for i in {1..3072}; do

# Construct the sumstats files for embeddings i and j (comma-separated)
if [[ ! -e /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/gpticd3072_${train_type}_EUR_UKB_Embedding${i}_ldsc.sumstats.gz ]]; then
echo "module load miniconda; conda activate ldsc; python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/munge_sumstats.py --sumstats /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/ldsc/gpticd3072_${train_type}_EUR_UKB_Embedding${i}_ldsc.txt --out /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/gpticd3072_${train_type}_EUR_UKB_Embedding${i}_ldsc" >> $job_file
fi

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_munge.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=1-00:00:00 --mail-type=ALL
sbatch dsq-gpticd3072_munge-$(date +%Y-%m-%d).sh

## Step3: calculate heritability
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_h2.txt"
> $job_file  # Empty the job file if it already exists

train_type=train

for i in {1..3072}; do

# Construct the sumstats files for embeddings i and j (comma-separated)
if [[ ! -e /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/gpticd3072_${train_type}_EUR_UKB_Embedding${i}_h2.log ]]; then
echo "module load miniconda; conda activate ldsc; python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/ldsc.py --h2 /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/gpticd3072_${train_type}_EUR_UKB_Embedding${i}_ldsc.sumstats.gz --ref-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ --w-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ --out /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/gpticd3072_${train_type}_EUR_UKB_Embedding${i}_h2" >> $job_file
fi
    
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_h2.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=1-00:00:00 --mail-type=ALL
sbatch dsq-gpticd3072_h2-$(date +%Y-%m-%d).sh


## Step4: Obtain h2 table
setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/")

get_embedding_number <- function(filename) {
  sub("^gpticd3072_train_EUR_UKB_Embedding([0-9]+)_h2\\.log$", "\\1", filename)
}

files <- list.files(
  pattern = "^gpticd3072_train_EUR_UKB_Embedding[0-9]+_h2\\.log$",
  full.names = FALSE
)

results <- data.frame(
  file = character(),
  embedding = integer(),
  h2 = numeric(),
  h2_se = numeric(),
  stringsAsFactors = FALSE
)

# Regex pattern for extracting h2 and SE from the log line
h2_pattern <- "Total Observed scale h2:\\s+([-0-9.]+) \\(([-0-9.]+)\\)"

for (f in files) {
  # Read the file contents
  lines <- readLines(f, warn = FALSE)
  
  # Locate line with "Total Observed scale h2: X (Y)"
  line_of_interest <- grep(h2_pattern, lines, value = TRUE)
  
  if (length(line_of_interest) == 1) {
    match_info <- regexec(h2_pattern, line_of_interest)
    captures   <- regmatches(line_of_interest, match_info)[[1]]
    h2_val     <- as.numeric(captures[2])  # the numeric for h2
    h2_se_val  <- as.numeric(captures[3])  # the numeric for h2_se
  } else {
    h2_val     <- NA_real_
    h2_se_val  <- NA_real_
  }
  
  # Parse the embedding number from the filename
  emb_num_str <- get_embedding_number(basename(f))
  emb_num     <- as.integer(emb_num_str)  # convert "99" to 99
  
  # Add one row to results
  results <- rbind(results, data.frame(
    file = f,
    embedding = emb_num,
    h2 = h2_val,
    h2_se = h2_se_val,
    stringsAsFactors = FALSE
  ))
}

results$embedding_approach = "GPT"
results

write.table(results, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/UKB_gpticd3072_train_h2.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

## Step5: Plot heritability
library(ggplot2)
library(dplyr)

results$h2_lower <- results$h2 - 1.96 * results$h2_se
results <- results %>% mutate(above_zero = (h2 - 1.96 * h2_se) > 0)

results <- results %>% 
  arrange(desc(h2)) %>%
  mutate(embedding_ordered = factor(as.character(embedding), levels = unique(as.character(embedding))))

# Create a ggplot
## distribution
ggplot(results, aes(x = h2, fill = above_zero)) +
  geom_histogram(binwidth = 0.005, color = "black") +
  scale_fill_manual(values = c("FALSE" = "gray80", "TRUE" = "skyblue")) +
  labs(
    title = "Distribution of Heritability Estimates",
    x = "Heritability (h2)",
    y = "Count"
  ) +
  theme_bw()

## top100
library(dplyr)
library(ggplot2)

top_n_to_plot <- 100
results_top <- results %>% 
  arrange(desc(h2)) %>% 
  slice_head(n = top_n_to_plot) %>%
  mutate(embedding = factor(embedding, levels = embedding))

ggplot(results_top, aes(x = embedding, y = h2)) +
  geom_point(color = "blue") +
  geom_errorbar(
    aes(ymin = h2 - 1.96 * h2_se, ymax = h2 + 1.96 * h2_se),
    width = 0.3
  ) +
  labs(
    x = "Embedding",
    y = "Heritability (h2)",
    title = paste("Top", top_n_to_plot, "Heritability Estimates")
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Width: 1800 Height: 800
