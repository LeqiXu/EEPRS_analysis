## Step1: LDSC format data
## PCA
library(data.table)

train_type = "train"

for (i in c(1:30)){
    # sumstat for clean format
    sumstat_clean = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/clean/word2vec100_PCA30_",train_type,"_EUR_UKB_PC",i,"_clean.txt"))

    # sumstat for ldsc format
    sumstat_ldsc = sumstat_clean[,c("SNP","A1","A2","N","P","Z")]
    write.table(sumstat_ldsc, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/ldsc/word2vec100_PCA30_",train_type,"_EUR_UKB_PC",i,"_ldsc.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}

## ICA
library(data.table)

train_type = "train"

for (i in c(1:30)){
    # sumstat for clean format
    sumstat_clean = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/clean/word2vec100_ICA30_",train_type,"_EUR_UKB_V",i,"_clean.txt"))

    # sumstat for ldsc format
    sumstat_ldsc = sumstat_clean[,c("SNP","A1","A2","N","P","Z")]
    write.table(sumstat_ldsc, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/ldsc/word2vec100_ICA30_",train_type,"_EUR_UKB_V",i,"_ldsc.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}


## Step2: run LDSC munge
## PCA
train_type=train
for i in {1..30}; do

python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/munge_sumstats.py \
--sumstats /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/ldsc/word2vec100_PCA30_${train_type}_EUR_UKB_PC${i}_ldsc.txt \
--out /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/word2vec100_PCA30_${train_type}_EUR_UKB_PC${i}_ldsc

done

## ICA
train_type=train
for i in {1..30}; do

python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/munge_sumstats.py \
--sumstats /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/ldsc/word2vec100_ICA30_${train_type}_EUR_UKB_V${i}_ldsc.txt \
--out /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/word2vec100_ICA30_${train_type}_EUR_UKB_V${i}_ldsc

done

## Step3: calculate heritability
## PCA
train_type=train
for i in {1..30}; do

python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/ldsc.py \
--h2 /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/word2vec100_PCA30_${train_type}_EUR_UKB_PC${i}_ldsc.sumstats.gz \
--ref-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ \
--w-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ \
--out /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/word2vec100_PCA30_${train_type}_EUR_UKB_PC${i}_h2

done

## ICA
train_type=train
for i in {1..30}; do

python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/ldsc.py \
--h2 /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/word2vec100_ICA30_${train_type}_EUR_UKB_V${i}_ldsc.sumstats.gz \
--ref-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ \
--w-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ \
--out /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/word2vec100_ICA30_${train_type}_EUR_UKB_V${i}_h2

done


## Step4: Obtain h2 table
## -----------------------------
## 1. PCA: read logs, store results
## -----------------------------
setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/")

get_embedding_number_pca <- function(filename) {
  sub("^word2vec100_PCA30_train_EUR_UKB_PC([0-9]+)_h2\\.log$", "\\1", filename)
}

files_pca <- list.files(
  pattern = "^word2vec100_PCA30_train_EUR_UKB_PC[0-9]+_h2\\.log$",
  full.names = FALSE
)

results_pca <- data.frame(
  file = character(),
  embedding = integer(),
  h2 = numeric(),
  h2_se = numeric(),
  stringsAsFactors = FALSE
)

h2_pattern <- "Total Observed scale h2:\\s+([-0-9.]+) \\(([-0-9.]+)\\)"

for (f in files_pca) {
  lines <- readLines(f, warn = FALSE)
  line_of_interest <- grep(h2_pattern, lines, value = TRUE)
  
  if (length(line_of_interest) == 1) {
    match_info <- regexec(h2_pattern, line_of_interest)
    captures   <- regmatches(line_of_interest, match_info)[[1]]
    h2_val     <- as.numeric(captures[2]) 
    h2_se_val  <- as.numeric(captures[3]) 
  } else {
    h2_val     <- NA_real_
    h2_se_val  <- NA_real_
  }
  
  emb_num_str <- get_embedding_number_pca(basename(f))
  emb_num     <- as.integer(emb_num_str)
  
  results_pca <- rbind(
    results_pca,
    data.frame(
      file = f,
      embedding = emb_num,
      h2 = h2_val,
      h2_se = h2_se_val,
      stringsAsFactors = FALSE
    )
  )
}

## -----------------------------
## 2. ICA: read logs, store results
## -----------------------------
get_embedding_number_ica <- function(filename) {
  sub("^word2vec100_ICA30_train_EUR_UKB_V([0-9]+)_h2\\.log$", "\\1", filename)
}

files_ica <- list.files(
  pattern = "^word2vec100_ICA30_train_EUR_UKB_V[0-9]+_h2\\.log$",
  full.names = FALSE
)

results_ica <- data.frame(
  file = character(),
  embedding = integer(),
  h2 = numeric(),
  h2_se = numeric(),
  stringsAsFactors = FALSE
)

for (f in files_ica) {
  lines <- readLines(f, warn = FALSE)
  line_of_interest <- grep(h2_pattern, lines, value = TRUE)
  
  if (length(line_of_interest) == 1) {
    match_info <- regexec(h2_pattern, line_of_interest)
    captures   <- regmatches(line_of_interest, match_info)[[1]]
    h2_val     <- as.numeric(captures[2])
    h2_se_val  <- as.numeric(captures[3])
  } else {
    h2_val     <- NA_real_
    h2_se_val  <- NA_real_
  }
  
  emb_num_str <- get_embedding_number_ica(basename(f))
  emb_num     <- as.integer(emb_num_str)
  
  results_ica <- rbind(
    results_ica,
    data.frame(
      file = f,
      embedding = emb_num,
      h2 = h2_val,
      h2_se = h2_se_val,
      stringsAsFactors = FALSE
    )
  )
}

## -----------------------------
## 3. Combine the results
## -----------------------------
library(dplyr)

results_pca <- results_pca %>% mutate(method = "PCA")
results_ica <- results_ica %>% mutate(method = "ICA")

results_pca$h2_lower <- results_pca$h2 - 1.96 * results_pca$h2_se
results_ica$h2_lower <- results_ica$h2 - 1.96 * results_ica$h2_se

results_pca <- results_pca %>% 
  arrange(desc(h2)) %>%
  mutate(embedding_ordered = factor(paste0("PCA_",as.character(embedding)), levels = unique(paste0("PCA_",as.character(embedding)))))

results_ica <- results_ica %>% 
  arrange(desc(h2)) %>%
  mutate(embedding_ordered = factor(paste0("ICA_",as.character(embedding)), levels = unique(paste0("ICA_",as.character(embedding)))))


results_pca$embedding_approach = "Word2Vec_PCA"
results_pca

write.table(results_pca, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/UKB_word2vec100_pca_train_h2.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

results_ica$embedding_approach = "Word2Vec_ICA"
results_ica

write.table(results_ica, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/UKB_word2vec100_ica_train_h2.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

## -----------------------------
## 4. Plot with ggplot2
## -----------------------------
library(ggplot2)
library(cowplot)
library(ggpubr)

results_pca <- results_pca %>%
  mutate(above_zero = (h2 - 1.96 * h2_se) > 0)
results_pca$above_zero <- factor(results_pca$above_zero, levels = c("TRUE","FALSE"))

results_ica <- results_ica %>%
  mutate(above_zero = (h2 - 1.96 * h2_se) > 0)
results_ica$above_zero <- factor(results_ica$above_zero, levels = c("TRUE","FALSE"))

p_pca_noy <- ggplot(results_pca, aes(x = embedding_ordered, y = h2, color = above_zero)) +
  geom_point() +
  geom_errorbar(aes(ymin = h2 - 1.96*h2_se, ymax = h2 + 1.96*h2_se), width=0.3) +
  facet_wrap(~method, scales="free_x") +
  labs(x="Embedding Number", y=NULL, title="Heritability Estimates (Word2Vec_PCA)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none") + 
  ylim(c(-0.01,0.1))

p_ica_noy <- ggplot(results_ica, aes(x=embedding_ordered, y=h2, color=above_zero)) +
  geom_point() +
  geom_errorbar(aes(ymin = h2 - 1.96*h2_se, ymax = h2 + 1.96*h2_se), width=0.3) +
  facet_wrap(~method, scales="free_x") +
  labs(x="Embedding Number", y=NULL, title="Heritability Estimates (Word2Vec_ICA)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none") + 
  ylim(c(-0.01,0.1))

p_pca_legend <- p_pca_noy + theme(legend.position="right")
shared_legend <- get_legend(p_pca_legend)

two_plots <- ggarrange(
  p_pca_noy, p_ica_noy, 
  ncol=2, nrow=1, 
  align="v",
  labels = c("A", "B")
)

combined_plots <- ggarrange(
  two_plots, shared_legend,
  ncol=2,
  widths = c(2, 0.2)  # adjust legend space as needed
)

final_plot <- annotate_figure(
  combined_plots,
  left = text_grob("Estimated Heritability (h2)", rot = 90)
)

print(final_plot)

pca_sig_embedding <- results_pca$embedding[results_pca$h2_lower > 0]
pca_sig_embedding <- c(1,2,7,9,4,3,5,18,11,8,20,15,6,10,16,17,14,13,21,25,26,12,27,24)

ica_sig_embedding <- results_ica$embedding[results_ica$h2_lower > 0]
ica_sig_embedding <- c(25,19,12,2,15,11,22,27,4,14,24,5,8,1,10,20,16,13,9,6,29,7)

## Width: 1800 Height: 800