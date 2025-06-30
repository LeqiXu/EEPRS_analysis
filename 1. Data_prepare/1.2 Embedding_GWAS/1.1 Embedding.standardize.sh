module load miniconda
conda activate r-tensorflow

# Step0: environment prepare & data prepare
## environment prepare
unlink("~/.virtualenvs/r-reticulate", recursive = TRUE, force = TRUE)

library(reticulate)

Sys.unsetenv("VIRTUALENV_NAME")
use_condaenv("r-tensorflow", required = TRUE)
py_config()

## data prepare
library(data.table)

setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/individual_embedding")

pd <- import("pandas")
train_embedding_word2vec <- pd$read_pickle("ind_eurtrain_embed_word2vec100.pkl")
train_embedding_gpticd <- pd$read_pickle("ind_eurtrain_embed_gpticd3072.pkl")

reshape_embeddings <- function(dt) {
  dt = as.data.table(dt)

  ID_vector = colnames(dt)  # Extract column names
  dt = transpose(dt)  # Transpose so individuals are rows
  colnames(dt) = paste0("Embedding_", seq_len(ncol(dt)))  # Rename embeddings
  
  dt = cbind(data.table(FID = ID_vector, IID = ID_vector),dt)
  
  return(dt)
}

train_embedding_word2vec = reshape_embeddings(train_embedding_word2vec)
train_embedding_gpticd = reshape_embeddings(train_embedding_gpticd)

write.table(train_embedding_word2vec, file=paste0("word2vec100_EUR_train_embed.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")
write.table(train_embedding_gpticd, file=paste0("gpticd3072_EUR_train_embed.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

# Step1: plot original embedding distribution
library(data.table)
library(ggplot2)

setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/individual_embedding")

train_embedding_word2vec = fread(paste0("word2vec100_EUR_train_embed.txt"))
train_embedding_gpticd = fread(paste0("gpticd3072_EUR_train_embed.txt"))

train_embedding_word2vec = na.omit(train_embedding_word2vec)
train_embedding_gpticd = na.omit(train_embedding_gpticd)

ggplot(train_embedding_word2vec, aes(sample = Embedding_1)) +
  stat_qq() + 
  stat_qq_line() +
  labs(title = "QQ-Plot of Word2Vec Embedding_1",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()

ggplot(train_embedding_gpticd, aes(sample = Embedding_1)) +
  stat_qq() + 
  stat_qq_line() +
  labs(title = "QQ-Plot of GPT-based Embedding_1",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()

# Step2: normalize the phenotype
inverse_normal_transform <- function(x) {
  rank_x <- rank(x, ties.method = "average")  # Compute ranks
  p <- (rank_x - 0.5) / length(x)             # Convert to percentiles
  qnorm(p)                                    # Apply inverse normal transformation
}

apply_inverse_normal_transform <- function(dt) {
  embedding_cols <- setdiff(names(dt), c("FID", "IID"))
  dt[, (embedding_cols) := lapply(.SD, inverse_normal_transform), .SDcols = embedding_cols]
  
  return(dt)
}

train_embedding_word2vec_transformed <- apply_inverse_normal_transform(train_embedding_word2vec)
train_embedding_gpticd_transformed <- apply_inverse_normal_transform(train_embedding_gpticd)

# Step3: plot standardized embedding distribution
## first embedding plot
ggplot(train_embedding_word2vec_transformed, aes(sample = Embedding_1)) +
  stat_qq() + 
  stat_qq_line() +
  labs(title = "QQ-Plot of Word2Vec Standardized Embedding_1",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()

write.table(train_embedding_word2vec_transformed, file=paste0("word2vec100_EUR_train_embed_normalize.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

ggplot(train_embedding_gpticd_transformed, aes(sample = Embedding_1)) +
  stat_qq() + 
  stat_qq_line() +
  labs(title = "QQ-Plot of GPT-based Standardized Embedding_1",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()

write.table(train_embedding_gpticd_transformed, file=paste0("gpticd3072_EUR_train_embed_normalize.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")
