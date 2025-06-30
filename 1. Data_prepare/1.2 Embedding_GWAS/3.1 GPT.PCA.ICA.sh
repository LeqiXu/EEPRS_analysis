library(data.table)
library(fastICA)
library(ggplot2)

setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/individual_embedding")

## Prepare data
train_embedding_transformed = fread(paste0("gpticd3072_EUR_train_embed_normalize.txt"))
train_embedding_matrix <- as.matrix(train_embedding_transformed[,-c("FID","IID")])

## Run PCA to find the most important embeddings
pca_res <- prcomp(train_embedding_matrix, center=FALSE, scale.=FALSE)
pve <- pca_res$sdev^2 / sum(pca_res$sdev^2)
cum_pve <- cumsum(pve)
threshold <- 0.80
idx80 <- which(cum_pve >= threshold)[1]

cat("Smallest number of PCs to reach >= 80% variance explained =", idx80, "\n")

df <- data.frame(
  PC = seq_along(pve),              # PC index (1..100)
  PVE = pve,                        # Proportion of variance explained by each PC
  CumulativePVE = cum_pve           # Cumulative proportion
)

library(ggplot2)

ggplot(df, aes(x = PC, y = CumulativePVE)) +
  geom_line() +
  geom_point(size = 1) +
  
  # Horizontal line at the threshold (80%)
  geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
  
  # Vertical line at idx80
  geom_vline(xintercept = idx80, color = "blue", linetype = "dashed") +
  
  # Single text annotation near the intersection of x = idx80 and y = threshold
  annotate(
    "text", 
    x = idx80, 
    y = threshold, 
    label = paste0("PC=", idx80, "\n(", round(threshold*100, 0), "%)"),
    color = "blue", 
    vjust = -0.5,     # adjusts vertical position slightly above the crossing
    hjust = 1.2       # adjusts horizontal position slightly left of the crossing
  ) +
  
  labs(
    title = "Scree Plot: Cumulative Variance Explained by PCs",
    x = "Principal Component",
    y = "Cumulative Proportion of Variance Explained"
  ) +
  theme_minimal()

## Obtain PCA embeddings
train_pca_embedding_matrix <- pca_res$x[,1:idx80]  # PCA components (samples x n.comp)
train_pca_embedding <- cbind(train_embedding_transformed[,c("FID","IID")],as.data.table(train_pca_embedding_matrix))
write.table(train_pca_embedding, file=paste0("gpticd3072_PCA",idx80,"_EUR_train_embed_normalize.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

## Run ICA to obtain independent embeddings
num_components <- idx80
ica_res <- fastICA(train_embedding_matrix, n.comp = num_components)

train_ica_embedding_matrix <- ica_res$S  # independent components (samples x n.comp)
train_ica_embedding <- cbind(train_embedding_transformed[,c("FID","IID")],as.data.table(train_ica_embedding_matrix))

write.table(train_ica_embedding, file=paste0("gpticd3072_ICA",idx80,"_EUR_train_embed_normalize.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")