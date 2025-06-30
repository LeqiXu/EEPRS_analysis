# Packaged prepare
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(gridExtra)
library(cowplot)
library(grid)
library(scales)
library(readr)
library(ggpattern)

my_theme <- theme(
  # shrink axis titles
  axis.title.x = element_text(size = 8, face = "bold"),
  axis.title.y = element_text(size = 8, face = "bold"),
  # shrink tick labels
  axis.text.x  = element_text(size = 5, face = "bold"),
  axis.text.y  = element_text(size = 6, face = "bold"),
  # leave all other text at your preferred size
  plot.title   = element_text(size = 10, face = "bold"),
  strip.text   = element_text(size = 10, face = "bold"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_line(color = "grey", size = 0.5),
  panel.grid.minor.y = element_line(color = "lightgrey", size = 0.25),
  legend.position  = "top",
  legend.title     = element_text(size = 7, face = "bold"),     # title size
  legend.text      = element_text(size = 7),
  plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "cm")
)

setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/plot")

library(data.table)
library(ggplot2)

# Read and format each file
Word2Vec_corr <- fread("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/trait_embedding_corr/UKB_word2vec100_train_trait_corr.txt")[, .(rg, rg_se, emd_i, trait)]
Word2Vec_corr[, emd_i := paste0("E_", emd_i)]
Word2Vec_corr[, embedding_approach := "Word2Vec"]

Word2Vec_PCA_corr <- fread("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/trait_embedding_corr/UKB_word2vec100_pca_train_trait_corr.txt")[, .(rg, rg_se, emd_i, trait)]
Word2Vec_PCA_corr[, emd_i := paste0("PCA_", emd_i)]
Word2Vec_PCA_corr[, embedding_approach := "Word2Vec_PCA"]
word2vec_pca_sig_embedding <- paste0("PCA_",c(1,2,7,9,4,3,5,18,11,8,20,15,6,10,16,17,14,13,21,25,26,12,27,24))
Word2Vec_PCA_corr <- Word2Vec_PCA_corr[which(Word2Vec_PCA_corr$emd_i %in% word2vec_pca_sig_embedding),]

Word2Vec_ICA_corr <- fread("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/trait_embedding_corr/UKB_word2vec100_ica_train_trait_corr.txt")[, .(rg, rg_se, emd_i, trait)]
Word2Vec_ICA_corr[, emd_i := paste0("ICA_", emd_i)]
Word2Vec_ICA_corr[, embedding_approach := "Word2Vec_ICA"]
word2vec_ica_sig_embedding <- paste0("ICA_",c(25,19,12,2,15,11,22,27,4,14,24,5,8,1,10,20,16,13,9,6,29,7))
Word2Vec_ICA_corr <- Word2Vec_ICA_corr[which(Word2Vec_ICA_corr$emd_i %in% word2vec_ica_sig_embedding),]

GPT_PCA_corr <- fread("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/trait_embedding_corr/UKB_gpticd3072_pca_train_trait_corr.txt")[, .(rg, rg_se, emd_i, trait)]
GPT_PCA_corr[, emd_i := paste0("PCA_", emd_i)]
GPT_PCA_corr[, embedding_approach := "GPT_PCA"]
gpt_pca_sig_embedding <- paste0("PCA_",c(2,3,4,12,6,1,9,15,14,5,22,21,38,34,10,18,11,20,16,42,40,32,19,24,23,33,26,37,17,7,44,55,56,8,43,54,45,29,48,39,25,35,65,46,51,31,28,67,13,41,30,60))
GPT_PCA_corr <- GPT_PCA_corr[which(GPT_PCA_corr$emd_i %in% gpt_pca_sig_embedding),]

GPT_ICA_corr <- fread("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/trait_embedding_corr/UKB_gpticd3072_ica_train_trait_corr.txt")[, .(rg, rg_se, emd_i, trait)]
GPT_ICA_corr[, emd_i := paste0("ICA_", emd_i)]
GPT_ICA_corr[, embedding_approach := "GPT_ICA"]
gpt_ica_sig_embedding <- paste0("ICA_",c(54,1,32,15,22,9,56,67,43,25,18,38,12,2,50,20,57,28,37,52,30,49,29,23,21,11,3,51,40,55,42,19,48,5,63,24,34,45,44,35,31,4,7,46))
GPT_ICA_corr <- GPT_ICA_corr[which(GPT_ICA_corr$emd_i %in% gpt_ica_sig_embedding),]

#--- Helper functions ---
add_significance <- function(df) {
  df %>%
    mutate(
      significant = ifelse(abs(rg) > 1.96 * rg_se, "*", ""),
      rg_label = paste0(round(rg, 2), significant)
    )
}

cluster_order <- function(df) {
  mat <- df %>%
    select(trait, emd_i, rg) %>%
    pivot_wider(names_from = emd_i, values_from = rg) %>%
    column_to_rownames("trait") %>%
    as.matrix()

  trait_order <- rownames(mat)[hclust(dist(mat))$order]
  emd_order <- colnames(mat)[hclust(dist(t(mat)))$order]

  df %>%
    mutate(
      trait = factor(trait, levels = trait_order),
      emd_i = factor(emd_i, levels = emd_order)
    )
}

heatmap_plot_clustered <- function(df) {
  ggplot(df, aes(x = emd_i, y = trait, fill = rg)) +
    geom_tile(color = "grey80", na.rm = FALSE) +
    geom_text(aes(label = significant), color = "black", size = 2, na.rm = TRUE) +
    scale_fill_gradient2(
      low = "#4575B4", mid = "white", high = "#D73027",
      midpoint = 0, limits = c(-1, 1), name = "Genetic\nCorrelation", na.value = "grey95"
    ) +
    facet_wrap(~embedding_approach, scales = "free_x", ncol = 1) +
    labs(x = "Embedding Dimension", y = "Traits") +
    theme_classic() +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
      axis.text.y = element_text(size = 5),
      strip.text = element_text(size = 7, face = "bold"),
      legend.position  = "right",
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 5),
      legend.key.size = unit(0.3, "cm") 
    )
}
#--- Figure 1: Word2Vec ---
fig1_data <- Word2Vec_corr %>%
  add_significance() %>%
  cluster_order()

fig1_heatmap_clustered <- heatmap_plot_clustered(fig1_data)
print(fig1_heatmap_clustered)


# Save the figure to match the journal guidelines
ggsave(filename = "Figure2.pdf",
       plot = fig1_heatmap_clustered,
       width = 6.5, 
       height = 5,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication

#--- Figure 2: Word2Vec PCA & ICA ---
fig2_data <- rbind(Word2Vec_PCA_corr, Word2Vec_ICA_corr) %>%
  add_significance() %>%
  mutate(
    embedding_approach = factor(
      embedding_approach, 
      levels = c("Word2Vec_PCA", "Word2Vec_ICA")
    )
  ) %>%
  group_by(embedding_approach) %>%
  group_modify(~ cluster_order(.x)) %>%
  ungroup()


fig2_heatmap_clustered <- heatmap_plot_clustered(fig2_data)
print(fig2_heatmap_clustered)

# Save the figure to match the journal guidelines
ggsave(filename = "FigureS3.pdf",
       plot = fig2_heatmap_clustered,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication

#--- Figure 3: GPT PCA & ICA ---
fig3_data <- rbind(GPT_PCA_corr, GPT_ICA_corr) %>%
  add_significance() %>%
  mutate(
    embedding_approach = factor(
      embedding_approach, 
      levels = c("GPT_PCA", "GPT_ICA")
    )
  ) %>%
  group_by(embedding_approach) %>%
  group_modify(~cluster_order(.x)) %>%
  ungroup()

fig3_heatmap_clustered <- heatmap_plot_clustered(fig3_data)
print(fig3_heatmap_clustered)

# Save the figure to match the journal guidelines
ggsave(filename = "FigureS4.pdf",
       plot = fig3_heatmap_clustered,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication