## auto and tuning method comparison
library(data.table)
library(dplyr)
library(tidyr)
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
  legend.title     = element_text(size = 9, face = "bold"),     # title size
  legend.text      = element_text(size = 7),
  plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "cm")
)

setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/plot")

# Read and format each file
Word2Vec_h2 <- fread("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/UKB_word2vec100_train_h2.txt")[, .(h2, h2_se, embedding)]
Word2Vec_h2[, embedding := paste0("E_", embedding)]
Word2Vec_h2[, embedding_approach := "Word2Vec"]

Word2Vec_PCA_h2 <- fread("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/UKB_word2vec100_pca_train_h2.txt")[, .(h2, h2_se, embedding)]
Word2Vec_PCA_h2[, embedding := paste0("PCA_", embedding)]
Word2Vec_PCA_h2[, embedding_approach := "Word2Vec_PCA"]

Word2Vec_ICA_h2 <- fread("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/UKB_word2vec100_ica_train_h2.txt")[, .(h2, h2_se, embedding)]
Word2Vec_ICA_h2[, embedding := paste0("ICA_", embedding)]
Word2Vec_ICA_h2[, embedding_approach := "Word2Vec_ICA"]

GPT_PCA_h2 <- fread("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/UKB_gpticd3072_pca_train_h2.txt")[, .(h2, h2_se, embedding)]
GPT_PCA_h2[, embedding := paste0("PCA_", embedding)]
GPT_PCA_h2[, embedding_approach := "GPT_PCA"]

GPT_ICA_h2 <- fread("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/h2/embedding_h2/UKB_gpticd3072_ica_train_h2.txt")[, .(h2, h2_se, embedding)]
GPT_ICA_h2[, embedding := paste0("ICA_", embedding)]
GPT_ICA_h2[, embedding_approach := "GPT_ICA"]

# Combine all data
combined_h2 <- rbindlist(list(Word2Vec_h2, Word2Vec_PCA_h2, Word2Vec_ICA_h2, GPT_PCA_h2, GPT_ICA_h2))
combined_h2$embedding_approach <- factor(combined_h2$embedding_approach, levels = c("Word2Vec","Word2Vec_PCA","Word2Vec_ICA","GPT_PCA","GPT_ICA"))

# Compute significance (e.g., h2 significantly > 0 at 95% CI)
combined_h2[, significant := h2 - 1.96 * h2_se > 0]
combined_h2[, significant := factor(significant, levels = c("TRUE", "FALSE"))]

# Plotting
# Define consistent colors for each method
method_colors <- c(
  "Word2Vec"     = "#8DA0CB",
  "Word2Vec_PCA" = "#E78AC3",
  "Word2Vec_ICA" = "#FC8D62",
  "GPT_PCA"      = "#66C2A5",
  "GPT_ICA"      = "#A6D854"
)

# Update your plot with fixed colors
combined_plot <- ggplot(combined_h2, aes(x = embedding_approach, y = h2, color = embedding_approach)) +
  geom_violin(trim = TRUE, position = position_dodge(width = 0.8), width = 0.8, adjust = 1.5, scale = "width", show.legend = FALSE) +
  geom_jitter(aes(shape = significant), width = 0.15, size = 1, alpha = 1, show.legend = TRUE) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, position = position_dodge(width = 0.4), 
               color = alpha("black", alpha = 0.6), size = 0.2, show.legend = FALSE) +
  scale_color_manual(values = method_colors, guide = "none") +  # remove color legend
  scale_shape_manual(values = c("TRUE" = 19,"FALSE" = 1), name = "Significant") +
  labs(x = "Embedding Approach",
       y = paste0("Heritability")) +
  theme_classic(base_size = 14) +
  my_theme

print(combined_plot)

# Save the figure to match the journal guidelines
ggsave(filename = "FigureS2.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication