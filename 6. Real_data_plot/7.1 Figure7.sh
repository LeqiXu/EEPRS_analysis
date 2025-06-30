## 1. Prepare
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(forcats)

my_theme <- theme(
  # shrink axis titles
  axis.title.x = element_text(size = 9, face = "bold"),
  axis.title.y = element_text(size = 9, face = "bold"),
  # shrink tick labels
  axis.text.x  = element_text(size = 6, face = "bold"),
  axis.text.y  = element_text(size = 6, face = "bold"),
  # leave all other text at your preferred size
  plot.title   = element_text(size = 10, face = "bold"),
  strip.text   = element_text(size = 10, face = "bold"),
  panel.grid.major.x = element_line(color = "grey", size = 0.5),
  panel.grid.minor.x = element_line(color = "lightgrey", size = 0.25),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  legend.position  = "top",
  legend.title     = element_text(size = 9, face = "bold"),     # title size
  legend.text      = element_text(size = 7),
  plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "cm")
)

setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/plot")

continuous_trait = c("HDL","LDL","TC","logTG","BMI","WHR","AgeSmk","CigDay","DrnkWk","Glu2h","HbA1c","eGFR","logCp")

## 2. Data
final_cv_results = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/evaluation/MTAG_EEPRS_optimal_R2_AUC.csv"))
comparison_continuous_table <- final_cv_results %>%
  filter(trait %in% continuous_trait) %>%
  pivot_wider(
    names_from = method, 
    values_from = mean_performance
  ) %>%
  mutate(
    difference      = MTAG_EEPRS - PRScs,
    relative_change = (MTAG_EEPRS - PRScs) / PRScs
  ) %>%
  arrange(desc(relative_change))

comparison_binary_table <- final_cv_results %>%
  filter(!(trait %in% continuous_trait)) %>%
  pivot_wider(
    names_from = method, 
    values_from = mean_performance
  ) %>%
  mutate(
    difference      = MTAG_EEPRS - PRScs,
    relative_change = (MTAG_EEPRS - PRScs) / (PRScs - 0.5)
  ) %>%
  arrange(desc(relative_change))
comparison_binary_table = comparison_binary_table[which(comparison_binary_table$PRScs > 0.5),]
comparison_binary_table = comparison_binary_table[which(comparison_binary_table$MTAG_EEPRS > 0.5),]

comparison_table = rbind(comparison_continuous_table,comparison_binary_table)
comparison_table <- comparison_table %>%
  mutate(
    trait_type = if_else(trait %in% continuous_trait, "Continuous Traits", "Binary Traits")
  )
comparison_table$trait_type <- factor(comparison_table$trait_type, levels = c("Continuous Traits", "Binary Traits"))

print(comparison_table)

df_avg <- comparison_table %>%
  group_by(trait_type) %>%
  summarize(avg_improve = mean(relative_change, na.rm = TRUE)) %>%
  ungroup()

## 3. Plot
combined_plot <- ggplot(comparison_table,
       aes(x = reorder(trait, relative_change),
           y = 100 * relative_change,
           fill = relative_change > 0)) +
  geom_col() +
  scale_fill_manual(values = c("grey","#FB9A99")) +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    x = "Trait",
    y = "Relative Improvement (%)"
  ) +
  theme_classic() +
  my_theme +
  guides(fill = FALSE) +
  facet_wrap(~ trait_type, ncol = 1, scales = "free") +
  # Draw one dashed line per facet for average improvement:
  geom_hline(
    data = df_avg,
    aes(yintercept = 100 * avg_improve),
    color = "purple",
    linetype = "dashed",
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = df_avg,
    aes(
      x = 1,
      y = 100 * avg_improve,
      label = paste0("Avg: ", sprintf("%.2f%%", 100 * avg_improve))
    ),
    color = "purple",
    fontface = "bold",
    vjust = -0.5, hjust = -0.1,
    inherit.aes = FALSE,
    size = 3
  )

print(combined_plot)

# Save the figure to match the journal guidelines
ggsave(filename = "Figure7.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication