#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript plot_manhattan.R <input> <threshold_file> <output_prefix>")

library(ggplot2)
library(dplyr)

df <- read.table(args[1], header=FALSE)
colnames(df) <- c("CHR", "Start", "End", "Prob", "logP")

thresh_lines <- readLines(args[2])
target_line <- thresh_lines[grep("Threshold_97.5", thresh_lines)]
sig_y <- as.numeric(gsub(".*: ", "", target_line))

chr_order <- unique(df$CHR[order(as.numeric(gsub("[a-zA-Z]", "", df$CHR)))])
df$CHR <- factor(df$CHR, levels = chr_order)

n_chr <- length(levels(df$CHR))
binary_colors <- rep(c("#34495E", "#3498DB"), length.out = n_chr)

p <- ggplot(df, aes(x = Start, y = logP, color = CHR)) +
  geom_point(alpha = 0.6, size = 0.8) +
  facet_grid(. ~ CHR, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = sig_y, color = "#E74C3C", linetype = "dashed", linewidth = 0.6) +
  scale_color_manual(values = binary_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0, "lines"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "white"),
    strip.text = element_text(size = 8, face = "bold")
  ) +
  labs(
    title = "Manhattan Plot of Cumulative Probabilities",
    subtitle = paste0("97.5% Threshold Line: ", round(sig_y, 4)),
    x = "Genomic Position",
    y = expression(-log[10](Probability))
  )

ggsave(paste0(args[3], ".png"), plot = p, width = 14, height = 6, dpi = 300)
ggsave(paste0(args[3], ".pdf"), plot = p, width = 14, height = 6)
