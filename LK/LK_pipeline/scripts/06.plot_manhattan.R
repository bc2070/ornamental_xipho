#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript 06.plot_manhattan.R config.ini method output_prefix")
}

suppressPackageStartupMessages(library(ggplot2))

config <- normalizePath(args[1])
method <- args[2]
out_prefix <- args[3]


# ==========================================================
# Read configuration file
# ==========================================================

read_cfg <- function(file) {

  x <- readLines(file, warn = FALSE)

  # Remove comments and empty lines
  x <- x[!grepl("^\\s*[#;]", x) & nzchar(trimws(x))]

  kv <- strsplit(x, "=", fixed = TRUE)

  out <- list()

  for (z in kv) {
    if (length(z) >= 2) {
      out[[trimws(z[1])]] <- trimws(z[2])
    }
  }

  out
}


cfg <- read_cfg(config)


# ==========================================================
# Define pipeline root
# ==========================================================

root <- dirname(config)

base <- file.path(
  root,
  "output",
  "baseline"
)

input <- file.path(
  base,
  "observed_vs_baseline.txt"
)

threshold_file <- file.path(
  base,
  "baseline_threshold.txt"
)


# ==========================================================
# Check input files
# ==========================================================

if (!file.exists(input)) {
  stop(
    paste0(
      "Input file not found: ",
      input
    )
  )
}

if (!file.exists(threshold_file)) {
  stop(
    paste0(
      "Threshold file not found: ",
      threshold_file
    )
  )
}


# ==========================================================
# Read data
# ==========================================================

df <- read.table(
  input,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

threshold_data <- read.table(
  threshold_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

th <- threshold_data$threshold_score[1]

df$chrom <- as.character(df$chrom)


# ==========================================================
# Sort chromosomes naturally
# ==========================================================

natural_key <- function(x) {

  suppressWarnings(v <- as.numeric(x))

  ifelse(
    is.na(v),
    1e9,
    v
  )
}


chroms <- unique(
  df$chrom[
    order(natural_key(df$chrom))
  ]
)

df$chrom <- factor(
  df$chrom,
  levels = chroms
)


# ==========================================================
# Define two-color palette
# ==========================================================

two_colors <- c(
  "#2F5597",
  "#B7C9E2"
)


# ==========================================================
# Manhattan plot
# ==========================================================

p <- ggplot(
  df,
  aes(
    x = pos,
    y = score,
    color = chrom
  )
) +

  geom_point(
    alpha = 0.65,
    size = 0.8
  ) +

  # Alternate the two colors across chromosomes
  scale_color_manual(
    values = rep(
      two_colors,
      length.out = length(chroms)
    )
  ) +

  facet_grid(
    . ~ chrom,
    scales = "free_x",
    space = "free_x"
  ) +

  geom_hline(
    yintercept = th,
    linetype = "dashed",
    linewidth = 0.6
  ) +

  theme_bw() +

  theme(

    # Remove x-axis labels and ticks
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),

    # Remove spacing between chromosomes
    panel.spacing = unit(0, "lines"),

    # Remove legend
    legend.position = "none",

    # Remove background grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    # Chromosome labels
    strip.background = element_rect(
      fill = "grey95",
      color = "white"
    ),

    strip.text = element_text(
      size = 8,
      face = "bold"
    )
  ) +

  labs(
    title = paste0(
      "Manhattan Plot of ",
      method
    ),

    subtitle = paste0(
      "97.5% Threshold: ",
      round(th, 4)
    ),

    x = "Genomic Position",

    y = paste0(
      "-log10(",
      method,
      ")"
    )
  )


# ==========================================================
# Create output directory
# ==========================================================

output_dir <- dirname(out_prefix)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ==========================================================
# Save figures
# ==========================================================

ggsave(
  paste0(out_prefix, ".png"),
  p,
  width = 14,
  height = 6,
  dpi = 300
)

ggsave(
  paste0(out_prefix, ".pdf"),
  p,
  width = 14,
  height = 6
)

cat(
  "Manhattan plot written to:\n",
  paste0(out_prefix, ".png"),
  "\n",
  paste0(out_prefix, ".pdf"),
  "\n"
)

