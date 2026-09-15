library(readr)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(tools)

# CLI args: <output PNG> <input tsv.gz> [<input tsv.gz> ...]
args <- commandArgs(trailingOnly = TRUE)
output_png <- args[1]
files <- args[-1]

# Chromosome ordering
chrom_levels <- c(paste0("chr", 1:22))

# Read one file and collapse it to per-(chr, super_bin) mean copy number, so we
# never hold every file's per-bin rows in memory at once.
summarise_file <- function(f) {
  message(paste("Reading", f))
  read_tsv(
    f,
    col_types = cols_only(
      chr   = col_character(),
      start = col_integer(),
      end   = col_integer(),
      copy  = col_double()
    )
  ) %>%
    filter(chr %in% chrom_levels) %>%
    mutate(
      chr = factor(chr, levels = chrom_levels),
      cell_line_name = file_path_sans_ext(file_path_sans_ext(basename(f))),
      copy = pmax(-2, pmin(2, copy))
    ) %>%
    arrange(chr, start) %>%
    group_by(chr) %>%
    mutate(super_bin = (row_number() - 1L) %/% 100L + 1L) %>%
    group_by(cell_line_name, chr, super_bin) %>%
    summarise(copy = if (all(is.na(copy))) {
      NA_real_
    } else {
      mean(copy, na.rm = TRUE)
    }, .groups = "drop")
}

message(paste("Reading", length(files), "files"))
cn_summary <- lapply(files, summarise_file) %>%
  bind_rows() %>%
  mutate(cell_line_name = factor(
    cell_line_name,
    levels = sort(unique(cell_line_name))
  ))

# Heatmap
message("Making plot")
p <-
  ggplot(cn_summary, aes(
    x = cell_line_name,
    y = -super_bin,
    fill = copy
  )) +
  geom_raster() +
  facet_grid(chr ~ .,
             scales = "free_y",
             space = "free_y",
             switch = "y") +
  scale_fill_distiller(
    palette = "RdBu",
    direction = -1,
    limits = c(-2, 2),
    oob = scales::squish,
    na.value = "grey85",
    name = "log2 relative\ncopy number"
  ) +
  labs(x = "Cell line", y = "Chromosome") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7
    ),
    panel.spacing.y = unit(0, "pt")
  )

message(paste("Writing plot to", output_png))
ggsave(
  filename = output_png,
  plot = p,
  width = 4 + length(files) / 2,
  height = 15
)
