library(readr)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
chrom_sizes_path <- args[1]
bins_path <- args[2]
segments_path <- args[3]
output_png <- args[4]

# read hg38 chromosome sizes
chrom_sizes <- read_tsv(
  chrom_sizes_path,
  col_names = c("CONTIG", "length"),
  show_col_types = FALSE
) %>%
  mutate(
    offset = lag(cumsum(length), default = 0),
    midpoint = offset + length / 2
  )

# read bin log2 CRs
bins <- read_tsv(
  bins_path,
  show_col_types = FALSE
) %>%
  filter(
    valid,
    !is.na(copy)
  ) %>%
  inner_join(
    select(chrom_sizes, CONTIG, offset),
    by = c("chr" = "CONTIG")
  ) %>%
  mutate(
    genome_pos = (start + end) / 2 + offset
  )

bins_plot <- bins %>%
  group_by(chr) %>%
  slice_sample(
    n = 10000
  ) %>%
  ungroup()

# read segment log2 CRs
segments <- read_tsv(
  segments_path,
  show_col_types = FALSE
) %>%
  inner_join(
    select(chrom_sizes, CONTIG, offset),
    by = "CONTIG"
  ) %>%
  mutate(
    genome_start = START + offset,
    genome_end = END + offset
  )

p <- ggplot() +
  geom_point(
    data = bins_plot,
    aes(
      x = genome_pos,
      y = copy,
      color = copy
    ),
    size = 0.15,
    alpha = 0.5
  ) +
  geom_segment(
    data = segments,
    aes(
      x = genome_start,
      xend = genome_end,
      y = LOG2_COPY_RATIO_POSTERIOR_50,
      yend = LOG2_COPY_RATIO_POSTERIOR_50
    ),
    linewidth = 0.8,
    color = "black",
    lineend = "round"
  ) +
  geom_vline(
    data = chrom_sizes,
    aes(xintercept = offset),
    color = "grey80",
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = 0,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_gradient2(
    low = "blue",
    mid = "grey60",
    high = "red",
    midpoint = 0,
    limits = c(-3, 3),
    oob = scales::squish,
    name = expression(log[2]("copy ratio"))
  ) +
  scale_x_continuous(
    breaks = chrom_sizes$midpoint,
    labels = sub("^chr", "", chrom_sizes$CONTIG),
    expand = c(0, 0)
  ) +
  coord_cartesian(
    ylim = c(-3, 3)
  ) +
  labs(
    x = NULL,
    y = expression(log[2]("copy ratio"))
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 8),
    panel.border = element_rect(linewidth = 0.5)
  )

ggsave(
  output_png,
  p,
  width = 14,
  height = 4,
  dpi = 300
)
