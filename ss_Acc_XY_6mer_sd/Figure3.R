rm(list = ls())
source("../code/kmers-function.R")
library(tidyverse)
load(file = "../data/sites_RC.rda")
XYmers <- paste0(kmer(2), 1)
modes <- c("Alternative", "Common", "Constitutive")
poss <- "Acceptor"
x_title <- "Distance from splicing site"
y_title <- "Dispersion ρ"
y_column <- "Conservation"
plot.file <- "figures/sd_alt_con_norm_XY1.pdf"
table(RC$ss_type)
# All  Alternative       Common Constitutive       Exonic
# 3360         3360         3360         3360         3360
# Intronic       Normal
# 3360         3360

RC %>%
  arrange(type, XYmer) %>%
  filter(pos == poss, ss_type %in% modes, XYmer %in% XYmers) %>%
  select(Separation) %>%
  range(na.rm = TRUE)
# [1] 0.0545 2.6704
disp_y_break <- seq(from = 0, to = 2.8, by = 0.5)
RC %>%
  arrange(type, XYmer) %>%
  filter(pos == poss, ss_type %in% modes, XYmer %in% XYmers) %>%
  select(Conservation) %>%
  range(na.rm = TRUE)
# [1] 0.0064 1.9743
sd_y_break <- seq(from = 0, to = 2.0, by = 0.5)

RC %>%
  arrange(type, XYmer) %>%
  filter(pos == poss, ss_type %in% modes, XYmer %in% XYmers) %>%
  xy_mer_disp_sd_plot_50(
    split_by = "samples",
    disp_y_breaks = disp_y_break,
    sd_y_breaks = sd_y_break,
    is_facet = TRUE
  ) -> gg

ggsave(
  gg$sd +
    ggtitle("") +
    ylab(y_title) +
    xlab(x_title),
  filename = plot.file,
  height = 12, width = 12
)
