rm(list = ls())
source("../code/kmers-function.R")
### Alternative Acceptor
load(file = "../data/sites_RC.rda")
XYmers <- paste0(kmer(2), 1)
modes <- c("Normal", "Exonic", "Intronic")
poss <- "Acceptor"
x_title <- "Distance from splicing site"
y_title <- "Dispersion ρ"
y_column <- "Conservation"
plot.file <- "figures/alt_subset_sd_XY1.pdf"

table(RC$ss_type)
# All  Alternative       Common Constitutive       Exonic 
# 6400         6400         6400         6400         6400 
# Intronic       Normal 
# 6400         6400 

RC %>%
  arrange(type, XYmer) %>%
  filter(pos == poss,
         (ss_type %in% modes),
         XYmer %in% XYmers
  ) %>%
  select(Separation) %>%
  range(na.rm = TRUE)
# [1] 0.0746 2.5494
disp_y_break <- seq(from = 0, to = 2.6, by = 0.5)

RC %>%
  arrange(type, XYmer) %>%
  filter(pos == poss,
    (ss_type %in% modes),
    XYmer %in% XYmers
  ) %>%
  select(Conservation) %>%
  range(na.rm = TRUE)
# [1] 0.0196 1.9475
sd_y_break <- seq(from = 0, to = 2, by = 0.5)


RC %>%
  arrange(type, XYmer) %>%
  filter(pos == poss,
    (ss_type %in% modes),
    XYmer %in% XYmers
  ) %>%
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
