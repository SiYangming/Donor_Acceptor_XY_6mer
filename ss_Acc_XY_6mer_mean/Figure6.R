rm(list = ls())
source("../code/kmers-function.R")

### Alternative Acceptor
load(file = "data/RC.rda")
XYmers <- paste0(kmer(2), 1)
table(RC$ss_type)
# All  Alternative       Common Constitutive       Exonic
# 3360         3360         3360         3360         3360
# Intronic       Normal
# 3360         3360

RC$ss_type <- tolower(RC$ss_type)
modes <- c("normal", "exonic", "intronic")
RC %>%
  arrange(type, XYmer) %>%
  filter(
    (ss_type %in% modes),
    XYmer %in% XYmers
  ) %>%
  select(Separation) %>%
  range(na.rm = TRUE)
# [1] 0.1034 2.4003
RC %>%
  arrange(type, XYmer) %>%
  filter(
    (ss_type %in% modes),
    XYmer %in% XYmers
  ) %>%
  select(Conservation) %>%
  range(na.rm = TRUE)
# [1] 0.0196 1.8804


RC %>%
  arrange(type, XYmer) %>%
  filter(
    (ss_type %in% modes),
    XYmer %in% XYmers
  ) %>%
  xy_mer_disp_sd_plot_50(
    split_by = "samples",
    disp_y_breaks = seq(from = 0, to = 2.5, by = 0.5),
    sd_y_breaks = seq(from = 0, to = 2, by = 0.2),
    is_facet = TRUE
  ) -> gg
ggsave(
  gg$disp +
    ggtitle("") +
    ylab("Separability δ") +
    xlab("Distance from splicing site"),
  filename = "figures/alt_acc_subset_disp_XY1.pdf",
  height = 12, width = 12
)

ggsave(
  gg$sd +
    ggtitle("") +
    ylab("Conservation  δ") +
    xlab("Distance from splicing site"),
  filename = "figures/alt_subset_sd_XY1.pdf",
  height = 12, width = 12
)
