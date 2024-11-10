rm(list = ls())
source("../code/kmers-function.R")
library(tidyverse)
load(file = "data/RC.rda")
RC%>%filter(pos=="Donor")%>%View()
RC <- RC %>%
  distinct()
save(RC, file = "data/RC.rda")
table(RC$ss_type)
# All  Alternative       Common Constitutive       Exonic
# 3200         3200         3200         3200         3200
# Intronic       Normal
# 3200         3200

XYmers <- paste0(kmer(2), 1)
RC$ss_type <- tolower(RC$ss_type)
modes <- c("alternative", "common", "constitutive")

RC %>%
  arrange(type, XYmer) %>%
  filter(ss_type %in% modes, XYmer %in% XYmers) %>%
  select(Separation) %>%
  range(na.rm = TRUE)
# [1] 0.0545 2.6805
RC %>%
  arrange(type, XYmer) %>%
  filter(ss_type %in% modes, XYmer %in% XYmers) %>%
  select(Conservation) %>%
  range(na.rm = TRUE)
# [1] 0.0064 1.9752


RC %>%
  arrange(type, XYmer) %>%
  filter(ss_type %in% modes, XYmer %in% XYmers) %>%
  xy_mer_disp_sd_plot_50(
    split_by = "samples",
    disp_y_breaks = seq(from = 0, to = 3.0, by = 0.5),
    sd_y_breaks = seq(from = 0, to = 2, by = 0.2),
    is_facet = TRUE
  ) -> gg
ggsave(
  gg$disp +
    ggtitle("") +
    ylab("Separability δ") +
    xlab("Distance from splicing site"),
  filename = "figures/acc_disp_alt_con_norm_XY1.pdf",
  height = 12, width = 12
)

ggsave(
  gg$sd +
    ggtitle("") +
    ylab("Conservation") +
    xlab("Distance from splicing site"),
  filename = "figures/sd_alt_con_norm_XY1.pdf",
  height = 12, width = 12
)
