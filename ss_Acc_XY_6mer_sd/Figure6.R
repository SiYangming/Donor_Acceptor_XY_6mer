rm(list = ls())
source("../code/kmers-function.R")
source("../code/plot.R")
library(dplyr)
library(ggpubr)
load(file = "../data/sites_RC.rda")
XYmers <- paste0(kmer(2), 1)
modes <- c("Alternative", "Common", "Constitutive")
poss <- "Acceptor"
regions <- "+2_+45"
x_title <- "Acceptor splicing mode"
y_title <- "Dispersion ρ"
y_column <- "Conservation"
stat.test.file <- "data/Ass_all_region.csv"
plot.file <- "figures/all_DownRegion.pdf"

stat.test.df <- performTtests(data = RC, modes = modes, 
                              position = poss, XYmers = XYmers, 
                              y_value_column = y_column, 
                              outputFileName = stat.test.file)

RC_lst <- RC %>%
  filter(
    ss_type %in% modes, pos == poss,
    XYmer %in% XYmers, region == regions
  )

G1_1 <- c("AA1", "AT1", "CA1", "CG1", "CT1", "GC1", "GT1", "TG1")
G1_2 <- c("TC1")
G2 <- c("AG1", "GA1", "TA1")
G3 <- c("AC1", "CC1", "GG1", "TT1")

lapply(list(G1_1, G1_2, G2, G3), function(x){
  BarChartWithSignificance(data = RC_lst, 
                           XYmers = x, region_filter = regions, 
                           x_axis_title = x_title, 
                           y_axis_title = y_title, 
                           y_value_column = y_column)
}) -> gg

library(cowplot)
gg_all <- plot_grid(gg[[1]],
  plot_grid(gg[[2]], gg[[3]], gg[[4]],
    rel_widths = c(1, 3, 4),
    labels = c("", "B", "C"), nrow = 1
  ),
  ncol = 1, labels = "A"
)
ggsave(gg_all, filename = plot.file, width = 12, height = 6)

