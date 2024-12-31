rm(list = ls())
library(dplyr)
library(ggpubr)
source("../code/kmers-function.R")
source("../code/plot.R")
load(file = "../data/sites_RC.rda")
XYmers <- paste0(kmer(2), 1)
modes<-c("Normal", "Exonic", "Intronic")
poss <- "Acceptor"
regions <- "+2_+45"
x_title <- "Acceptor splicing mode"
y_title <- "ispersion  ρ"
y_column <- "Conservation"
stat.test.file <- "data/AltAss_all_region.csv"
plot.file <- "figures/alt_subset_DownRegion.pdf"

RC_lst <- RC %>%
  filter(
    ss_type %in% modes, pos == poss,
    XYmer %in% XYmers, region == regions
  )
stat.test.df <- performTtests(data = RC, modes = modes, 
                              position = poss, XYmers = XYmers, 
                              y_value_column = y_column, 
                              outputFileName = stat.test.file)

G1 <- c("AA1", "AT1", "CT1", "TA1", "TG1")
G2 <- c("AG1", "CA1", "CC1", "CG1", "GA1", "GT1", "TC1")
G3 <- c("AC1", "GC1", "GG1", "TT1")

lapply(list(G1, G2, G3), function(x){
  BarChartWithSignificance(data = RC_lst, 
                           XYmers = x, region_filter = regions, 
                           x_axis_title = x_title, 
                           y_axis_title = y_title, 
                           y_value_column = y_column)
}) -> gg

library(cowplot)
gg_all <- plot_grid(
  plot_grid(gg[[1]], gg[[3]],
    rel_widths = c(5, 4),
    labels = c("A", "C"), nrow = 1
  ), gg[[2]],
  ncol = 1, labels = c("", "B")
)
ggsave(gg_all,
  filename = plot.file,
  width = 12, height = 6
)
