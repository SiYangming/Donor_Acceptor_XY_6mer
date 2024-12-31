rm(list = ls())
source("../code/function.R")
source("../code/kmers-function.R")
source("../code/plot.R")
library(dplyr)
library(ggpubr)
load(file = "../data/sites_RC.rda")
XYmers <- paste0(kmer(2), 1)
modes <- c("Alternative", "Common", "Constitutive")
poss <- "Acceptor"
regions <- "-6_+1"
x_title <- "Acceptor splicing mode"
y_title <- "Dispersion ρ"
y_column <- "Conservation"
stat.test.file <- "data/Ass_all_region.csv"
plot.file <- "figures/all_UpRegion.pdf"


RC_lst <- RC %>%
  filter(
    ss_type %in% modes, pos == poss,
    XYmer %in% XYmers, region == regions
  )

stat.test.df <- performTtests(data = RC, modes = modes, 
                              position = poss, XYmers = XYmers, 
                              y_value_column = y_column, 
                              outputFileName = stat.test.file)

G1 <- c("CA1", "TC1")
G2 <- c("AC1", "CT1", "TT1", "GA1")
G3 <- setdiff(XYmers, c(G1, G2))
lapply(list(G1, G2, G3), function(x){
  BarChartWithSignificance(data = RC_lst, 
                           XYmers = x, region_filter = regions, 
                           x_axis_title = x_title, 
                           y_axis_title = y_title, 
                           y_value_column = y_column)
}) -> gg

library(cowplot)
gg_all <- plot_grid(
  plot_grid(gg[[1]], gg[[2]],
            rel_widths = c(2, 4),
            labels = c("A", "B"), nrow = 1
                ), 
  gg[[3]], ncol = 1, labels = c("", "C")
)

library(patchwork)
ggsave(gg_all, filename = "figures/all_CenterRegion.pdf", 
       width = 12, height = 6)
