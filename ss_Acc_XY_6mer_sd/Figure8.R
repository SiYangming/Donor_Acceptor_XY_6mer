rm(list = ls())
library(dplyr)
library(ggpubr)
source("../code/kmers-function.R")
source("../code/plot.R")
load(file = "../data/sites_RC.rda")
XYmers <- paste0(kmer(2), 1)
modes<-c("Normal", "Exonic", "Intronic")
poss <- "Acceptor"
regions <- "-50_-7"
x_title <- "Acceptor splicing mode"
y_title <- "Dispersion ρ"
y_column <- "Conservation"
stat.test.file <- "data/AltAss_all_region.csv"
plot.file <- "figures/alt_subset_UpRegion.pdf"


RC_lst <- RC %>%
  filter(
    ss_type %in% modes, pos == poss,
    XYmer %in% XYmers, region == regions
  )

### T-test
stat.test.df <- performTtests(data = RC, modes = modes, 
                              position = poss, XYmers = XYmers, 
                              y_value_column = y_column, 
                              outputFileName = stat.test.file)
table(stat.test.df$region, stat.test.df$method)
# T-test
# -50_-7     48
# -6_+1      48
# +2_+45     48

g1_1 <- c("AC1", "CA1", "CC1", "CG1", "AG1", "GC1", "GG1", "GT1")
g1_2 <- c("CT1", "TC1", "TT1")
g2 <- c("TA1", "GA1", "TG1")
g3 <- c("AA1", "AT1")

lapply(list(g1_1, g1_2, g2, g3), function(x){
  BarChartWithSignificance(data = RC_lst, 
                           XYmers = x, region_filter = regions, 
                           x_axis_title = x_title, 
                           y_axis_title = y_title, 
                           y_value_column = y_column)
}) -> gg

library(cowplot)
gg_all <- plot_grid(gg[[1]],
  plot_grid(gg[[2]], gg[[3]], gg[[4]],
    rel_widths = c(3, 3, 2),
    labels = c("", "B", "C"), nrow = 1
  ),
  ncol = 1, labels = "A"
)
ggsave(gg_all,
  filename = plot.file,
  width = 12, height = 6
)
