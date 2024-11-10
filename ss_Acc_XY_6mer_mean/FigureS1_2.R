rm(list = ls())
library(tidyverse)

###### 1,2,3,6-mer ######
freq_files <- list.files("data",
  pattern = "*_freqs.csv.gz",
  full.names = TRUE
)
freq <- lapply(freq_files, function(freq_file) {
  rc <- read_csv(freq_file, col_types = as.list(c("c", rep("d", 95))), )
  rc[is.na(rc)] <- 0
  outfile <- gsub(".gz", "", freq_file)
  outfile <- gsub("1mer", "Xmer", outfile)
  outfile <- gsub("2mer", "XYmer", outfile)
  outfile <- gsub("3mer", "XYZmer", outfile)
  colnames(rc) <- "kmer"
  write_csv(rc, file = outfile)

  return(rc)
})

##### relative frequence #####

if (!file.exists("data/freq.rda")) {
  freq_files <- list.files("data",
    pattern = "[XYZ]{1,3}mer_freqs.csv",
    full.names = TRUE
  )
  freq <- lapply(freq_files, function(freq_file) {
    rc <- read_csv(freq_file)
    colnames(rc) <- c("kmer", 1:(ncol(rc) - 1))
    rc %>%
      reshape2::melt(variable.name = "site", value.name = "count") %>%
      dplyr::mutate(
        type = gsub("data/|_[XYZ]{1,3}mer_freqs.csv", "", freq_file),
        xyzmer = str_split(freq_file, "_")[[c(1, 3)]],
        total_count = unique(colSums(rc[-1])),
        freq = count / total_count,
        bg_freq = 1 / nrow(rc), rf_freq = freq / bg_freq
      )
  })
  freq_df <- do.call(rbind, freq)
  save(freq_df, file = "data/freq.rda")
} else {
  load(file = "data/freq.rda")
}
plotlst <- freq_df %>%
  filter(type == "All_Acceptor") %>%
  group_by(xyzmer) %>%
  group_split()

library(ggplot2)
library(ggbreak)
x_breaks <- (c(seq(1, 49, 10), 50:51, seq(60, 100, 10)) - 50)
x_labels <- if_else(x_breaks <= 0, x_breaks - 1, x_breaks)
x_labels <- if_else(x_labels > 0,
  paste0("+", as.character(x_labels)),
  as.character(x_labels)
)
x_labels <- if_else(x_labels %in% c("-1", "+1"), "", x_labels)
y_breaks <- seq(0, 16, 1)

source("../code/kmers-function.R")
source("../code/layers.R")
plot_lst[[2]] %>%
  group_by(kmer) %>%
  group_split(.keep = TRUE) -> plotlsts
lapply(plotlsts, function(plotlst) {
  plotlst %>%
    dplyr::mutate(
      site0 = as.numeric(site) - 50,
      site1 = if_else(site0 <= 0, site0 - 1, site0),
      XY = xyzmer
    ) %>%
    ggplot(mapping = aes(
      x = site0, y = rf_freq # , color = kmer
    ), color = "black") +
    geom_line(linewidth = 0.2, show.legend = FALSE) +
    geom_point(size = 0.1, show.legend = FALSE) +
    geom_vhline(
      v = c(0, 1), h = c(1),
      vlinetype = c("dashed", "dotted"),
      hlinetype = c("dashed"), vlinewidth = rep(0.2, 2),
      hlinewidth = c(0.2), color = rep("black", 2)
    ) +
    scale_fill_brewer(palette = "Set1") +
    annotate_sites(
      geom = rep("curve", 2),
      x = c(-13, 14), y = c(16, 16),
      xend = c(0, 1), yend = c(16.2, 16.2),
      curvature = c(-0, 0), colour = rep("orange", 2),
      units = list(c(1, "mm"), c(1, "mm")),
      label = c("-1", "+1")
    ) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(breaks = y_breaks) +
    facet_rep_wrap(~kmer,
      ncol = 4, scales = "fixed",
      repeat.tick.labels = TRUE
    ) +
    scale_y_break(c(7.5, 15.5), space = 0.1, scales = 0.1) +
    ggtitle("") +
    ylab("") +
    xlab("") +
    # ylab("Relative frequency") +
    # xlab("Distance from splicing site") +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        size = 10,
        family = "Times",
        face = "bold",
        vjust = 0.4, hjust = 1
      ),
      axis.text.y = element_text(
        family = "Times", size = 10,
        face = "bold"
      ),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right = element_blank()
    )
}) -> gg_lst
for (i in 1:length(gg_lst)) {
  ggsave(gg_lst[[i]],
    filename = paste0("figures/XYmer", i, ".pdf"),
    height = 4, width = 4
  )
}
##### 16 in one plot #####
plot_lst[[2]] %>%
  dplyr::mutate(
    site0 = as.numeric(site) - 50,
    site1 = if_else(site0 <= 0, site0 - 1, site0),
    XY = xyzmer
  ) %>%
  ggplot(mapping = aes(x = site0, y = rf_freq, color = kmer)) +
  geom_line(
    linewidth = 0.2, # show.legend = FALSE
  ) +
  geom_point(aes(x = site0, y = rf_freq),
    size = 0.1
  ) +
  geom_vhline(
    v = c(0, 1), h = c(1),
    vlinetype = c("dashed", "dotted"),
    hlinetype = c("dashed"), vlinewidth = rep(0.2, 2),
    hlinewidth = c(0.2), color = rep("black", 2)
  ) +
  scale_fill_brewer(palette = "Set1") +
  annotate_sites(
    geom = rep("curve", 2),
    x = c(-13, 14), y = c(4, 4),
    xend = c(0, 1), yend = c(4.2, 4.2),
    curvature = c(-0, 0), colour = rep("orange", 2),
    units = list(c(1, "mm"), c(1, "mm")),
    label = c("-1", "+1")
  ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(breaks = y_breaks) +
  scale_y_break(c(7.5, 15.5), space = 0.1, scales = 0.1) +
  ggtitle("") +
  ylab("Relative frequency") +
  xlab("Distance from splicing site") +
  theme_classic2() +
  # theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      size = 10,
      # size=c(rep(10,49),rep(6,2),rep(10,49)),
      family = "Times",
      face = "bold",
      vjust = 0.4, hjust = 1
    ),
    axis.text.y = element_text(
      family = "Times", size = 10,
      face = "bold"
    ),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank()
  ) -> gg1
ggsave(gg1,
  filename = "figures/XYmer_merged.pdf",
  height = 6, width = 6
)

##### 1mer plot #####
y_breaks <- seq(0, 4, 1)
plot_lst[[3]] %>%
  dplyr::mutate(
    site0 = as.numeric(site) - 50,
    site1 = if_else(site0 <= 0, site0 - 1, site0),
    XY = xyzmer
  ) %>%
  ggplot(mapping = aes(
    x = site0, y = rf_freq # , color = kmer
  ), color = "black") +
  geom_line(linewidth = 0.2) +
  geom_point(aes(x = site0, y = rf_freq),
    size = 0.1, show.legend = FALSE
  ) +
  geom_vhline(
    v = c(0, 1), h = c(1),
    vlinetype = c("dashed", "dotted"),
    hlinetype = c("dashed"), vlinewidth = rep(0.2, 2),
    hlinewidth = c(0.2), color = rep("black", 2)
  ) +
  scale_fill_brewer(palette = "Set1") +
  annotate_sites(
    geom = rep("curve", 2),
    x = c(-13, 14), y = c(4, 4),
    xend = c(0, 1), yend = c(4.2, 4.2),
    curvature = c(-0, 0), colour = rep("orange", 2),
    units = list(c(1, "mm"), c(1, "mm")),
    label = c("-1", "+1")
  ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(breaks = y_breaks) +
  facet_rep_wrap(~kmer,
    ncol = 4, scales = "fixed",
    repeat.tick.labels = TRUE
  ) +
  ggtitle("") +
  ylab("Relative frequency") +
  xlab("Distance from splicing site") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      size = 10,
      # size=c(rep(10,49),rep(6,2),rep(10,49)),
      family = "Times",
      face = "bold",
      vjust = 0.4, hjust = 1
    ),
    axis.text.y = element_text(
      family = "Times", size = 10,
      face = "bold"
    )
  ) -> gg1
ggsave(gg1,
  filename = "figures/Xmer.pdf",
  height = 3, width = 12
)

##### 4 in one plot #####
plot_lst[[3]] %>%
  dplyr::mutate(
    site0 = as.numeric(site) - 50,
    site1 = if_else(site0 <= 0, site0 - 1, site0),
    XY = xyzmer
  ) %>%
  ggplot(mapping = aes(x = site0, y = rf_freq, color = kmer)) +
  geom_line(
    linewidth = 0.2, # show.legend = FALSE
  ) +
  geom_point(aes(x = site0, y = rf_freq),
    size = 0.1
  ) +
  geom_vhline(
    v = c(0, 1), h = c(1),
    vlinetype = c("dashed", "dotted"),
    hlinetype = c("dashed"), vlinewidth = rep(0.2, 2),
    hlinewidth = c(0.2), color = rep("black", 2)
  ) +
  scale_fill_brewer(palette = "Set1") +
  annotate_sites(
    geom = rep("curve", 2),
    x = c(-13, 14), y = c(4, 4),
    xend = c(0, 1), yend = c(4.2, 4.2),
    curvature = c(-0, 0), colour = rep("orange", 2),
    units = list(c(1, "mm"), c(1, "mm")),
    label = c("-1", "+1")
  ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(breaks = y_breaks) +
  ggtitle("") +
  ylab("Relative frequency") +
  xlab("Distance from splicing site") +
  theme_classic2() +
  # theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      size = 10,
      # size=c(rep(10,49),rep(6,2),rep(10,49)),
      family = "Times",
      face = "bold",
      vjust = 0.4, hjust = 1
    ),
    axis.text.y = element_text(
      family = "Times", size = 10,
      face = "bold"
    )
  ) -> gg1
ggsave(gg1,
  filename = "figures/Xmer_merged.pdf",
  height = 4, width = 4
)
