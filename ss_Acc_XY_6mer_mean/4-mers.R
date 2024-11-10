rm(list = ls())
source("../code/kmers-function.R")
library(tidyverse)
paths <- "data"
freqs_filenames <- list.files(paths,
  pattern = "4mer_freqs.csv.gz",
  full.names = TRUE
)

mer <- "4mer"
for (filename in freqs_filenames) {
  rt <- read_csv(filename,
    col_types = as.list(c("c", rep("d", 95))),
    col_names = T
  )
  colnames(rt)[1] <- mer
  mer_classify <- lapply(kmer(2), function(x) {
    kmer_classification(rt[[1]], kmer = x, group = 2)
  })

  mer_classify_df <- rt %>%
    add_column(do.call(cbind, mer_classify), .after = "4mer")
  outfile <- gsub("_freqs.csv.gz", "_classification.csv.gz", filename)
  write_csv(mer_classify_df, file = outfile, col_names = T)
}

clf_filenames <- list.files(paths,
  pattern = "4mer_classification.csv.gz",
  full.names = TRUE
)
### Run
pattern <- "data\\/|_4mer_classification.csv.gz"
plot_filename <- "data/4-mer_plotdata.rda"
if (!file.exists(plot_filename)) {
  clf_df_lst <- lapply(clf_filenames, function(x) {
    type <- gsub(pattern, "", x)
    rt <- read_csv(x,
      col_types = as.list(c(rep("c", 17), rep("d", 95))),
      col_names = T
    )
    GG_kmer_RidgesData(rt, prefix = "4mer") %>%
      mutate(type = type)
  })

  plotdata <- do.call(rbind, clf_df_lst)
  save(plotdata, file = plot_filename)
} else {
  load(file = plot_filename)
}

RC_filename <- "data/4-mer_RC_matrix.rda"

if (!file.exists(RC_filename)) {
  RC_matrix <- plotdata %>%
    dplyr::group_by(type, site, value, .drop = FALSE) %>%
    # group_split()
    dplyr::summarise(
      mean = mean(frequence, na.rm = TRUE),
      sd = sd(frequence, na.rm = TRUE)
    )
  save(RC_matrix, file = RC_filename)
} else {
  load(RC_filename)
}

##### Calc Separation and Conservation
lst <- RC_matrix %>%
  group_by(type, site, .drop = FALSE) %>%
  group_map(~ calc_Ratio(.x, digits = 4), .keep = TRUE)

RC <- do.call(rbind, lst) %>%
  select(-mean, -sd)
RC[is.na(RC)] <- 0
colnames(RC) <- c("type", "site", "XYmer", "Separation", "Conservation")

range(RC$Separation, na.rm = TRUE)
# [1] 0.0000 3.3125
range(RC$Conservation, na.rm = TRUE)
# [1] 0.0000 2.0282

RC_filename <- "data/4-mer_RC.rda"
k=4
if (!file.exists(RC_filename)) {
  modes <- do.call(rbind, str_split(RC$type, "_"))
  RC <- RC %>%
    mutate(pos = unique(modes[, 2]), ss_type = unique(modes[, 1]))

  RC <- RC %>%
    add_row(
      type = rep(rep(unique(RC$type), each = 32), 3),
      site = rep(98:100, each = 32),
      XYmer = rep(c(paste0(kmer(2), 0), paste0(kmer(2), 1)), 3),
      Separation = NA, Conservation = NA,
      pos = "Acceptor",
      ss_type = rep(rep(unique(RC$ss_type),
        each = 32
      ), 3)
    )
  save(RC, file = RC_filename)
} else {
  load(file = RC_filename)
}

##### Figure 2:  Overall Acceptor Splicing Signal 背景值XY1&XY0 #####
table(RC$ss_type)
# All  Alternative       Common Constitutive       Exonic
# 3200         3200         3200         3200         3200
# Intronic       Normal
# 3200         3200
table(RC$pos)
# Acceptor
# 22400

### Overall Alternative Acceptor 背景值XY1
source("../code/layers.R")
palette <- rep(c("#F8786F", "#609DFF"), 16)
x_breaks <- (c(seq(1, 49, 10), 45:46, 50:51, seq(60, 100, 10)) - 50)
x_labels <- if_else(x_breaks <= 0, x_breaks - 1, x_breaks)
x_labels <- if_else(x_labels > 0,
  paste0("+", as.character(x_labels)),
  as.character(x_labels)
)
x_labels <- if_else(x_labels %in% c("-6", "-5", "-1", "+1"), "", x_labels)
y_breaks <- seq(0, 3.5, 0.5)
RC %>%
  arrange(type, XYmer) %>%
  dplyr::filter(ss_type == "All") %>%
  group_by(XYmer, site) %>%
  dplyr::mutate(
    site0 = site - 50,
    site1 = if_else(site0 <= 0, site0 - 1, site0),
    XY = substr(XYmer, 1, 2)
  ) %>%
  # filter(XYmer == "AA1") %>%
  ggplot(mapping = aes(x = site0, y = Separation)) +
  geom_line(aes(color = XYmer)) +
  # geom_line(aes(color=XYmer), show.legend = TRUE) +
  geom_point(
    mapping = aes(x = site0, y = Separation),
    size = 0.2, color = "black",
    data = ~ subset(., site0 %in% c(-5, -4, 0, 1)), show.legend = FALSE
  ) +
  geom_vhline(
    v = c(-5, -4, 0, 1), h = c(1),
    vlinetype = rep(c("dashed", "dotted"), 2),
    hlinetype = c("dashed"), vlinewidth = rep(0.2, 4),
    hlinewidth = c(0.2), color = rep("black", 4)
  ) +
  # scale_fill_brewer(palette = "Set1") +
  scale_color_manual(values = palette) +
  annotate_sites(
    geom = rep("curve", 4),
    x = c(-19, -15, 10, 14), y = c(
      max(y_breaks) + 0.2, max(y_breaks) + 0.4,
      max(y_breaks) + 0.4, max(y_breaks) + 0.2
    ),
    xend = c(-5, -4, 0, 1), yend = c(
      max(y_breaks) + 0.2, max(y_breaks) + 0.4,
      max(y_breaks) + 0.4, max(y_breaks) + 0.2
    ),
    curvature = c(0, 0, .3, .3), colour = rep("orange", 4),
    units = list(c(2, "mm"), c(2, "mm"), c(2, "mm"), c(2, "mm")),
    label = c("-6", "-5", "-1", "+1")
  ) +
  scale_x_continuous(
    breaks = x_breaks, # limits = limits
    labels = x_labels
  ) +
  scale_y_continuous(
    breaks = y_breaks, # limits = c(min(disp), max(disp))
  ) +
  facet_rep_wrap(~XY,
    ncol = 4, scales = "fixed",
    repeat.tick.labels = TRUE
  ) +
  # facet_rep_wrap(  ~ XYmer, ncol = 4, scales = "fixed",
  #                  repeat.tick.labels = TRUE) +
  ggtitle("") +
  ylab("Separability δ") +
  xlab("Distance from splicing site") +
  # guides(fill=none) +
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
  ) -> gg
directlabels::direct.label(gg)

ggsave(gg + theme(legend.position = "none"),
  filename = "figures/acc_disp_XY0_XY1_4mer.pdf",
  height = 12, width = 12
)
