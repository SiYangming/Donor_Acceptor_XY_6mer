rm(list = ls())
source("../code/kmers-function.R")
paths <- "data"
freqs_filenames <- list.files(paths,
  pattern = "6mer_freqs.csv.gz",
  full.names = TRUE
)
#### Test
rt <- read_csv(freqs_filenames[2])
freqs <- kmer_classification(rt[[1]], kmer = "CG")
mer_classify <- lapply(kmer(2), function(x) {
  kmer_classification(rt[[1]], kmer = x)
})
mer_classify_df <- do.call(cbind, mer_classify)

mer <- "6mer"
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
    add_column(do.call(cbind, mer_classify), .after = "6mer")
  outfile <- gsub("_freqs.csv.gz", "_classification.csv.gz", filename)
  write_csv(mer_classify_df, file = outfile, col_names = T)
}

clf_filenames <- list.files(paths,
  pattern = "6mer_classification.csv.gz",
  full.names = TRUE
)
## Test
clf_df <- read_csv(clf_filenames[1], col_names = T)
### Run
pattern <- "data\\/|_6mer_classification.csv.gz"
if (!file.exists("data/plotdata.rda")) {
  clf_df_lst <- lapply(clf_filenames, function(x) {
    type <- gsub(pattern, "", x)
    rt <- read_csv(x,
      col_types = as.list(c(rep("c", 17), rep("d", 95))),
      col_names = T
    )
    GG_kmer_RidgesData(rt) %>%
      mutate(type = type)
  })

  plotdata <- do.call(rbind, clf_df_lst)
  save(plotdata, file = "data/plotdata.rda")
} else {
  load(file = "data/plotdata.rda")
}

# rm(list = ls())
# source("code/kmers-function.R")
# load(file = "data/sel/plotdata.rda")
# sites = 46:55
# plotlist <- plotdata %>% 
#   filter(site %in% sites) %>%
#   group_split(site)
# 
# 
# ggplotobj <- GG_kmer_RidgesPlotBatch(plotlist)
# names(ggplotobj) <- sites
# # save(ggplotobj, file = '20221024/ggplotobj.rda')
# dir.create("figs/")
# for (site in names(ggplotobj)) {
#   ggsave(ggplotobj[[site]], filename = paste0("figs/alt_con_norm/",
#                                               site, ".pdf"), 
#          width = 16, height = 26)
# }

if (!file.exists("data/RC_matrix.rda")) {
  RC_matrix <- plotdata %>%
    dplyr::group_by(type, site, value, .drop = FALSE) %>%
    # group_split()
    dplyr::summarise(
      mean = mean(frequence, na.rm = TRUE),
      sd = sd(frequence, na.rm = TRUE)
    )
  save(RC_matrix, file = "data/RC_matrix.rda")
} else {
  load("data/RC_matrix.rda")
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
# [1] 0.0000 2.6805
range(RC$Conservation, na.rm = TRUE)
# [1] 0.0000 1.9752
if (!file.exists("data/RC.rda")) {
  modes <- do.call(rbind, str_split(RC$type, "_"))
  RC <- RC %>%
    mutate(pos = modes[, 2], ss_type = modes[, 1])

  RC <- RC %>%
    add_row(
      type = rep(rep(unique(RC$type), each = 32), 5),
      site = rep(96:100, each = 224),
      XYmer = rep(c(paste0(kmer(2), 0), paste0(kmer(2), 1)), 35),
      Separation = NA, Conservation = NA,
      pos = unique(RC$pos),
      ss_type = rep(rep(unique(RC$ss_type),
        each = 32
      ), 5)
    )
  save(RC, file = "data/RC.rda")
} else {
  load(file = "data/RC.rda")
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
y_breaks <- seq(0, 2.8, 0.5)
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

ggsave(directlabels::direct.label(gg),
  filename = "figures/acc_disp_XY0_XY1.pdf",
  height = 12, width = 12
)
