rm(list = ls())
source("../code/kmers-function.R")
load(file = "../data/sites_RC.rda")
source("../code/layers.R")
modes <- "All"
poss <- "Acceptor"
x_title <- "Distance from splicing site"
y_title <- "Dispersion ρ"
y_column <- "Conservation"
plot.file <- "figures/acc_sd_XY0_XY1.pdf"

RC %>%
  dplyr::filter(pos == poss &
                  ss_type == modes) %>%
  View()

##### Figure 2:  Overall Acceptor Splicing Signal 背景值XY1&XY0 #####
table(RC$ss_type)
table(RC$pos)

### Overall Alternative Acceptor 背景值XY1

palette <- rep(c("#F8786F", "#609DFF"), 16)
x_breaks <- (c(seq(1, 49, 10), 45:46, 50:51, seq(60, 100, 10)) - 50)
x_labels <- if_else(x_breaks <= 0, x_breaks - 1, x_breaks)
x_labels <- if_else(x_labels > 0,
  paste0("+", as.character(x_labels)),
  as.character(x_labels)
)
x_labels <- if_else(x_labels %in% c("-6", "-5", "-1", "+1"), "", x_labels)
RC %>%
  arrange(type, XYmer) %>%
  filter(pos == poss, ss_type == modes) %>%
  select(all_of(y_column)) %>%
  range(na.rm = TRUE)
# [1] 0.0000 1.9074
y_breaks <- seq(0, 2, 0.5)
RC %>%
  arrange(type, XYmer) %>%
  dplyr::filter(pos == poss & ss_type == modes) %>%
  group_by(XYmer, site) %>%
  dplyr::mutate(
    site0 = site - 50,
    site1 = if_else(site0 <= 0, site0 - 1, site0),
    XY = substr(XYmer, 1, 2)
  ) %>%
  # filter(XYmer == "AA1") %>%
  ggplot(mapping = aes(x = site0, y = .data[[y_column]])) +
  geom_line(aes(color = XYmer)) +
  # geom_line(aes(color=XYmer), show.legend = TRUE) +
  geom_point(
    mapping = aes(x = site0, y = .data[[y_column]]),
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
  ylab(y_title) +
  xlab(x_title) +
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
  filename = plot.file,
  height = 12, width = 12
)

# lack_XY <- c("AA1", "GG1", "CC1", "AG1", "GA1", "AT1", "TA1", "AC1", "CA1", "GC1")
# rich_XY <- c("TC1", "TT1")
# other_XY <- c("GT1", "TG1", "CG1", "CT1")
# 
# ms <- function(x) {
#   paste(
#     unique(x[["motifs"]]),
#     paste0(
#       round(mean(x[["Conservation"]], na.rm = TRUE), digits = 2),
#       "±",
#       round(sd(x[["Conservation"]], na.rm = TRUE), digits = 2)
#     )
#   )
# }
# 
# RC %>%
#   mutate(
#     motifs = "",
#     motifs = if_else(XYmer %in% lack_XY, "lack", motifs),
#     motifs = if_else(XYmer %in% rich_XY, "rich", motifs),
#     motifs = if_else(XYmer %in% other_XY, "other", motifs)
#   ) %>%
#   group_by(motifs) %>%
#   group_split() %>%
#   map(~ ms(.x))
# # [[1]]
# # [1] " 1.01±0.15"
# #
# # [[2]]
# # [1] "lack 0.77±0.33"
# #
# # [[3]]
# # [1] "other 0.77±0.36"
# #
# # [[4]]
# # [1] "rich 1.05±0.37"
# 
# 
# (RC %>%
#   dplyr::filter(XYmer == "CG1") %>%
#   dplyr::select(Conservation))[[1]] %>%
#   mean(na.rm = TRUE) %>%
#   round(digits = 2)
# 
# map(~ ms(.x))
