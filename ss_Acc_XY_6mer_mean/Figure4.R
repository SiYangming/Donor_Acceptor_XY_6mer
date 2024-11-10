rm(list = ls())
source("../code/plot.R")

##### site stat位点统计指标 ######
if (!file.exists("data/RC_filered.rda")) {
  load(file = "data/RC.rda")
  XYmers <- RC$XYmer[grep("0", RC$XYmer, invert = TRUE)]
  XYmer <- unique(XYmers)

  RC_filered <- RC %>%
    filter(ss_type != "All" & (XYmer %in% XYmers)) %>%
    mutate(region = "") %>%
    mutate(region = if_else(between(site, 1, 30), "-50_-21", region)) %>%
    mutate(region = if_else(between(site, 31, 44), "-20_-7", region)) %>%
    mutate(region = if_else(between(site, 45, 51), "-6_+1", region)) %>%
    mutate(region = if_else(between(site, 52, 70), "+2_+20", region)) %>%
    mutate(region = if_else(between(site, 71, 95), "+21_+45", region)) %>%
    na.omit() -> RC_filered
  save(RC_filered, file = "data/RC_filered.rda")
} else {
  load("data/RC_filered.rda")
}

library(dplyr)
library(ggpubr)
modes <- c("Alternative", "Common", "Constitutive")
RC_lst <- RC_filered %>%
  filter(
    ss_type %in% modes,
    region == "-20_-7"
  )
### T-test
if (!file.exists("data/Ass_all_region.csv")) {
  RC_filered %>%
    filter(ss_type %in% modes) %>%
    group_by(XYmer, region) %>%
    group_split() -> lst
  stat.test.lst <- lapply(lst, function(x) {
    Group <- unique(x[["ss_type"]])
    res <- tibble()
    for (i in 1:(length(Group) - 1)) {
      for (j in 2:length(Group)) {
        if (i != j) {
          print(paste(unique(x[["XYmer"]]), paste(Group[i], Group[j], sep = "-"), sep = ": "))
          g1 <- x[["Separation"]][x[["ss_type"]] == Group[i]]
          g2 <- x[["Separation"]][x[["ss_type"]] == Group[j]]
          # print(check_var_equality(g1, g2))
          if (check_var_equality(g1, g2)) {
            p <- t.test(g1, g2, paired = TRUE)$p.value
            method <- "T-test"
          } else {
            p <- pairwise.t.test(c(g1, g2), rep(c(Group[i], Group[j]), each = length(g1)), paired = TRUE)$p.value[, 1]
            method <- "Pairwise-T-test"
          }
          res <- bind_rows(
            res,
            tibble(
              region = unique(x[["region"]]),
              XYmer = unique(x[["XYmer"]]),
              group1 = Group[i], group2 = Group[j],
              p = p, method = method
            )
          )
        }
      }
    }
    # return(res %>%
    #          mutate(p.adj = p.adjust(p, method = "holm")))
    return(res)
  })
  stat.test.df <- do.call(bind_rows, stat.test.lst) %>%
    mutate(
      p.format = sapply(p, pValueFormat),
      p.adj = p.adjust(p, method = "fdr"),
      p.signif = symnum(p,
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
        symbols = c("****", "***", "**", "*", "ns")
      )
    )
  readr::write_csv(stat.test.df, file = "data/Ass_all_region.csv")
} else {
  stat.test.df <- readr::read_csv("data/Ass_all_region.csv")
}
table(stat.test.df$region, stat.test.df$method)
# Pairwise-T-test T-test
# -20_-7                5     43
# -50_-21              14     34
# -6_+1                 0     48
# +2_+20               16     32
# +21_+45              34     14

RC_lst$ss_type <- tolower(RC_lst$ss_type)
my_comparisons <- list(
  c("alternative", "constitutive"),
  c("alternative", "common"),
  c("constitutive", "common")
)
g1 <- c("AA1", "GG1", "AG1", "GA1", "AT1", "TA1", "AC1", "CA1")
g2 <- c("CT1", "TC1")
g3 <- c("CG1", "GC1", "GT1", "TG1", "CC1", "TT1")

gg1 <- ggplot(
  RC_lst %>%
    filter(XYmer %in% g1),
  aes(ss_type, Separation)
) +
  # 添加柱
  stat_summary(
    mapping = aes(fill = ss_type),
    fun = mean, geom = "bar",
    fun.args = list(mult = 1), width = 0.7
  ) +
  # 添加误差线
  stat_summary(
    fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "errorbar", width = 0.2
  ) +
  facet_grid(~XYmer) +
  # stat_compare_means(comparisons = my_comparisons,
  #                    label = "p.signif",
  #                    label.y = c(1.7, 1.9, 2.1),
  #                    method = "t.test",
  #                    paired = TRUE,
  #                    p.adjust.method = "fdr") +
  stat_pvalue_manual(
    stat.test.df %>%
      filter(region == "-20_-7" &
        XYmer %in% g1),
    y.position = c(1.15, 1.35, 1.55),
    label = "p.signif"
  ) +
  # 添加x，y轴名
  labs(
    x = "Acceptor splicing mode",
    y = "Separability δ"
  ) +
  # 坐标轴延伸，确保图形元素覆盖至坐标
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.6)) +
  # 主题类型
  theme_classic() +
  # 设置主题
  theme(
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25), # 填充框内主题颜色，边框颜色和边框线条粗细
    axis.line = element_line(colour = "black", linewidth = 0.25), # x,y轴颜色，粗细
    axis.title = element_text(size = 13, color = "black"), # x,y轴名设置
    axis.text = element_text(size = 12, color = "black"), # x,y轴文本设置
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# ggsave(gg, filename = "figures/all_BestXYGroup_UpRegion.pdf", width = 12, height = 6)

gg2 <- ggplot(
  RC_lst %>%
    filter(XYmer %in% g2),
  aes(ss_type, Separation)
) +
  # 添加柱
  stat_summary(
    mapping = aes(fill = ss_type),
    fun = mean, geom = "bar",
    fun.args = list(mult = 1), width = 0.7
  ) +
  # 添加误差线
  stat_summary(
    fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "errorbar", width = 0.2
  ) +
  facet_grid(~XYmer) +
  stat_pvalue_manual(
    stat.test.df %>%
      filter(region == "-20_-7" &
        XYmer %in% g2),
    y.position = c(1.8, 2.0, 2.2),
    label = "p.signif"
  ) +
  # 添加x，y轴名
  labs(x = "Acceptor splicing mode", y = "Separability δ") +
  # 坐标轴延伸，确保图形元素覆盖至坐标
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.3)) +
  # 主题类型
  theme_classic() +
  # 设置主题
  theme(
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25), # 填充框内主题颜色，边框颜色和边框线条粗细
    axis.line = element_line(colour = "black", linewidth = 0.25), # x,y轴颜色，粗细
    axis.title = element_text(size = 13, color = "black"), # x,y轴名设置
    axis.text = element_text(size = 12, color = "black"), # x,y轴文本设置
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
# ggsave(gg, filename = "figures/all_GoodXYGroup_UpRegion.pdf", width = 6, height = 6)

gg3 <- ggplot(
  RC_lst %>%
    filter(XYmer %in% g3),
  aes(ss_type, Separation)
) +
  # 添加柱
  stat_summary(
    mapping = aes(fill = ss_type),
    fun = mean, geom = "bar",
    fun.args = list(mult = 1), width = 0.7
  ) +
  # 添加误差线
  stat_summary(
    fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "errorbar", width = 0.2
  ) +
  facet_grid(~XYmer) +
  stat_pvalue_manual(
    stat.test.df %>%
      filter(region == "-20_-7" &
        XYmer %in% g3),
    y.position = c(2.35, 2.55, 2.75),
    label = "p.signif"
  ) +
  # 添加x，y轴名
  labs(x = "Acceptor splicing mode", y = "Separability δ") +
  # 坐标轴延伸，确保图形元素覆盖至坐标
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.9)) +
  # 主题类型
  theme_classic() +
  # 设置主题
  theme(
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25), # 填充框内主题颜色，边框颜色和边框线条粗细
    axis.line = element_line(colour = "black", linewidth = 0.25), # x,y轴颜色，粗细
    axis.title = element_text(size = 13, color = "black"), # x,y轴名设置
    axis.text = element_text(size = 12, color = "black"), # x,y轴文本设置
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
# ggsave(gg, filename = "figures/all_BadXYGroup_UpRegion.pdf", width = 6, height = 6)
# library(patchwork)
# gg <- gg1/(gg2 + gg3) + plot_annotation(tag_levels = "A")
library(cowplot)
up_row <- plot_grid(gg1, labels = "A")
bottom_row <- plot_grid(gg2, gg3, rel_widths = c(2, 6), labels = c("B", "C"))
gg <- plot_grid(up_row, bottom_row, ncol = 1)
ggsave(gg, filename = "figures/all_UpRegion.pdf", width = 20, height = 10)
