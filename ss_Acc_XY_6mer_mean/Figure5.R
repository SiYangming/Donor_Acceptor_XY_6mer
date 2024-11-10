library(dplyr)
library(ggpubr)
load(file = "data/RC_filered.rda")
modes <- c("Alternative", "Common", "Constitutive")
RC_lst <- RC_filered %>%
  filter(
    ss_type %in% modes,
    region %in% c("+2_+20", "+21_+45")
  ) %>%
  mutate(region = "+2_+45")

RC_lst %>%
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
# readr::write_csv(stat.test.df, file = "data/Ass_downstream_region.csv")

RC_lst$ss_type <- tolower(RC_lst$ss_type)
my_comparisons <- list(
  c("alternative", "constitutive"),
  c("alternative", "common"),
  c("constitutive", "common")
)

G1 <- c("AA1", "AC1", "AG1", "CG1", "CT1", "GA1", "GC1", "TC1", "TG1")
G2 <- c("AT1", "CA1", "CC1", "GG1", "GT1", "TA1", "TT1")
gg1 <- ggplot(
  RC_lst %>%
    filter(XYmer %in% G1),
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
      filter(XYmer %in% G1),
    y.position = c(1.4, 1.6, 1.8),
    label = "p.signif"
  ) +
  # 添加x，y轴名
  labs(x = "Acceptor splicing mode", y = "Separability δ") +
  # 坐标轴延伸，确保图形元素覆盖至坐标
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.9)) +
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
gg2 <- ggplot(
  RC_lst %>%
    filter(XYmer %in% G2),
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
      filter(XYmer %in% G2),
    y.position = c(1.4, 1.6, 1.8),
    label = "p.signif"
  ) +
  # 添加x，y轴名
  labs(x = "Acceptor splicing mode", y = "Separability δ") +
  # 坐标轴延伸，确保图形元素覆盖至坐标
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.9)) +
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
library(patchwork)
gg <- gg1 / gg2 + plot_annotation(tag_levels = "A")
ggsave(gg, filename = "figures/all_DownRegion.pdf", width = 20, height = 10)
