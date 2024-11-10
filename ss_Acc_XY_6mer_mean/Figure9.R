library(dplyr)
library(ggpubr)
load(file = "data/RC_filered.rda")
downregion_2_20 <- c("GA1", "GC1", "TC1", "TT1")
RC_filered$ss_type <- tolower(RC_filered$ss_type)
RC_lst <- RC_filered %>%
  filter((ss_type %in% c("normal", "exonic", "intronic")) &
    region %in% c("+2_+20", "+21_+45")) %>%
  dplyr::filter(!(region == "+21_+45" & (XYmer %in% downregion_2_20))) %>%
  mutate(region = if_else(region %in% c("+2_+20", "+21_+45"), "+2_+45", region))
stat.test.df <- readr::read_csv("data/AltAss_all_region.csv")
my_comparisons <- list(
  c("exonic", "intronic"),
  c("exonic", "normal"),
  c("intronic", "normal")
)
G1 <- c("AG1", "CA1", "CC1", "CG1", "GA1", "GC1")
G2 <- c("GT1", "TG1")
G3 <- c("TA1", "TC1", "TT1")
G4 <- c("AA1", "AC1", "AT1", "CT1", "GG1")

gg1 <- ggplot(
  RC_lst %>%
    dplyr::filter(region == "+2_+45" &
      XYmer %in% G1),
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
      dplyr::filter(region == "+2_+45" &
        XYmer %in% G1),
    y.position = c(1.3, 1.5, 1.7),
    label = "p.signif"
  ) +
  # 添加x，y轴名
  labs(x = "Acceptor splicing mode", y = "Separability δ") +
  # 坐标轴延伸，确保图形元素覆盖至坐标
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.8)) +
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
    dplyr::filter(region == "+2_+45" &
      XYmer %in% G2),
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
      dplyr::filter(region == "+2_+45" &
        XYmer %in% G2),
    y.position = c(1.3, 1.5, 1.7),
    label = "p.signif"
  ) +
  # 添加x，y轴名
  labs(x = "Acceptor splicing mode", y = "Separability δ") +
  # 坐标轴延伸，确保图形元素覆盖至坐标
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.8)) +
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

gg3 <- ggplot(
  RC_lst %>%
    dplyr::filter(region == "+2_+45" &
      XYmer %in% G3),
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
      dplyr::filter(region == "+2_+45" &
        XYmer %in% G3),
    y.position = c(1.3, 1.5, 1.7),
    label = "p.signif"
  ) +
  # 添加x，y轴名
  labs(x = "Acceptor splicing mode", y = "Separability δ") +
  # 坐标轴延伸，确保图形元素覆盖至坐标
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.8)) +
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

gg4 <- ggplot(
  RC_lst %>%
    dplyr::filter(region == "+2_+45" &
      XYmer %in% G4),
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
      dplyr::filter(region == "+2_+45" &
        XYmer %in% G4),
    y.position = c(1.3, 1.5, 1.7),
    label = "p.signif"
  ) +
  # 添加x，y轴名
  labs(x = "Acceptor splicing mode", y = "Separability δ") +
  # 坐标轴延伸，确保图形元素覆盖至坐标
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.8)) +
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

library(cowplot)
up_row <- plot_grid(gg1, gg2, rel_widths = c(6, 2), labels = c("A", "B"))
bottom_row <- plot_grid(gg3, gg4, rel_widths = c(3, 5), labels = c("C", "D"))
gg <- plot_grid(up_row, bottom_row, ncol = 1)
# gg <- gg1/gg2 + plot_annotation(tag_levels = "A")
ggsave(gg,
  filename = "figures/alt_subset_DownRegion.pdf",
  width = 18, height = 9
)
