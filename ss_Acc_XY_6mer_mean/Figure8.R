library(dplyr)
library(ggpubr)
load(file = "data/RC_filered.rda")
RC_filered$ss_type <- tolower(RC_filered$ss_type)
RC_lst <- RC_filered %>%
  filter((ss_type %in% c("normal", "exonic", "intronic")) &
    region == "-6_+1")
stat.test.df <- readr::read_csv("data/AltAss_all_region.csv")
my_comparisons <- list(
  c("exonic", "intronic"),
  c("exonic", "normal"),
  c("intronic", "normal")
)

G1 <- c("CG1", "AA1", "CA1", "CC1", "GC1")
G2 <- c("TT1", "AT1", "TA1", "TC1")
G3 <- c("AC1", "AG1", "CT1", "GA1", "GG1", "GT1", "TG1")

gg1 <- ggplot(
  RC_lst %>%
    filter(region == "-6_+1" &
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
      filter(region == "-6_+1" &
        XYmer %in% G1),
    y.position = c(2.0, 2.2, 2.4),
    label = "p.signif"
  ) +
  # 添加x，y轴名
  labs(x = "Acceptor splicing mode", y = "Separability δ") +
  # 坐标轴延伸，确保图形元素覆盖至坐标
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5)) +
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
    filter(region == "-6_+1" &
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
      filter(region == "-6_+1" &
        XYmer %in% G2),
    y.position = c(2.0, 2.2, 2.4),
    label = "p.signif"
  ) +
  # 添加x，y轴名
  labs(x = "Acceptor splicing mode", y = "Separability δ") +
  # 坐标轴延伸，确保图形元素覆盖至坐标
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5)) +
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

# gg3 <- ggplot(RC_lst %>%
#                 filter(region == "-6_+1" &
#                          XYmer %in% G3),
#               aes(ss_type, Separation)) +
#   # 添加柱
#   stat_summary(mapping=aes(fill = ss_type),
#                fun=mean, geom = "bar",
#                fun.args = list(mult=1), width=0.7) +
#   # 添加误差线
#   stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
#                geom="errorbar", width=0.2) +
#   facet_grid(~XYmer) +
#   stat_pvalue_manual(stat.test.df %>%
#                        filter(region == "-6_+1" &
#                                 XYmer %in% G3),
#                      y.position = c(1.6, 1.8, 2.0),
#                      label = "p.signif") +
#   # 添加x，y轴名
#   labs(x = "Acceptor splicing mode",y = "Separability δ")+
#   # 坐标轴延伸，确保图形元素覆盖至坐标
#   scale_y_continuous(expand = c(0,0), limits = c(0,2.5)
#   ) +
#   # 主题类型
#   theme_classic() +
#   # 设置主题
#   theme(panel.background=element_rect(fill="white",colour="black",linewidth=0.25), # 填充框内主题颜色，边框颜色和边框线条粗细
#         axis.line=element_line(colour="black",linewidth=0.25), # x,y轴颜色，粗细
#         axis.title=element_text(size=13,color="black"), # x,y轴名设置
#         axis.text = element_text(size=12,color="black"), # x,y轴文本设置
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position="none")

library(cowplot)
# up_row <- plot_grid(gg1, gg2, rel_widths = c(2,7), labels = c('A', 'B'))
# bottom_row <- plot_grid(gg3, labels = c('C'))
# gg <- plot_grid(up_row, bottom_row, ncol = 1)
gg <- plot_grid(gg1, gg2, ncol = 1, labels = c("A", "B"))
ggsave(gg,
  filename = "figures/alt_subset_CenterRegion.pdf",
  width = 12, height = 7
)
