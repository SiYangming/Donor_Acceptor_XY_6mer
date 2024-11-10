#' check var equality
#'
#' @param x A vector
#' @param y A vector
#' @return Logistic
#' @author Yangming Si
#' @export
check_var_equality <- function(x, y = NULL) {
  # 确保x是一个向量
  if (!is.vector(x)) {
    stop("x must be a vector.")
  }
  
  # 如果y不为空，也确保y是一个向量
  if (y != NULL && !is.vector(y)) {
    stop("y must be a vector.")
  }
  
  # 进行F检验
  f_test_result <- var.test(x, y)
  
  # 判断p值是否大于0.05
  if (f_test_result$p.value > 0.05) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# 函数参数：
# data：数据集
# modes, pos, XYmers：过滤条件
# conservationField：要用于统计测试的保守性分数字段
# outputFileName：输出文件名
performTtests <- function(data, modes, position, XYmers, y_value_column = "Conservation", outputFileName) {
  library(dplyr)
  library(readr)
  
  # 检查输出文件是否已存在
  if (!file.exists(outputFileName)) {
    # 数据过滤和处理
    data %>%
      filter(ss_type %in% modes, 
             pos == position, 
             XYmer %in% XYmers) %>%
      na.omit() %>%
      group_by(XYmer, region) %>%
      group_split() -> lst
    
    # 对每个分组应用统计检验
    stat.test.lst <- lapply(lst, function(x) {
      Group <- unique(x[["ss_type"]])
      res <- tibble()
      
      # 对每对Group进行比较
      for (i in 1:(length(Group) - 1)) {
        for (j in (i + 1):length(Group)) {
          print(paste(unique(x[["XYmer"]]), paste(Group[i], Group[j], sep = "-"), sep = ": "))
          
          g1 <- x[[y_value_column]][x[["ss_type"]] == Group[i]]
          g2 <- x[[y_value_column]][x[["ss_type"]] == Group[j]]
          p <- t.test(g1, g2, paired = TRUE)$p.value
          method <- "T-test"
          
          # 根据样本数据的方差相等性选择检验方法
          # if (check_var_equality(g1, g2)) {
          #   p <- t.test(g1, g2, paired = TRUE)$p.value
          #   method <- "T-test"
          # } else {
          #   p <- pairwise.t.test(c(g1, g2), rep(c(Group[i], Group[j]), each = length(g1)), paired = TRUE)$p.value[, 1]
          #   method <- "Pairwise-T-test"
          # }
          
          # 将结果行添加到结果数据框
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
      return(res)
    })
    
    # 合并结果并调整p值
    stat.test.df <- do.call(bind_rows, stat.test.lst) %>%
      mutate(
        p.format = sapply(p, function(p) format.pval(p, method = "auto")),
        p.adj = p.adjust(p, method = "fdr"),
        p.signif = symnum(round(p.adj,2),
                          cutpoints = c(0, 0.0001, 0.001, 0.01, 0.049, Inf),
                          symbols = c("****", "***", "**", "*", "ns")
        )
      )
    
    # 输出结果到CSV文件
    readr::write_csv(stat.test.df, file = outputFileName)
  } else {
    # 从文件读取之前的结果
    stat.test.df <- readr::read_csv(outputFileName)
  }
  
  return(stat.test.df)
}



library(dplyr)
library(ggpubr)
BarChartWithSignificance <-  function(data, XYmers, region_filter, x_axis_title = "X Axis", y_axis_title = "Y Axis", y_value_column = "Conservation") {
  filtered_data <- data %>%
    filter(XYmer %in% XYmers, region == region_filter)
  
  # 动态计算y_pos的第一个值和y_lim_upper
  y_pos_first <- round(max(tapply(
    filtered_data[[y_value_column]], 
    paste0(filtered_data$XYmer, filtered_data$ss_type), 
    function(x) {mean(x) + sd(x)}
  )) * 1.05, 1)
  
  y_lim_upper <- y_pos_first + 0.6
  
  # 动态确定y_pos的其它值
  y_pos <- seq(from = y_pos_first, by = 0.2, length.out = 3)
  
  # 检查以确保y_lim_upper足以容纳所有y_pos值
  if (max(y_pos) > y_lim_upper) y_lim_upper <- max(y_pos) + 0.2
  
  # 构造图表
  plot <- ggplot(filtered_data, aes(x = .data[["ss_type"]], y = .data[[y_value_column]])) +
    stat_summary(
      mapping = aes(fill = ss_type),
      fun = mean, geom = "bar",
      fun.args = list(mult = 1), width = 0.7
    ) +
    stat_summary(
      fun.data = mean_sdl, fun.args = list(mult = 1),
      geom = "errorbar", width = 0.2
    ) +
    facet_grid(~XYmer) +
    stat_pvalue_manual(
      stat.test.df %>%
        filter(region == region_filter & XYmer %in% XYmers),
      y.position = y_pos,
      label = "p.signif"
    ) +
    labs(
      x = x_axis_title,
      y = y_axis_title
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, y_lim_upper)) +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),
      axis.line = element_line(colour = "black", linewidth = 0.25),
      axis.title = element_text(size = 13, color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
