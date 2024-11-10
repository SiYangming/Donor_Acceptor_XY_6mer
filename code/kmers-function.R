library(tidyverse)
library(reshape2)
# library(ggTimeSeries)
library(ggseas)
library(WriteXLS)
library(directlabels)
library(lemon)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(randomcoloR)
library(lattice)


kmer <- function(k = 2) {
  unique(sapply(data.frame(combn(rep(c("A", "T", "C", "G"), k), k)),
    FUN = function(x) {
      paste(x, collapse = "")
    }
  ))
}

kmer_classification <- function(strings, kmer = "CG", group = 3, collapse = TRUE, ROWNAME = NULL) {
  if (collapse) {
    kmer_counts <- stringr::str_count(strings, pattern = kmer)
    kmer_counts[kmer_counts >= (group - 1)] <- group - 1
    clf <- data.frame(stringr::str_c(kmer, kmer_counts))
  } else {
    clf <- data.frame(stringr::str_c(
      kmer,
      stringr::str_count(strings,
        pattern = kmer
      )
    ))
  }

  if (is.null(ROWNAME)) {
    rownames(clf) <- strings
  } else {
    rownames(clf) <- ROWNAME
  }

  colnames(clf) <- kmer
  return(clf)
}

WriteXYmerPPM <- function(Sixmer, Twomer, XY = "CG",
                          group = "1", path = "data/XY1/") {
  freqs_lst <- lapply(1:ncol(Sixmer), function(x) {
    logi <- str_detect(
      (kmer_classification(Sixmer[[x]],
        kmer = XY,
        group = 2,
        ROWNAME = 1:length(Sixmer[[x]])
      ))[[1]],
      pattern = paste0(XY, group)
    )
    if (sum(logi) == 0) {
      res <- tibble(kmer = kmer(2), 0)
    } else {
      res <- Twomer[[x]][logi] %>%
        table() %>%
        as.data.frame()
    }
    colnames(res) <- c("kmer", x)
    res
  })

  for (i in 1:length(freqs_lst)) {
    if (i == 1) {
      freqs_df <- freqs_lst[[i]]
    } else {
      freqs_df <- freqs_df %>%
        full_join(freqs_lst[[i]])
    }
  }
  pfm2ppm(tidyToy::DF2matrix(freqs_df)) %>%
    as.data.frame() %>%
    rownames_to_column(var = "kmer") %>%
    write_csv(file = paste0(path, XY, group, ".csv"))
}


atcg2num <- function(df) {
  as.data.frame(do.call(rbind, lapply(df, function(x) {
    as.numeric(gsub("[ATCG]", "", x))
  })))
}

GG_kmer_RidgesData <- function(clf_df, prefix = "6mer") {
  clf_df_melted <- reshape2::melt(clf_df, , value.name = "frequence")
  colnames(clf_df_melted)[ncol(clf_df_melted) - 1] <- "site"
  lst <- clf_df_melted %>%
    add_column(kmer = "kmer", .after = prefix) %>%
    group_by(site, .drop = FALSE) %>%
    # nest() %>%
    group_map(~ reshape2::melt(.x[-1], id = "frequence"))

  plotdata_lst <- lapply(1:length(lst), function(x) {
    lst[[x]] %>%
      mutate(site = x)
  })
  PlotData <- do.call(rbind, plotdata_lst)

  return(PlotData)
}

GG_kmer_RidgesPlot <- function(PlotData, facet = TRUE) {
  set.seed(123456789)
  palette <- randomColor(count = 60) # 随机生成60种颜色，其实里面有重复的
  palette <- distinctColorPalette(60)
  ridges <- ggplot(PlotData, aes(x = frequence, y = variable, fill = value, group = value)) +
    scale_fill_manual(values = palette) +
    geom_density_ridges_gradient() +
    # scale_x_continuous(limits = c(-500, 10000)) +
    xlab("") +
    ylab("") +
    guides(fill = guide_legend(title = "XY-mer")) +
    theme_ridges()
  if (facet) {
    ridges + facet_wrap(~type, ncol = 2, scales = "free_y")
    # facet_grid(. ~ type)
  } else {
    ridges
  }
}

GG_kmer_RidgesPlotBatch <- function(plotlist) {
  RidgesPlotList <- list()
  for (i in 1:length(plotlist)) {
    RidgesPlotList[[i]] <- GG_kmer_RidgesPlot(plotlist[[i]])
  }
  return(RidgesPlotList)
}

calc_Ratio <- function(sites_RC, digits = 2, invert = TRUE) {
  kmer <- sites_RC %>% filter(value == "kmer")
  if (invert) {
    # Observations divided by overall
    sites_RC %>%
      filter(value != "kmer") %>%
      mutate(
        Separation = round(mean / kmer[["mean"]], digits),
        Conservation = round(sd / kmer[["sd"]], digits)
      )
  } else {
    # Overall divided by observations
    sites_RC %>%
      filter(value != "kmer") %>%
      mutate(
        Separation = round(kmer[["mean"]] / mean, digits),
        Conservation = round(kmer[["sd"]] / sd, digits)
      )
  }
}

plot_connect <- function(plotdata, name = "disp", range = NULL) {
  library(patchwork)
  p <- NULL
  if (is.null(range)) {
    range <- 1:length(plotdata)
  }
  for (i in range) {
    if (i == 1) {
      p <- plotdata[[i]][[name]]
    } else {
      p <- p + plotdata[[i]][[name]]
    }
  }
  return(p)
}

plot_connect2 <- function(plotdata, range = NULL) {
  library(patchwork)
  p <- NULL
  if (is.null(range)) {
    range <- 1:length(plotdata)
  }
  for (i in range) {
    if (i == 1) {
      p <- plotdata[[i]]
    } else {
      p <- p + plotdata[[i]]
    }
  }
  return(p)
}

# xy_mer_disp_sd_plot1 <- function(plotdata){
#   set.seed(123456789)
#   palette <- randomColor(count = 60)  #随机生成60种颜色，其实里面有重复的
#   palette <- distinctColorPalette(60)
#   #plotdata_melted <- melt(plotdata, id = "value")
#   #colnames(plotdata_melted) <- c("XYmer", "site", "disp")
#   #plotdata_plot[[name]] <- as.numeric(plotdata_plot[[name]])
#   plotdata[["site"]] <- as.numeric(plotdata[["site"]])
#   plotdata <- plotdata %>%
#     mutate(type = substr(XYmer, 1, 2), site0 = site -300)
#
#   disp_p = ggplot(data = plotdata, mapping = aes(x = site0, y = log2(disp), fill = XYmer)) +
#     geom_line(aes(color=XYmer), show.legend = TRUE) +
#     stat_seas() +
#     scale_color_manual(values = palette) +
#     facet_wrap( ~ type, ncol =4) +
#     #scale_y_continuous(breaks = seq(0,8,2), limits = c(-2, 8), labels = c(2^seq(0,8,2))) +
#     scale_x_continuous(breaks = seq(-300,300,50),limits = c(-300, 350)) +
#     #scale_y_log10() +
#     ggtitle("") +
#     ylab("log2(disp)") +
#     xlab("site") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, size=6,face="bold",
#                                      vjust=1, hjust=1))
#
#   sd_p = ggplot(data = plotdata, mapping = aes(x = site0, y = log2(sd), fill = XYmer)) +
#     geom_line(aes(color=XYmer), show.legend = TRUE) +
#     stat_seas() +
#     scale_color_manual(values = palette) +
#     facet_wrap( ~ type, ncol =4) +
#     #scale_y_continuous(breaks = seq(0,10,2.5), limits = c(-2.0, 10.0), labels = round(c(2^seq(0,10,2.5)))) +
#     scale_x_continuous(breaks = seq(-300,300,50),limits = c(-300, 350)) +
#     #scale_y_log10() +
#     ggtitle("") +
#     ylab("log2(sd)") +
#     xlab("site") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, size=6, face="bold",
#                                      vjust=1, hjust=1))
#
#   return(list(disp = direct.label(disp_p, method = list("right.polygons", cex=0.4)),
#               sd = direct.label(sd_p, method = list("right.polygons", cex=0.3))))
# }

# xy_mer_disp_sd_plot2 <- function(plotdata){
#   set.seed(123456789)
#   palette <- randomColor(count = 60)  #随机生成60种颜色，其实里面有重复的
#   palette <- distinctColorPalette(60)
#   plotdata[["site"]] <- as.numeric(plotdata[["site"]])
#   max_site <- max(plotdata$site)
#   plotdata <- plotdata %>%
#     mutate(type = substr(XYmer, 1, 2), site0 = site -(max_site/2)) %>%
#     mutate(site = if_else(site0 <= 0, site0-1, site0))
#   start <- min(plotdata$site)
#   stop <- max(plotdata$site)
#
#   disp_p = ggplot(data = plotdata, mapping = aes(x = site, y = log2(disp), fill = XYmer)) +
#     geom_line(aes(color=XYmer), show.legend = TRUE) +
#     stat_seas() +
#     scale_color_manual(values = palette) +
#     facet_rep_wrap( ~ type, ncol =4,  repeat.tick.labels = TRUE) +
#     #facet_wrap( ~ type, ncol =4) +
#     scale_x_continuous(breaks = seq(start,stop,5)[-c((stop-start)/(5*2)+1)], limits = c(start,stop+15)) +
#     ggtitle("") +
#     ylab("log2(disp)") +
#     xlab("site") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, size=8,face="bold",
#                                      vjust=1, hjust=1))
#
#   sd_p = ggplot(data = plotdata, mapping = aes(x = site0, y = log2(sd), fill = XYmer)) +
#     geom_line(aes(color=XYmer), show.legend = TRUE) +
#     stat_seas() +
#     scale_color_manual(values = palette) +
#     facet_rep_wrap( ~ type, ncol =4,  repeat.tick.labels = TRUE) +
#     #facet_wrap( ~ type, ncol =4) +
#     scale_x_continuous(breaks = seq(start,stop,5)[-c((stop-start)/(5*2)+1)], limits = c(start,stop+15)) +
#     ggtitle("") +
#     ylab("log2(sd)") +
#     xlab("site") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, size=8, face="bold",
#                                      vjust=1, hjust=1))
#
#   return(list(disp = direct.label(disp_p, method = list("right.polygons")),
#               sd = direct.label(sd_p, method = list("right.polygons"))))
# }

# xy_mer_disp_sd_plot_20 <- function(plotdata){
#   set.seed(123456789)
#   palette <- randomColor(count = 60)  #随机生成60种颜色，其实里面有重复的
#   palette <- distinctColorPalette(60)
#   plotdata[["site"]] <- as.numeric(plotdata[["site"]])
#   stream <- round(median(plotdata[["site"]]))
#   plotdata <- plotdata %>%
#     mutate(type = substr(XYmer, 1, 2), site0 = site - stream) %>%
#     mutate(site = if_else(site0 <= 0, site0 - 1, site0))
#   start <- min(plotdata$site)
#   stop <- max(plotdata$site)
#
#   disp_p = ggplot(data = plotdata, mapping = aes(x = site, y = log2(disp), fill = XYmer)) +
#     geom_line(aes(color=XYmer), show.legend = TRUE) +
#     geom_vline(xintercept = -1,color = "black", linetype="dashed", linewidth = 0.2) +
#     stat_seas() +
#     scale_color_manual(values = palette) +
#     facet_rep_wrap( ~ type, ncol =4,  repeat.tick.labels = TRUE) +
#     #facet_wrap( ~ type, ncol =4) +
#     scale_x_continuous(breaks = seq(start,stop,2)[-c((stop-start)/(2*2)+1)], limits = c(start,stop+10)) +
#     ggtitle("") +
#     ylab("log2(disp)") +
#     xlab("site") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, size=5,face="bold",
#                                      vjust=1, hjust=1))
#
#   sd_p = ggplot(data = plotdata, mapping = aes(x = site0, y = log2(sd), fill = XYmer)) +
#     geom_line(aes(color=XYmer), show.legend = TRUE) +
#     geom_vline(xintercept = -1,color = "black", linetype="dashed", linewidth = 0.2) +
#     stat_seas() +
#     scale_color_manual(values = palette) +
#     facet_rep_wrap( ~ type, ncol =4,  repeat.tick.labels = TRUE) +
#     #facet_wrap( ~ type, ncol =4) +
#     scale_x_continuous(breaks = seq(start,stop,2)[-c((stop-start)/(2*2)+1)], limits = c(start,stop+10)) +
#     ggtitle("") +
#     ylab("log2(sd)") +
#     xlab("site") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, size=5, face="bold",
#                                      vjust=1, hjust=1))
#
#   return(list(disp = direct.label(disp_p, method = list("right.polygons")),
#               sd = direct.label(sd_p, method = list("right.polygons"))))
# }

xy_mer_disp_sd_plot_10 <- function(plotdata, split_by = "samples", disp_y_breaks = NULL, sd_y_breaks = NULL) {
  set.seed(123456789)
  palette <- randomColor(count = 60) # 随机生成60种颜色，其实里面有重复的
  palette <- distinctColorPalette(60)
  plotdata[["site"]] <- as.numeric(plotdata[["site"]])
  stream <- round(median(plotdata$site))
  if (split_by == "samples") {
    plotdata <- plotdata %>%
      mutate(type = substr(XYmer, 1, 2), site0 = site - stream) %>%
      mutate(site = if_else(site0 <= 0, site0 - 1, site0))
  } else if (split_by == "XYmer") {
    plotdata <- plotdata %>%
      mutate(site0 = site - stream) %>%
      mutate(site = if_else(site0 <= 0, site0 - 1, site0))
  } else {
    stop("Only support samples or XYmer!")
  }

  start <- -10
  stop <- 10

  # breaks <- c(-20,seq(start,stop,1),20)[-c((stop-start+2)/2+1)]
  x_breaks <- plotdata[["site0"]]
  if (is.null(disp_y_breaks)) {
    disp_y_breaks <- seq(from = 0, to = 2.4, by = 0.2)
  }
  disp_y_limits <- c(min(disp_y_breaks), max(disp_y_breaks))
  # x_limits <- c(min(plotdata[["site"]]),max(plotdata[["site"]]))
  # y_limits <- c(min(y_breaks), max(y_breaks))
  x_labels <- as.character(plotdata[["site"]])
  x_labels[grep("-", x_labels, invert = TRUE)] <- paste0("+", x_labels[grep("-", x_labels, invert = TRUE)])

  disp_p <- ggplot(data = plotdata, mapping = aes(x = site0, y = Separation, fill = XYmer)) +
    geom_point(aes(color = XYmer), show.legend = FALSE) +
    geom_line(aes(color = XYmer), show.legend = FALSE) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.2) +
    geom_vline(xintercept = 1, color = "black", linetype = "dashed", linewidth = 0.2) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.2) +
    scale_color_manual(values = palette) +
    facet_rep_wrap(~type,
      ncol = 4, scales = "fixed",
      repeat.tick.labels = TRUE
    ) +
    scale_x_continuous(
      breaks = x_breaks, # limits = limits
      labels = x_labels
    ) +
    scale_y_continuous(breaks = disp_y_breaks, limits = disp_y_limits) +
    ggtitle("") +
    ylab("Separation") +
    xlab("site") +
    # guides(fill=none) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 90, # size=5,
      face = "bold",
      vjust = 0.4, hjust = 1
    ))
  if (is.null(sd_y_breaks)) {
    sd_y_breaks <- seq(from = 0, to = 2.4, by = 0.2)
  }
  sd_y_limits <- c(min(sd_y_breaks), max(sd_y_breaks))
  sd_p <- ggplot(data = plotdata, mapping = aes(x = site0, y = Conservation, fill = XYmer)) +
    geom_point(aes(color = XYmer), show.legend = FALSE) +
    geom_line(aes(color = XYmer), show.legend = FALSE) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.2) +
    geom_vline(xintercept = 1, color = "black", linetype = "dashed", linewidth = 0.2) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.2) +
    scale_color_manual(values = palette) +
    facet_rep_wrap(~type,
      ncol = 4, scales = "fixed",
      repeat.tick.labels = TRUE
    ) +
    scale_x_continuous(
      breaks = x_breaks, # limits = x_limits
      labels = x_labels
    ) +
    scale_y_continuous(breaks = sd_y_breaks, limits = sd_y_limits) +
    ggtitle("") +
    ylab("Conservation") +
    xlab("site") +
    # guides(fill=none) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 90, # size=5,
      face = "bold",
      vjust = 0.4, hjust = 1
    ))

  return(list(
    disp = direct.label(disp_p, method = list("right.polygons")),
    sd = direct.label(sd_p, method = list("right.polygons"))
  ))
}

xy_mer_disp_sd_plot_101 <- function(plotdata, split_by = "samples", disp_y_breaks = NULL, sd_y_breaks = NULL) {
  set.seed(123456789)
  palette <- randomColor(count = 60) # 随机生成60种颜色，其实里面有重复的
  palette <- distinctColorPalette(60)
  plotdata[["site"]] <- as.numeric(plotdata[["site"]])
  stream <- round(median(plotdata$site))
  if (split_by == "samples") {
    plotdata <- plotdata %>%
      mutate(type = substr(XYmer, 1, 2), site0 = site - stream) %>%
      mutate(site = if_else(site0 <= 0, site0 - 1, site0))
  } else if (split_by == "XYmer") {
    plotdata <- plotdata %>%
      mutate(site0 = site - stream) %>%
      mutate(site = if_else(site0 <= 0, site0 - 1, site0))
  } else {
    stop("Only support samples or XYmer!")
  }

  start <- -10
  stop <- 10

  # breaks <- c(-20,seq(start,stop,1),20)[-c((stop-start+2)/2+1)]
  x_breaks <- plotdata[["site0"]]
  if (is.null(disp_y_breaks)) {
    disp_y_breaks <- seq(from = 0, to = 2.4, by = 0.2)
  }
  disp_y_limits <- c(min(disp_y_breaks), max(disp_y_breaks))
  # x_limits <- c(min(plotdata[["site"]]),max(plotdata[["site"]]))
  # y_limits <- c(min(y_breaks), max(y_breaks))
  x_labels <- as.character(plotdata[["site"]])
  x_labels[grep("-", x_labels, invert = TRUE)] <- paste0("+", x_labels[grep("-", x_labels, invert = TRUE)])

  disp_p <- ggplot(data = plotdata, mapping = aes(x = site0, y = Separation, fill = ss_type)) +
    geom_point(aes(color = ss_type), show.legend = FALSE) +
    geom_line(aes(color = ss_type), show.legend = FALSE) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.2) +
    geom_vline(xintercept = 1, color = "black", linetype = "dashed", linewidth = 0.2) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.2) +
    scale_color_manual(values = palette) +
    facet_rep_wrap(~XYmer,
      ncol = 4, scales = "fixed",
      repeat.tick.labels = TRUE
    ) +
    scale_x_continuous(
      breaks = x_breaks, # limits = limits
      labels = x_labels
    ) +
    scale_y_continuous(breaks = disp_y_breaks, limits = disp_y_limits) +
    ggtitle("") +
    ylab("Separation") +
    xlab("site") +
    # guides(fill=none) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 90, # size=5,
      face = "bold",
      vjust = 0.4, hjust = 1
    ))
  if (is.null(sd_y_breaks)) {
    sd_y_breaks <- seq(from = 0, to = 2.4, by = 0.2)
  }
  sd_y_limits <- c(min(sd_y_breaks), max(sd_y_breaks))
  sd_p <- ggplot(data = plotdata, mapping = aes(x = site0, y = Conservation, fill = ss_type)) +
    geom_point(aes(color = ss_type), show.legend = FALSE) +
    geom_line(aes(color = ss_type), show.legend = FALSE) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.2) +
    geom_vline(xintercept = 1, color = "black", linetype = "dashed", linewidth = 0.2) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.2) +
    scale_color_manual(values = palette) +
    facet_rep_wrap(~XYmer,
      ncol = 4, scales = "fixed",
      repeat.tick.labels = TRUE
    ) +
    scale_x_continuous(
      breaks = x_breaks, # limits = x_limits
      labels = x_labels
    ) +
    scale_y_continuous(breaks = sd_y_breaks, limits = sd_y_limits) +
    ggtitle("") +
    ylab("Conservation") +
    xlab("site") +
    # guides(fill=none) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 90, # size=5,
      face = "bold",
      vjust = 0.4, hjust = 1
    ))

  return(list(
    disp = direct.label(disp_p, method = list("right.polygons")),
    sd = direct.label(sd_p, method = list("right.polygons"))
  ))
}

xy_mer_disp_sd_plot_50 <- function(plotdata, split_by = "samples", disp_y_breaks = NULL, sd_y_breaks = NULL, is_facet = TRUE) {
  set.seed(123456789)
  palette <- randomColor(count = 60) # 随机生成60种颜色，其实里面有重复的
  palette <- distinctColorPalette(60)

  plotdata[["site"]] <- as.numeric(plotdata[["site"]])
  # stream <- round(median(plotdata$site))
  stream <- 50
  if (split_by == "samples") {
    plotdata <- plotdata %>%
      mutate(type = substr(XYmer, 1, 2), site0 = site - stream) %>%
      mutate(site = if_else(site0 <= 0, site0 - 1, site0))
  } else if (split_by == "XYmer") {
    plotdata <- plotdata %>%
      mutate(site0 = site - stream) %>%
      mutate(site = if_else(site0 <= 0, site0 - 1, site0))
  } else {
    stop("Only support samples or XYmer!")
  }

  start <- -10
  stop <- 10

  # breaks <- c(-20,seq(start,stop,1),20)[-c((stop-start+2)/2+1)]
  x_breaks <- plotdata[["site0"]][c(
    seq(1, stream - 9, 10), 45:46, 50:51,
    seq(stream + 10, stream * 2, 10)
  )]
  if (is.null(disp_y_breaks)) {
    disp_y_breaks <- seq(from = 0, to = 2.4, by = 0.2)
  }
  disp_y_limits <- c(min(disp_y_breaks), max(disp_y_breaks))
  # x_limits <- c(min(plotdata[["site"]]),max(plotdata[["site"]]))
  # y_limits <- c(min(y_breaks), max(y_breaks))
  x_labels <- as.character(plotdata[["site"]][c(
    seq(1, stream - 9, 10), 45:46, 50:51,
    seq(stream + 10, stream * 2, 10)
  )])
  sel <- x_labels %in%
    as.character(c(
      -49:-41, -39:-31, -29:-21,
      # seq(-18,-2,2), seq(2,18,2),
      -19:-11, -9:-7, -6:-5, -4:-2, -1, 1,
      2:4, 5, 6, 7:9, 11:19,
      21:29, 31:39, 41:49
    ))
  x_labels[grep("-", x_labels, invert = TRUE)] <- paste0("+", x_labels[grep("-", x_labels, invert = TRUE)])
  x_labels[sel] <- ""

  point_x <- c(-5, -4, 0, 1)

  disp_p <- ggplot(data = plotdata, mapping = aes(x = site0, y = Separation, fill = ss_type)) +
    geom_line(aes(color = ss_type), show.legend = TRUE, linewidth = 0.2) +
    geom_point(
      mapping = aes(x = site0, y = Separation, color = ss_type), size = 0.2, # color="black",
      # data = (plotdata %>% filter(site0 %in% point_x)),
      show.legend = FALSE
    ) +
    geom_vline(xintercept = -5, color = "black", linetype = "dashed", linewidth = 0.2) +
    geom_vline(xintercept = -4, color = "black", linetype = "dashed", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "black", linetype = "dotted", linewidth = 0.2) +
    geom_vline(xintercept = 1, color = "black", linetype = "dotted", linewidth = 0.2) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.2) +
    # scale_color_manual(values = palette)
    scale_fill_brewer(palette = "Set1") +
    annotate( # 添加箭头
      geom = "curve", x = -19, y = max(disp_y_breaks) - 0.2, xend = -5,
      yend = max(disp_y_breaks) - 0.2,
      colour = "orange",
      curvature = 0, arrow = arrow(length = unit(2, "mm"))
    ) +
    annotate(
      geom = "text", x = -20, y = max(disp_y_breaks) - 0.2,
      label = "-6", hjust = "right"
    ) +
    annotate( # 添加箭头
      geom = "curve", x = -15, y = max(disp_y_breaks), xend = -4, yend = max(disp_y_breaks),
      colour = "orange",
      curvature = 0, arrow = arrow(length = unit(2, "mm"))
    ) +
    annotate(
      geom = "text", x = -16, y = max(disp_y_breaks),
      label = "-5", hjust = "right"
    ) +
    annotate( # 添加箭头
      geom = "curve", x = 10, y = max(disp_y_breaks),
      xend = 0, yend = max(disp_y_breaks),
      colour = "orange",
      curvature = .3, arrow = arrow(length = unit(2, "mm"))
    ) +
    annotate(
      geom = "text", x = 11, y = max(disp_y_breaks),
      label = "-1", hjust = "left"
    ) +
    annotate( # 添加箭头
      geom = "curve", x = 14, y = max(disp_y_breaks) - 0.2,
      xend = 1, yend = max(disp_y_breaks) - 0.2,
      colour = "orange",
      curvature = .3, arrow = arrow(length = unit(2, "mm"))
    ) +
    annotate(
      geom = "text", x = 15, y = max(disp_y_breaks) - 0.2,
      label = "+1", hjust = "left"
    )

  if (is.null(sd_y_breaks)) {
    sd_y_breaks <- seq(from = 0, to = 2.4, by = 0.2)
  }
  sd_y_limits <- c(min(sd_y_breaks), max(sd_y_breaks))
  sd_p <- ggplot(data = plotdata, mapping = aes(x = site0, y = Conservation, fill = ss_type)) +
    geom_line(aes(colour = ss_type), show.legend = FALSE, linewidth = 0.2) +
    geom_point(
      mapping = aes(x = site0, y = Conservation, color = ss_type), size = 0.2, # color="black",
      # data = (plotdata %>% filter(site0 %in% point_x)),
      show.legend = FALSE
    ) +
    geom_vline(xintercept = -5, color = "black", linetype = "dashed", linewidth = 0.2) +
    geom_vline(xintercept = -4, color = "black", linetype = "dashed", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "black", linetype = "dotted", linewidth = 0.2) +
    geom_vline(xintercept = 1, color = "black", linetype = "dotted", linewidth = 0.2) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.2) +
    # scale_color_manual(values = palette)
    scale_fill_brewer(palette = "Set1") +
    annotate( # 添加箭头
      geom = "curve", x = -19, y = max(sd_y_breaks) - 0.2,
      xend = -5, yend = max(sd_y_breaks) - 0.2,
      colour = "orange",
      curvature = 0, arrow = arrow(length = unit(2, "mm"))
    ) +
    annotate(
      geom = "text", x = -20, y = max(sd_y_breaks) - 0.2,
      label = "-6", hjust = "right"
    ) +
    annotate( # 添加箭头
      geom = "curve", x = -15, y = max(sd_y_breaks),
      xend = -4, yend = max(sd_y_breaks),
      colour = "orange",
      curvature = 0, arrow = arrow(length = unit(2, "mm"))
    ) +
    annotate(
      geom = "text", x = -16, y = max(sd_y_breaks),
      label = "-5", hjust = "right"
    ) +
    annotate( # 添加箭头
      geom = "curve", x = 10, y = max(sd_y_breaks),
      xend = 0, yend = max(sd_y_breaks),
      colour = "orange",
      curvature = .3, arrow = arrow(length = unit(2, "mm"))
    ) +
    annotate(
      geom = "text", x = 11, y = max(sd_y_breaks),
      label = "-1", hjust = "left"
    ) +
    annotate( # 添加箭头
      geom = "curve", x = 14, y = max(sd_y_breaks) - 0.2,
      xend = 1, yend = max(sd_y_breaks) - 0.2,
      colour = "orange",
      curvature = .3, arrow = arrow(length = unit(2, "mm"))
    ) +
    annotate(
      geom = "text", x = 15, y = max(sd_y_breaks) - 0.2,
      label = "+1", hjust = "left"
    )

  if (is_facet) {
    disp_p <- disp_p +
      # facet_wrap( ~ XYmer, ncol = 4) +
      facet_rep_wrap(~XYmer,
        ncol = 4, scales = "fixed",
        repeat.tick.labels = TRUE
      ) +
      scale_x_continuous(
        breaks = x_breaks, # limits = limits
        labels = x_labels
      ) +
      scale_y_continuous(breaks = disp_y_breaks, limits = disp_y_limits) +
      ggtitle("") +
      # ylab("Separation") +
      ylab("") +
      # xlab("site") +
      xlab("") +
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
      )
    sd_p <- sd_p +
      # facet_wrap( ~ XYmer, ncol = 4) +
      facet_rep_wrap(~XYmer,
        ncol = 4, scales = "fixed",
        repeat.tick.labels = TRUE
      ) +
      scale_x_continuous(
        breaks = x_breaks, # limits = x_limits
        labels = x_labels
      ) +
      scale_y_continuous(breaks = sd_y_breaks, limits = sd_y_limits) +
      ggtitle("") +
      # ylab("Conservation") +
      ylab("") +
      # xlab("site") +
      xlab("") +
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
      )
  } else {
    disp_p <- disp_p +
      scale_x_continuous(
        breaks = x_breaks, # limits = limits
        labels = x_labels
      ) +
      scale_y_continuous(breaks = disp_y_breaks, limits = disp_y_limits) +
      ggtitle("") +
      # ylab("Separation") +
      ylab("") +
      # xlab("site") +
      xlab("") +
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
      )
    sd_p <- sd_p +
      scale_x_continuous(
        breaks = x_breaks, # limits = x_limits
        labels = x_labels
      ) +
      scale_y_continuous(breaks = sd_y_breaks, limits = sd_y_limits) +
      ggtitle("") +
      # ylab("Conservation") +
      ylab("") +
      # xlab("site") +
      xlab("") +
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
      )
  }

  return(list(
    disp = direct.label(
      disp_p,
      method = list("smart.grid")
    ),
    sd = direct.label(
      sd_p,
      method = list("smart.grid")
    )
  ))

  # return(list(disp = disp_p +
  #               directlabels::geom_dl(aes(label = ss_type, colour = ss_type),  ,method = "smart.grid"),
  #             sd = sd_p +
  #               directlabels::geom_dl(aes(label = ss_type, colour = ss_type), method = "smart.grid")))
}

xy_mer_disp_sd_plot_batch <- function(plotdata_lst, FUN = "xy_mer_disp_sd_plot10", split_by = "samples") {
  mapply(FUN, plotdata = plotdata_lst, split_by = split_by, SIMPLIFY = FALSE)
  # lapply(plotdata_lst, FUN)
}
