# geom = rep("curve", 4)
# x = c(-19, -15, 10, 14)
# y = c(2.4, 2.6, 2.6, 2.4)
# xend = c(-5, -4, 0, 1)
# yend = c(2.4, 2.6, 2.6, 2.4)
# curvature = c(0, 0, .3, .3)
# labels = c("-6", "-5", "-1", "+1")
# colour = rep("orange", 4)
# units <- list(c(2,"mm"), c(2,"mm"), c(2,"mm"), c(2,"mm"))

annotate_sites <- function(geom, x = NULL, y = NULL, xmin = NULL, xmax = NULL,
                           ymin = NULL, ymax = NULL, xend = NULL, yend = NULL,
                           colour = NULL, curvature = NULL, units = NULL,
                           label = NULL, ...,
                           na.rm = FALSE) {
  cnt <- 1
  layer_lst <- list()
  for (i in 1:length(x)) {
    layer_lst[[cnt]] <- annotate(
      geom = geom[i],
      x = x[i], y = y[i],
      xend = xend[i], yend = yend[i],
      colour = colour[i], curvature = curvature[i],
      arrow = arrow(length = unit(units[[i]][1], units[[i]][2]))
    )
    x1 <- ifelse(x <= 0, x - 1, x + 1)
    hjusts <- ifelse(x <= 0, "right", "left")
    layer_lst[[cnt + 1]] <- annotate(
      geom = "text", x = x1[i], y = y[i],
      label = label[i], hjust = hjusts[i]
    )
    cnt <- cnt + 2
  }
  return(layer_lst)
}

# v = c(-5, -4, 0, 1)
# h = c(1)
# color = rep("black", 4)
# vlinetype = rep(c("dashed", "dotted"), 2)
# hlinetype = c("dashed")
# vlinewidth = rep(0.2 ,4)
# hlinewidth = c(0.2)
geom_vhline <- function(mapping = NULL, data = NULL, ...,
                        v, h,
                        vlinetype, hlinetype,
                        vlinewidth, hlinewidth,
                        color,
                        na.rm = FALSE,
                        show.legend = NA) {
  v_cnt <- 1
  h_cnt <- 1
  vlayer_lst <- list()
  hlayer_lst <- list()
  for (i in 1:length(v)) {
    vlayer_lst[[v_cnt]] <- geom_vline(
      xintercept = v[i], color = color[i],
      linetype = vlinetype[i], linewidth = vlinewidth[i]
    )
    v_cnt <- v_cnt + 1
  }
  for (i in 1:length(h)) {
    hlayer_lst[[h_cnt]] <- geom_hline(
      yintercept = h[i], color = color[i],
      linetype = hlinetype[i], linewidth = hlinewidth[i]
    )
    h_cnt <- h_cnt + 1
  }
  return(c(vlayer_lst, hlayer_lst))
}
