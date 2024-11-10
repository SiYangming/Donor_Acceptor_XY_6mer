#' random select data frame
#'
#' @param df A data frame
#' @param num How many to select
#' @param by "row" or "col"
#' @return A tibble
#' @author Yangming Si
#' @export
random_select <- function(df, num, by = "row") {
  if (by == "row") {
    random_ids <- sample(1:nrow(df), num)
    sel <- df[random_ids, ]
  } else if (by == "col") {
    random_ids <- sample(1:ncol(df), num)
    sel <- df[, random_ids]
  } else {
    print("Please select by = \"col\" for column, by = \"row\" for row.")
  }
  return(sel)
}

#' Calculate the distance between adjacent sites or between canonical sites and other sites
#'
#' @param df A data frame
#' @return A tibble
#' @author Yangming Si
#' @export
calc_distance <- function(df) {
  chr <- unique(df[["chr"]])
  strand <- unique(df[["strand"]])
  pos <- unique(df[["pos"]])
  annotated <- unique(df[["annotated"]])

  if (strand == "+" & pos == "acceptor" |
    strand == "-" & pos == "donor") {
    up <- unique(df[["up"]])
    print(paste("-----Procces", pos, chr, up, strand, annotated, "-----", sep = " "))
    # if (annotated == "One") {
    #   ref_site_df <- df %>% filter(source == "annotated") %>% select(down)
    #   ref_site <- ref_site_df[[1]]
    # }
    #
    # if (annotated == "More") {
    #   ref_site_df <- df %>% filter(source == "annotated") %>% select(down)
    #   ref_site <- min(ref_site_df[[1]])
    # }

    if (sum(df[["ss_type"]] == "canonical") == 1) {
      ref_site_df <- df %>%
        filter(ss_type == "canonical") %>%
        select(down)
      ref_site <- ref_site_df[[1]]
    } else {
      ref_site <- min(df[["down"]])
    }

    # return(df %>% mutate(distance = (down - ref_site)))
    return(df %>% mutate(
      distance = (down - ref_site),
      distance2 = abs(c(df$down[-1], 0) - df$down),
      distance3 = down - up
    ))
  }
  if (strand == "-" & pos == "acceptor" |
    strand == "+" & pos == "donor") {
    down <- unique(df[["down"]])
    print(paste("-----Procces", pos, chr, down, strand, annotated, "-----", sep = " "))
    # if (annotated == "One") {
    #   ref_site_df <- df %>% filter(source == "annotated") %>% select(up)
    #   ref_site <- ref_site_df[[1]]
    # }
    #
    # if (annotated == "More") {
    #   ref_site_df <- df %>% filter(source == "annotated") %>% select(up)
    #   ref_site <- min(ref_site_df[[1]])
    # }

    if (sum(df[["ss_type"]] == "canonical") == 1) {
      ref_site_df <- df %>%
        filter(ss_type == "canonical") %>%
        select(up)
      ref_site <- ref_site_df[[1]]
    } else {
      ref_site <- min(df[["up"]])
    }

    # return(df %>% mutate(distance = (up - ref_site)))

    return(df %>% mutate(
      distance = (up - ref_site),
      distance2 = abs(c(df$up[-1], 0) - df$up),
      distance3 = down - up
    ))
  }
}

#' Batch calculate the distance between canonical sites and other sites
#'
#' @param df A data frame
#' @return A tibble
#' @author Yangming Si
#' @export
calc_distance_batch <- function(df) {
  strand <- unique(df[["strand"]])
  pos <- unique(df[["pos"]])
  if (strand == "+" & pos == "acceptor" |
    strand == "-" & pos == "donor") {
    lst <- df %>%
      group_by(up, .drop = FALSE) %>%
      group_map(~ calc_distance(.x), .keep = TRUE)
    return(do.call(rbind, lst))
  }
  if (strand == "-" & pos == "acceptor" |
    strand == "+" & pos == "donor") {
    lst <- df %>%
      group_by(down, .drop = FALSE) %>%
      group_map(~ calc_distance(.x), .keep = TRUE)
    return(do.call(rbind, lst))
  }
}

#' Calculate the distance between adjacent sites
#'
#' @param df A data frame
#' @return A tibble
#' @author Yangming Si
#' @export
# calc_distance2 <- function(df){
#   chr <- unique(df[["chr"]])
#   strand <- unique(df[["strand"]])
#   pos <- unique(df[["pos"]])
#   annotated <- unique(df[["annotated"]])
#   print(paste( "-----Procces", pos, chr, strand, annotated,  "-----", sep = " "))
#   if (strand == "+" & pos == "acceptor" |
#       strand == "-" & pos == "donor") {
#     return(df %>% mutate(distance = df$down - df$up,
#                          distance2 = abs(c(df$down[-1],0) - df$down)))
#   }
#   if (strand == "-" & pos == "acceptor" |
#       strand == "+" & pos == "donor") {
#     return(df %>% mutate(distance = df$down - df$up,
#                          distance2 = abs(c(df$up[-1],0) - df$up)))
#   }
# }

#' Batch calculate the distance between adjacent sites
#'
#' @param df A data frame
#' @return A tibble
#' @author Yangming Si
#' @export
# calc_distance2_batch <- function(df){
#   strand <- unique(df[["strand"]])
#   pos <- unique(df[["pos"]])
#   if (strand == "+" & pos == "acceptor" |
#       strand == "-" & pos == "donor") {
#     lst <- df %>% group_by(up, .drop = FALSE) %>%
#       group_map(~calc_distance2(.x), .keep = TRUE)
#     # print()
#     return(do.call(rbind, lst))
#   }
#   if (strand == "-" & pos == "acceptor" |
#       strand == "+" & pos == "donor") {
#     lst <- df %>% group_by(down, .drop = FALSE) %>%
#       group_map(~calc_distance2(.x), .keep = TRUE)
#     return(do.call(rbind, lst))
#   }
# }

#' recognize sites type
#'
#' @param df A data frame
#' @param pos acceptor or donor
#' @return A tibble
#' @author Yangming Si
#' @export
count_con_sites <- function(df, pos = "acc") {
  chr <- unique(df[["chr"]])
  annotated <- unique(df[["annotated"]])
  if (pos == "acc") {
    minus <- df %>%
      filter(strand == "-") %>%
      group_by(down) %>%
      dplyr::summarise(n = n()) %>%
      nrow()
    plus <- df %>%
      filter(strand == "+") %>%
      group_by(up) %>%
      dplyr::summarise(n = n()) %>%
      nrow()
    return(tibble(
      chr = chr, strand = c("-", "+"),
      n = c(minus, plus),
      annotated = annotated
    ))
  } else if (pos == "donor") {
    minus <- df %>%
      filter(strand == "-") %>%
      group_by(up) %>%
      dplyr::summarise(n = n()) %>%
      nrow()
    plus <- df %>%
      filter(strand == "+") %>%
      group_by(down) %>%
      dplyr::summarise(n = n()) %>%
      nrow()
    return(tibble(
      chr = chr, strand = c("-", "+"),
      n = c(minus, plus),
      annotated = annotated
    ))
  }
}

#' PValue format
#'
#' @param pValue a vector of p value
#' @return formatted p value
#' @author Yangming si
#' @examples
#' library(qvalue)
#' data(hedenfalk)
#' p <- hedenfalk$p[1]
#' pValueFormat(p)
#' @export

pValueFormat <- function(pValue) {
  pval <- 0
  if (pValue > 0.05) {
    pval <- round(as.numeric(pValue), 3)
  }
  if (pValue < 0.05) {
    pval <- signif(as.numeric(pValue), 4)
    pval <- format(pval, scientific = TRUE)
  }
  return(pval)
}
