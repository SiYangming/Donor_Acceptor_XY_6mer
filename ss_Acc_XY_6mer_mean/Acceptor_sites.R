##### save to file: three splicing mode #####
load("../data/DAssBackup.rda")
library(tidyverse)
con_donor_num <- rt_acc %>%
  dplyr::distinct(chr, con_up, con_down,
    score, strand,
    .keep_all = TRUE
  ) %>%
  nrow()
con_donor_num
# [1] 55557

con_acc_num <- rt_donor %>%
  dplyr::distinct(chr, con_up, con_down,
    score, strand,
    .keep_all = TRUE
  ) %>%
  nrow()
con_acc_num
# [1] 48414

nrows <- min(nrow(rt_norm), con_acc_num, nrow(rt_acc))
nrows
# [1] 48414

##### Acceptor Splicing mode #####
### Constitutive Acceptor
if (!file.exists("data/Constitutive_Acceptor.bed")) {
  file.copy(
    "../bed/Constitutive_Acceptor.bed",
    "data/Constitutive_Acceptor.bed"
  )
}

### Common Acceptor
if (!file.exists("data/Common_Acceptor.bed")) {
  read_tsv("../bed/Common_Acceptor.bed", col_names = F) %>%
    random_select(num = nrows, by = "row") %>%
    write_tsv(
      file = "data/Common_Acceptor.bed",
      col_names = F
    )
}


### Alternative Acceptor
if (!file.exists("data/Alternative_Acceptor.bed")) {
  read_tsv("../bed/Alternative_Acceptor.bed", col_names = F) %>%
    random_select(num = nrows, by = "row") %>%
    write_tsv(
      file = "data/Alternative_Acceptor.bed",
      col_names = F
    )
}

acc_num <- rt_acc %>%
  filter(ss_type != "unknow") %>%
  select(ss_type) %>%
  table()
acc_num
# canonical    exonic  intronic
#     12356     50497    144098
nrows <- min(acc_num)
nrows
# 12356
### Constitutive Acceptor
# read_tsv("bed/Constitutive_Acceptor.bed", col_names = F) %>%
#   random_select(num = nrows, by = "row") %>%
#   mutate(label = "Constitutive Acceptor") %>%
#   write_tsv(file = "data/Constitive_Acceptor.bed",
#             col_names = F)

### Common Acceptor
# read_tsv("bed/Common_Acceptor.bed", col_names = F) %>%
#   random_select( num = nrows, by = "row") %>%
#   mutate(label = "Common Donor") %>%
#   write_tsv(file = "data/Common_Acceptor.bed",
#             col_names = F)

### Normal Acceptor
rt_acc %>%
  filter(ss_type == "canonical") %>%
  random_select(num = nrows, by = "row") %>%
  select(chr, alt_up, alt_down, score, ss_type, strand) %>%
  mutate(label = "Normal Acceptor") %>%
  distinct(.keep_all = TRUE) %>%
  write_tsv(
    file = "data/Normal_Acceptor.bed",
    col_names = F
  )

### Intronic Acceptor
rt_acc %>%
  filter(ss_type == "intronic") %>%
  random_select(num = nrows, by = "row") %>%
  select(chr, alt_up, alt_down, score, ss_type, strand) %>%
  mutate(label = "Intronic Acceptor") %>%
  write_tsv(intronic_sel,
    file = "data/Intronic_Acceptor.bed",
    col_names = F
  )

### Exonic Acceptor
rt_acc %>%
  filter(ss_type == "exonic") %>%
  random_select(num = nrows, by = "row") %>%
  select(chr, alt_up, alt_down, score, ss_type, strand) %>%
  mutate(label = "Exonic Acceptor") %>%
  write_tsv(
    file = "data/Exonic_Acceptor.bed",
    col_names = F
  )
