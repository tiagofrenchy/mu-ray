# Utiliza a bigbigtable para filtrar somente os genes que possuem
# "secreted" em sua descrição.

library(readr)
library(dplyr)

# Vittorio
setwd("~/IC Alexandre")

# Tiago
# setwd("~/mu-ray_data/")

bigbigtable_only_fdr <- read_csv("Bigbigtable (FDR rejected and ABSOLUTE expression).csv")
A <- dplyr::filter(bigbigtable_only_fdr_log, grepl('secret', description))

bigbigtable_only_fdr_log <- read_csv("Bigbigtable (FDR rejected and LOG expression).csv")
B <- dplyr::filter(bigbigtable_only_fdr_log, grepl('secret', description))

secreted <- rbind(A, B)
secreted <- secreted[!duplicated(secreted),]

write.csv(secreted, "Genes with potentially secreted products.csv", row.names = FALSE)
