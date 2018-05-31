# Utiliza a bigbigtable para filtrar somente os genes que possuem
# "secreted" em sua descrição.

library(readr)
library(dplyr)

# Vittorio
setwd("~/IC Alexandre")

# Tiago
# setwd("~/mu-ray_data/")

bigbigtable_only_fdr_log <- read_csv("Bigbigtable (FDR rejected and LOG expression.csv")
secreted <- dplyr::filter(bigbigtable_only_fdr_log, grepl('secret', description))


write.csv(secreted, "Genes with potentially secreted products (LOG expression).csv", row.names = FALSE)
