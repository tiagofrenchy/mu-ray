library(readr)
library(dplyr)
library(tidyr)
library(poibin)
library(sgof)

# Vittorio
setwd("~/IC Alexandre")

# Tiago
# setwd("~/mu-ray_data/")

# Abre os dados
clean_data <- read_csv("clean_data_without_duplicates.csv")

# Iniciamos o data frame da regulacao
bigbigtable <- data.frame(transcript_cluster_id = double(), 
                         description = character(),
                         reference = character(),
                         control_avg = double(),
                         control_sd = double(),
                         sepsis_avg = double(),
                         sepsis_sd = double(),
                         wilcoxon = double(),
                         difference_avg = double(),
                         ratio_avg = double(),
                         control_avg_log = double(),
                         control_sd_log = double(),
                         sepsis_avg_log = double(),
                         sepsis_sd_log= double(),
                         wilcoxon_log = double(),
                         difference_avg_log= double(),
                         ratio_avg_log = double(),
                         stringsAsFactors = FALSE)

# Definimos os TCID para serem usados no loop
listUniqueTranscripts <- unique(clean_data$transcript_cluster_id)
numberTranscripts     <- length(listUniqueTranscripts)

#DEBUG
inicio <- Sys.time()

for (i in seq(numberTranscripts)) {
  
  # Para o gene atual, selecionamos a sua expressão em todos os dias do grupo com sepse e controle.
  control <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                   Dia == 'D0') %>% select(Dia, expression, Paciente)
  sepsis  <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                  Dia == 'D1') %>% select(Dia, expression, Paciente)
  
  # Calculamos o log de todas as expressões.
  control["log_expression"] <- log(control$expression)
  sepsis["log_expression"]  <- log(sepsis$expression)
  
  # Calculamos o wilcoxon-mann-whitney entre as medias de expressão
  wilcoxon     <- wilcox.test(control$expression, sepsis$expression, paired = FALSE)
  wilcoxon_log <- wilcox.test(control$log_expression, sepsis$log_expression, paired = FALSE)
  
  # Calculamos a media de todas as expressões.
  control_avg <- mean(control$expression)
  sepsis_avg  <- mean(sepsis$expression)
  
  # Calculamos o desvio padrao das medias.
  control_sd <- sd(control$expression)
  sepsis_sd  <- sd(sepsis$expression)
  
  # Calculamos a media dos logs de todas as expressões.
  control_avg_log <- mean(control$log_expression)
  sepsis_avg_log  <- mean(sepsis$log_expression)
  
  # Calculamos o desvio padrao das medias dos logs.
  control_sd_log   <- sd(control$log_expression)
  sepsis_sd_log <- sd(sepsis$log_expression)
  
  # Calculamos as diferenças entre as médias.
  difference_avg     <- (sepsis_avg - control_avg)
  difference_avg_log <- (sepsis_avg_log - control_avg_log)
  
  # Calculamos as razões entre as médias.
  ratio_avg     <- (sepsis_avg     / control_avg)
  ratio_avg_log <- (sepsis_avg_log / control_avg_log)
  
  # Obtemos a descrição e a referência do gene do gene.
  descr <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i]) %>% 
    select(description)
  
  ref <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i]) %>% 
    select(gene_ref)
  
  # Agora, inserimos todos os valores calculados no data frame.
  bigbigtable[nrow(bigbigtable) + 1,] = c(listUniqueTranscripts[i],
                                      c(descr[1][1,]),
                                      c(ref[1][1,]),
                                      control_avg,
                                      control_sd,
                                      sepsis_avg,
                                      sepsis_sd,
                                      wilcoxon$p.value,
                                      difference_avg,
                                      ratio_avg,
                                      control_avg_log,
                                      control_sd_log,
                                      sepsis_avg_log,
                                      sepsis_sd_log,
                                      wilcoxon_log$p.value,
                                      difference_avg_log,
                                      ratio_avg_log)
  
  # DEBUG
  print(sprintf('%5.1f%%', i/length(listUniqueTranscripts)*100), quote = FALSE)
  
}

########## FDR ##########

# Aplica FDR à comparação entre o controle e o D1 da sepse (i.e., valores de p derivados do Wilcoxon).
# Ordena os resultados em função do p value: os genes com os menores ficaram no topo da tabela.
# Então, selecionamos as n rejeições do topo da tabela.

# A partir das rejeições (i.e., diferenças aceitas pelo FDR), filtramos os nossos resultados para
# que possamos visualizar só estes genes

## Primeiro, sem logaritmo.
# Ordena os valores de p do menor para o maior
bigbigtable <- bigbigtable[order(bigbigtable$wilcoxon),]

# Calcula o FDR, seleciona os N genes rejeitados (pegando o head, i.e., os menores valores de P)
# pelo transcript_id e filtra.
fdr <- BH(bigbigtable$wilcoxon, alpha = .001)
id_fdr <- head(bigbigtable$transcript_cluster_id, n = fdr$Rejections)
bigbigtable_only_fdr <- bigbigtable %>% filter(transcript_cluster_id %in% id_fdr)


## Agora, com logaritmo.
bigbigtable <- bigbigtable[order(bigbigtable$wilcoxon_log),]
fdr_log <- BH(bigbigtable$wilcoxon_log, alpha = .001)
id_fdr_log <- head(bigbigtable$transcript_cluster_id, n = fdr_log$Rejections)
bigbigtable_only_fdr_log <- bigbigtable %>% filter(transcript_cluster_id %in% id_fdr_log)


write.csv(bigbigtable, "Bigbigtable.csv", row.names = FALSE)
write.csv(bigbigtable_only_fdr, "Bigbigtable (FDR rejected and ABSOLUTE expression).csv", row.names = FALSE)
write.csv(bigbigtable_only_fdr_log, "Bigbigtable (FDR rejected and LOG expression).csv", row.names = FALSE)


# DEBUG
fim <- Sys.time()
print(fim-inicio)
print("fuckYEAHHH ... #VaiPatente", quote = FALSE)
