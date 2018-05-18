library(readr)
library(dplyr)
library(tidyr)
library(poibin)
library(sgof)

inicio <- Sys.time()

# Vittorio
# setwd("~/IC Alexandre")

# Tiago
setwd("~/mu-ray_data/")

# Abre os dados
clean_data <- read_csv("clean_data_without_duplicates.csv")

# Iniciamos o data frame da regulacao
regulation <- data.frame(transcript_cluster_id = double(), 
                         description = character(),
                         media_controle = double(),
                         sd_controle = double(),
                         media_sepse = double(),
                         sd_sepse = double(),
                         wilcoxon = double(),
                         diferenca_medias = double(),
                         razao_medias = double(),
                         media_log_controle = double(),
                         sd_log_controle = double(),
                         media_log_sepse = double(),
                         sd_log_sepse = double(),
                         wilcoxon_log = double(),
                         diferenca_log_medias = double(),
                         razao_log_medias = double(),
                         stringsAsFactors = FALSE)

# Definimos os TCID para serem usados no loop
listUniqueTranscripts <- unique(clean_data$transcript_cluster_id)
numberTranscripts     <- length(listUniqueTranscripts)

for (i in seq(numberTranscripts)) {
  
  # Para o gene atual, selecionamos a sua expressão em todos os dias do grupo com sepse e controle.
  control <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                   Dia == 'D0') %>% select(Dia, expression, Paciente)
  sepsis  <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                  Dia == 'D1') %>% select(Dia, expression, Paciente)
  
  # criando coluna com o nome do dataframe para fazer merge
  control["expression_controle"] <- control$expression
  sepsis["expression_sepsis"]    <- sepsis$expression
  
  # Calculamos o log de todas as expressões.
  control["log_expression_controle"] <- log(control$expression)
  sepsis["log_expression_sepsis"]    <- log(sepsis$expression)
  
  # Calculamos o wilcoxon-mann-whitney entre as medias de expressão
  wilcoxon     <- wilcox.test(control$expression_controle, sepsis$expression_sepsis)
  wilcoxon_log <- wilcox.test(control$log_expression_controle, sepsis$log_expression_sepsis)
  
  # Calculamos a media de todas as expressões.
  control_avg <- mean(control$expression_controle)
  sepsis_avg  <- mean(sepsis$expression_sepsis)
  
  # Calculamos o desvio padrao das medias.
  sd_ctrl   <- sd(control$expression_controle)
  sd_sepsis <- sd(sepsis$expression_sepsis)
  
  # Calculamos a media dos logs de todas as expressões.
  control_log_avg <- mean(control$log_expression_controle)
  sepsis_log_avg  <- mean(sepsis$log_expression_sepsis)
  
  # Calculamos o desvio padrao das medias dos logs.
  sd_log_ctrl   <- sd(control$log_expression_controle)
  sd_log_sepsis <- sd(sepsis$log_expression_sepsis)
  
  # Calculamos as diferenças entre as médias.
  diferenca_medias     <- (sepsis_avg - control_avg)
  diferenca_log_medias <- (sepsis_log_avg - control_log_avg)
  
  # Calculamos as razões entre as médias.
  razao_medias     <- (sepsis_avg / control_avg)
  razao_log_medias <- (sepsis_log_avg / control_log_avg)
  
  # Obtemos a descrição do gene.
  descr <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i]) %>% 
    select(description)
  
  # Agora, inserimos todos os valores calculados no data frame.
  regulation[nrow(regulation) + 1,]     = c(listUniqueTranscripts[i],
                                      c(descr[1][1,]),
                                      control_avg,
                                      sd_ctrl,
                                      sepsis_avg,
                                      sd_sepsis,
                                      wilcoxon$p.value,
                                      diferenca_medias,
                                      razao_medias,
                                      control_log_avg,
                                      sd_log_ctrl,
                                      sepsis_log_avg,
                                      sd_log_sepsis,
                                      wilcoxon_log$p.value,
                                      diferenca_log_medias,
                                      razao_log_medias)
  
  # DEBUG
  print(sprintf('%5.1f%%', i/length(listUniqueTranscripts)*100), quote = FALSE)
  
}

########## FDR ##########

# Aplica FDR à comparação entre o controle e o D1 da sepse (i.e., valores de p derivados do Wilcoxon).
# Ordena os resultados em função do p value: os genes com os menores ficaram no topo da tabela.
# Então, selecionamos as n rejeições do topo da tabela.

# A partir das rejeições (i.e., diferenças aceitas pelo FDR), filtramos os nossos resultados para
# que possamos visualizar só estes genes

fdr <- BH(regulation$wilcoxon, alpha = .001)
regulation <- regulation[order(regulation$wilcoxon),]
id_fdr <- head(regulation$transcript_cluster_id, n = fdr$Rejections)
results_fdr <- regulation %>% filter(transcript_cluster_id %in% id_fdr)

fdr_log <- BH(regulation$wilcoxon_log, alpha = .001)
regulation <- regulation[order(regulation$wilcoxon_log),]
id_fdr_log <- head(regulation$transcript_cluster_id, n = fdr_log$Rejections)
results_fdr_log <- regulation %>% filter(transcript_cluster_id %in% id_fdr_log)

write.csv(regulation, "Regulação dos genes.csv", row.names = FALSE)

write.csv(results_fdr, "Regulação dos genes - fdr valores absolutos.csv", row.names = FALSE)

write.csv(results_fdr_log, "Regulação dos genes - fdr logs.csv", row.names = FALSE)

fim <- Sys.time()

print(fim-inicio)

print("fuckYEAHHH ... #VaiPatente", quote = FALSE)
