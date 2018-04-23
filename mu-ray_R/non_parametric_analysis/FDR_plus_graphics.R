library(readr)
library(sgof)
library(dplyr)
library(ggplot2)
library(stringr)


# Vittorio
#setwd("~/IC Alexandre")

# Tiago
setwd("~/mu-ray_data/")

# Abre os dados
clean_data <- read_csv("clean_data_without_duplicates.csv")
clean_data$expression <- log(clean_data$expression)

results <- read_csv("Resultado dos testes estatísticos NON-PARAMETRIC.csv")

# # Desenha gráfico com valores de p para
# # visualizacao dos resultados dos testes nao parametricos
# ggplot() +
#   geom_freqpoly(data = results, aes(x=p_SMACK), binwidth = NULL, color = "red", alpha = 1) +
#   geom_freqpoly(data = results, aes(x=p_ctrlSepsis), binwidth = NULL, color = "blue", alpha = 1) +
#   geom_freqpoly(data = results, aes(x=p_dun12), binwidth = NULL, color = "green", alpha = .5) +
#   geom_freqpoly(data = results, aes(x=p_dun13), binwidth = NULL, color = "green", alpha = .5) +
#   geom_freqpoly(data = results, aes(x=p_dun14), binwidth = NULL, color = "green", alpha = .5) +
#   geom_freqpoly(data = results, aes(x=p_dun23), binwidth = NULL, color = "green", alpha = .5) +
#   geom_freqpoly(data = results, aes(x=p_dun24), binwidth = NULL, color = "green", alpha = .5) +
#   geom_freqpoly(data = results, aes(x=p_dun34), binwidth = NULL, color = "green", alpha = .5) +
#   geom_freqpoly(data = results, aes(x=p_ctrlDia1), binwidth = NULL, color = "yellow", alpha = .5) +
#   geom_freqpoly(data = results, aes(x=p_ctrlDia2), binwidth = NULL, color = "yellow", alpha = .5) +
#   geom_freqpoly(data = results, aes(x=p_ctrlDia3), binwidth = NULL, color = "yellow", alpha = .5) +
#   geom_freqpoly(data = results, aes(x=p_ctrlDia4), binwidth = NULL, color = "yellow", alpha= .5) +
#   labs(title = "Curva de Frequencia dos p-values", x = "p-value", y = "count")

  
  
# Aplica fdr a comparacao entre os genes da sepse
fdr_smack  <- BH(results$p_SMACK, alpha = 0.001)
results <- results[order(results$p_SMACK),]
id_fdr <- head(results$transcript_cluster_id, n = fdr_smack$Rejections)

# Aplica fdr a comparacao entre controle e sepse
fdr_ctrl   <- BH(results$p_ctrlSepsis, alpha = 0.001)
results <- results[order(results$p_ctrlSepsis),]
id_ctrl <- head(results$transcript_cluster_id, n = fdr_ctrl$Rejections)

# Encontra e seleciona todos os genes que sao significantes segundo os dois fdr precedentes
coolGenesId <- intersect(id_fdr, id_ctrl)
coolGenes <- clean_data %>% filter(transcript_cluster_id %in% coolGenesId)

# Loop que cria os plots para todos os genes selecionados precedentemente
listUniqueTranscripts <- unique(coolGenes$transcript_cluster_id)
numberTranscripts <- length(listUniqueTranscripts)

for (i in seq(numberTranscripts)) {
  reallyCoolGenes <- coolGenes %>% filter(transcript_cluster_id == listUniqueTranscripts[i])
  
  nome <- reallyCoolGenes$description[1]
  if (nchar(nomeTit)<90){
    nomeTit <- paste(substr(nome, 1, 44), "\n", substr(nome, 45, nchar(nome)), sep = "")
  }
  if (nchar(nomeTit)>=90){
    nomeTit <- paste(substr(nome, 1, 44), "\n", substr(nome, 45, 89),
                     "\n", substr(nome, 90, nchar(nome)), sep = "")
  }
  
  id <- factor(listUniqueTranscripts[i])
  reallyCoolGenes$Dia <- gsub("D0", "Controle", reallyCoolGenes$Dia)
  
  todos <- subset(reallyCoolGenes, Dia != 'D0')
  todos$Dia <- gsub("D1", "Sepse", todos$Dia)
  todos$Dia <- gsub("D2", "Sepse", todos$Dia)
  todos$Dia <- gsub("D3", "Sepse", todos$Dia)
  todos$Dia <- gsub("D4", "Sepse", todos$Dia)
  
  
  ggplot() +
    geom_boxplot(data = reallyCoolGenes, aes(x=Dia, y=expression, fill=Dia), alpha=0.5) +
    geom_boxplot(data = todos, aes(x=Dia, y=expression, fill=Dia), alpha=0.5) +
    theme(legend.position="none") +
    ggtitle(nomeTit) +
    scale_fill_brewer(palette = "Set1") +
    ggsave(sprintf("[%s] - %s.png",id, nome), path = paste(getwd(), "/boxplots", sep=""))
  
}


##############################################
##############################################

sepsisExpression <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i], Dia != 'D0') %>% 
  select(Dia, expression, Paciente)


ID <- 7974851 

# Gráfico
reallyCoolGenes <- clean_data %>% filter(transcript_cluster_id == ID)
ggplot(reallyCoolGenes, aes(x=Dia, y=expression, fill=Dia)) +
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Dark2")

 
d <- clean_data %>% filter(transcript_cluster_id == ID) %>% 
  select(Dia, expression, Paciente)

d <- clean_data %>% filter(transcript_cluster_id == ID, Dia == 'D1') %>% 
  select(Dia, expression, Paciente)

d2 <- clean_data %>% filter(transcript_cluster_id == ID, Dia == 'D4') %>% 
  select(Dia, expression, Paciente)