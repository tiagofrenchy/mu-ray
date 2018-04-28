library(readr)
library(dplyr)
library(sgof)
library(ggplot2)

# Vittorio
setwd("~/IC Alexandre")

# Tiago
#setwd("~/mu-ray_data/")

# Abre os dados
clean_data <- read_csv("clean_data_without_duplicates.csv")
clean_data$expression <- log(clean_data$expression)

# Iniciamos o data frame dos resultados
wilcoxon_results <- data.frame(transcript_cluster_id = double(), 
                               description   = character(),
                               p_ctrlD1 = double(),
                               stringsAsFactors = FALSE)


# Seleciona o número de genes (transcritos) distintos dentre os nossos dados em clean_data
listUniqueTranscripts <- unique(clean_data$transcript_cluster_id)
numberTranscripts <- length(listUniqueTranscripts)

# Para cada um desses genes, comparamos sua expressão no controle com o dia 1 da sepse.
for (i in seq(numberTranscripts)) {
  # Para o gene atual, selecionamos a sua expressão no dia D1.
  sepsisD1 <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                            Dia == 'D1') %>% select(expression)
  
  # Para o gene atual, selecionamos a sua expressão nos controles
  control <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                             Dia == 'D0') %>% select(expression)
  
  # Nota: aqui os testes NÃO são pareados.
  controleDia1 = wilcox.test(control$expression, 
                             sepsisD1$expression, paired = FALSE,
                             exact = TRUE)

  # Primeiro, obtemos a descrição do gene.
  descr <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i]) %>% 
    select(description)
  
  # Agora, inserimos todos os valores calculados no data frame.
  wilcoxon_results[nrow(wilcoxon_results) + 1,]     = c(listUniqueTranscripts[i],
                                                        c(descr[1][1,]),
                                                        controleDia1$p.value)
  
  # DEBUG
  print(sprintf('%5.1f%%', i/length(listUniqueTranscripts)*100), quote = FALSE)
}


# Aplica FDR à comparação entre o controle e o D1 da sepse (i.e., valores de p derivados do Wilcoxon).
# Ordena os resultados em função do p value: os genes com os menores ficaram no topo da tabela.
# Então, selecionamos as n rejeições do topo da tabela.
fdr  <- BH(wilcoxon_results$p_ctrlD1, alpha = 0.001)
wilcoxon_results <- wilcoxon_results[order(wilcoxon_results$p_ctrlD1),]
id_fdr <- head(wilcoxon_results$transcript_cluster_id, n = fdr$Rejections)


# A partir das rejeições (i.e., diferenças aceitas pelo FDR), filtramos os nossos resultados para
# que possamos visualizar só estes genes
results_fdr <- wilcoxon_results %>% filter(transcript_cluster_id %in% id_fdr)
listUniqueTranscripts <- unique(results_fdr$transcript_cluster_id)
numberTranscripts <- length(listUniqueTranscripts)

## Criamos coluan que apontará TRUE se o gene for up_regualted, falso se não.
results_fdr["up_regulated"] <- FALSE

# Agora, fazemos um loop com o intuito de, para cada gene com diferença significativa, ver se ele
# é up regulated ou down regualted.
for (i in seq(numberTranscripts)) {
  # Para o gene atual, selecionamos a sua expressão no dia D1.
  sepsisD1 <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                    Dia == 'D1') %>% select(expression)
  
  # Para o gene atual, selecionamos a sua expressão nos controles
  control <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                   Dia == 'D0') %>% select(expression)
  
  averageControl  = mean(control$expression)
  averageSepsisD1 = mean(sepsisD1$expression)
  if (averageControl < averageSepsisD1)
    results_fdr[i,]$up_regulated = TRUE
  
  print(i)
}



####
# MÉTODOS GRÁFICOS
####
#
# Encontra e seleciona todos os genes que sao significantes segundo os dois fdr precedentes.
# Importamos o dado do clean_data porque precisamos da expressão em cada um dos pacientes para montar
# os boxplots (que é o que o loop faz).
coolGenes <- clean_data %>% filter(transcript_cluster_id %in% id_fdr)

listUniqueTranscripts <- unique(coolGenes$transcript_cluster_id)
numberTranscripts <- length(listUniqueTranscripts)

for (i in seq(numberTranscripts)) {
  reallyCoolGenes <- coolGenes %>% filter(transcript_cluster_id == listUniqueTranscripts[i], (Dia == 'D1' | Dia == 'D0'))

  # Nos dados, reescrevemos "D0" como "Controle" para que apareça assim no gráfico. É equivalente -- só faz mais sentido.
  reallyCoolGenes$Dia <- gsub("D0", "Controle", reallyCoolGenes$Dia)
  
  # Retiramos potencial slash ("/") do nome do gene para não termos erro na hora de salvar o arquivo.
  reallyCoolGenes$description <- gsub("/", "-", reallyCoolGenes$description, fixed = TRUE)
  
  # Salvamos nome (descrição) e transcript_cluster_id para utilizar na criação do gráfico.
  nome <- reallyCoolGenes$description[i]
  id <- factor(listUniqueTranscripts[i])

  tamTit <- 7
  tamResto <- 8
  ggplot() +
    geom_boxplot(data = reallyCoolGenes, aes(x=Dia, y=expression, fill=Dia), alpha=0.5) +
    xlab("Grupos") +
    ylab("Expressão") +
    ggtitle(nome) +
    theme(legend.position="none",
          title = element_text(size = tamTit),
          axis.text = element_text(size = tamResto),
          axis.title = element_text(size = tamResto)) +
    scale_fill_brewer(palette = "Set1") +
    ggsave(sprintf("[%s] - %s.png",id, nome), path = paste(getwd(), "/diagnostic_boxplots", sep=""))
}

write.csv(wilcoxon_results, "Diagnostic genes: control vs D1.csv", row.names = FALSE)
write.csv(results_fdr, "Diagnostic FDR genes: control vs D1.csv", row.names = FALSE)
