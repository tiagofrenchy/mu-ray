library(readr)
library(dplyr)
library(tidyr)
library(dunn.test)
library(PMCMRplus)


# Vittorio
#setwd("~/IC Alexandre")

# Tiago
setwd("~/mu-ray_data/")

# Abre os dados
clean_data <- read_csv("clean_data_without_duplicates.csv")
clean_data$expression <- log(clean_data$expression)

# Iniciamos o data frame dos resultados
results <- data.frame(transcript_cluster_id = double(), 
                      description   = character(),
                      p_SMACK = double(),
                      p_dun12 = double(),
                      p_dun13 = double(),
                      p_dun23 = double(),
                      p_dun14 = double(),
                      p_dun24 = double(),
                      p_dun34 = double(),
                      p_ctrlDia1 = double(),
                      p_ctrlDia2 = double(),
                      p_ctrlDia3 = double(),
                      p_ctrlDia4 = double(),
                      p_ctrlSepsis = double(),
                      stringsAsFactors = FALSE)


listUniqueTranscripts <- unique(clean_data$transcript_cluster_id)
numberTranscripts <- length(listUniqueTranscripts)

for (i in seq(numberTranscripts)) {
  # Para o gene atual, selecionamos a sua expressão em todos os dias do grupo com sepse.
  sepsisExpression <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                            Dia != 'D0') %>% select(Dia, expression, Paciente)
  sepsisExpression$Dia <- ordered(sepsisExpression$Dia, levels = c("D1", "D2", "D3", "D4"))
  
  # Para o gene atual, selecionamos a sua expressão nos controles
  controlExpression <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                             Dia == 'D0') %>% select(expression, Paciente)
  
  # Por que pegamos o Paciente no select() para depois tirar? Pois é necesário ao spread()
  x <- sepsisExpression %>% spread(Dia, expression)
  x$Paciente <- NULL
  
  a <- skillingsMackTest(as.matrix(x))
  p_mack <- a$p.value
  
  
  # Computamos o DunnTest para fazer várias comparações entre pares e localizar aonde
  # está a diferença (ou mais de uma, se houver).
  sink("/dev/null")
  DUNN <- dunn.test(sepsisExpression$expression, sepsisExpression$Dia)
  dun12 <- DUNN$P[1]
  dun13 <- DUNN$P[2]
  dun23 <- DUNN$P[3]
  dun14 <- DUNN$P[4]
  dun24 <- DUNN$P[5]
  dun34 <- DUNN$P[6]
  sink()
  # Realiza teste de Wilcoxon tanto entre a média da expressão dos controles e cada dia da sepse
  # quanto entre a primeira e a média da sepse inteira (todos os dias).
  #
  # Nota: aqui os testes NÃO são pareados.
  controleDia1 = wilcox.test(controlExpression$expression, 
                        (sepsisExpression %>% filter(Dia == 'D1'))$expression, paired = FALSE,
                        exact = TRUE)
  
  controleDia1 = wilcox.test(controlExpression$expression, 
                        (sepsisExpression %>% filter(Dia == 'D1'))$expression, paired = FALSE,
                        exact = TRUE)
  
  controleDia2 = wilcox.test(controlExpression$expression, 
                        (sepsisExpression %>% filter(Dia == 'D2'))$expression, paired = FALSE,
                        exact = TRUE)
  
  controleDia3 = wilcox.test(controlExpression$expression, 
                        (sepsisExpression %>% filter(Dia == 'D3'))$expression, paired = FALSE,
                        exact = TRUE)
  
  controleDia4 = wilcox.test(controlExpression$expression, 
                        (sepsisExpression %>% filter(Dia == 'D4'))$expression, paired = FALSE,
                        exact = TRUE)
  
  controleSepse = wilcox.test(controlExpression$expression, sepsisExpression$expression, 
                              paired = FALSE, exact = TRUE)
  
  
  # Primeiro, obtemos a descrição do gene.
  descr <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i]) %>% 
    select(description)
  
  # Agora, inserimos todos os valores calculados no data frame.
  results[nrow(results) + 1,]     = c(listUniqueTranscripts[i],
                                      c(descr[1][1,]),
                                      p_mack,
                                      dun12,
                                      dun13,
                                      dun23,
                                      dun14,
                                      dun24,
                                      dun34,
                                      controleDia1$p.value,
                                      controleDia2$p.value,
                                      controleDia3$p.value,
                                      controleDia4$p.value,
                                      controleSepse$p.value)
  
  
  # DEBUG
  print(sprintf('%5.1f%%', i/length(listUniqueTranscripts)*100), quote = FALSE)
  
  
}

write.csv(results, "Resultado dos testes estatísticos NON-PARAMETRIC.csv", row.names = FALSE)
