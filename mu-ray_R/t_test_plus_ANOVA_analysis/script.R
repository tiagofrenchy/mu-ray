library(readr)
library(dplyr)
library(car)


# Vittorio
setwd("~/IC Alexandre")

# Tiago
#setwd("~/mu-ray_data/")

# Abre os dados
clean_data <- read_csv("clean_data_without_duplicates.csv")
clean_data$expression <- log(clean_data$expression)

# Iniciamos o data frame dos resultados
results <- data.frame(transcript_cluster_id = double(), 
                      description   = character(),
                      normD1 = double(),
                      normD2 = double(),
                      normD3 = double(),
                      normD4 = double(),
                      normSepsis = double(),
                      bartSepsis = double(),
                      normControl = double(),
                      normResidual = double(),
                      p_ANOVA=double(),
                      p_tuk21=double(),
                      p_tuk31=double(),
                      p_tuk41=double(),
                      p_tuk32=double(),
                      p_tuk42=double(),
                      p_tuk43=double(),
                      p_ctrlDia1=double(),
                      p_ctrlDia2=double(),
                      p_ctrlDia3=double(),
                      p_ctrlDia4=double(),
                      p_ctrlSepsis=double(),
                      stringsAsFactors = FALSE)


listUniqueTranscripts <- unique(clean_data$transcript_cluster_id)
numberTranscripts <- length(listUniqueTranscripts)

for (i in seq(numberTranscripts)) {
  # Para o gene atual, selecionamos a sua expressão em todos os dias do grupo com sepse.
  groupExpression <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i], Dia != 'D0') %>% 
    select(Dia, expression)
  groupExpression$Dia <- ordered(groupExpression$Dia, levels = c("D1", "D2", "D3", "D4"))
  
  
  # Para o gene atual, selecioanmos a sua expressão nos controles
  controlExpression <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i], Dia == 'D0') %>% 
    select(expression)
  
  
  # Antes de computar o ANOVA, verificamos a normalidade dos dados e a homogeneidade das variâncias.
  s1 = shapiro.test((groupExpression %>% filter(Dia == 'D1'))$expression)
  s2 = shapiro.test((groupExpression %>% filter(Dia == 'D2'))$expression)
  s3 = shapiro.test((groupExpression %>% filter(Dia == 'D3'))$expression)
  s4 = shapiro.test((groupExpression %>% filter(Dia == 'D4'))$expression)
  b  = bartlett.test(groupExpression$expression, groupExpression$Dia)
  
  
  resultadoANOVA <- aov(expression ~ Dia, data = groupExpression)
  pANOVA <- summary(resultadoANOVA)[[1]][["Pr(>F)"]][[1]]
    
  residualANOVA <- residuals(object = resultadoANOVA )
  pResidual <- shapiro.test(x = residualANOVA )$p.value
  
  # Computamos o TukeyHSD para fazer várias comparações entre pares e localizar aonde
  # está a diferença (ou mais de uma, se houver).
  tuk <- data.frame(c(TukeyHSD(resultadoANOVA)))
  tuk21 <- tuk$Dia.p.adj[1]
  tuk31 <- tuk$Dia.p.adj[2]
  tuk41 <- tuk$Dia.p.adj[3]
  tuk32 <- tuk$Dia.p.adj[4]
  tuk42 <- tuk$Dia.p.adj[5]
  tuk43 <- tuk$Dia.p.adj[6]

  # Realiza testes T entre a média da expressão dos controles e cada dia da sepse.
  pNormControl <- shapiro.test(controlExpression$expression)$p.value
  pNormSepsis <- shapiro.test(groupExpression$expression)$p.value
 
  controleDia1 = t.test(controlExpression$expression, 
                              (groupExpression %>% filter(Dia == 'D1'))$expression)

  controleDia2 = t.test(controlExpression$expression, 
                              (groupExpression %>% filter(Dia == 'D2'))$expression)

  controleDia3 = t.test(controlExpression$expression, 
                              (groupExpression %>% filter(Dia == 'D3'))$expression)

  controleDia4 = t.test(controlExpression$expression, 
                              (groupExpression %>% filter(Dia == 'D4'))$expression)

  controleSepse = t.test(controlExpression$expression, groupExpression$expression)
  
  # Primeiro, obtemos a descrição do gene.
  descr <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i]) %>% 
    select(description)
  
  # Agora, inserimos todos os valores calculados no data frame.
  results[nrow(results) + 1,]     = c(listUniqueTranscripts[i],
                                      c(descr[1][1,]),
                                      s1$p.value,
                                      s2$p.value,
                                      s3$p.value,
                                      s4$p.value,
                                      pNormSepsis,
                                      b$p.value,
                                      pNormControl,
                                      pResidual,
                                      pANOVA,
                                      tuk21,
                                      tuk31,
                                      tuk41,
                                      tuk32,
                                      tuk42,
                                      tuk43,
                                      controleDia1$p.value,
                                      controleDia2$p.value,
                                      controleDia3$p.value,
                                      controleDia4$p.value,
                                      controleSepse$p.value)
  
  
  # DEBUG
  print(sprintf('%5.1f%%', i/length(listUniqueTranscripts)*100), quote = FALSE)
  

}

write.csv(results, "Resultado dos testes estatísticos.csv", row.names = FALSE)
