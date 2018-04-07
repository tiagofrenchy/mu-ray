library(readr)
library(dplyr)
library(car)
library(ggplot2)
library(ggpubr)

#library(ggsignif)


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
  
  
  # Cria boxplot das espressões
  # ggboxplot(dado, x = "Dia", y = "expression",
  #           color = "Dia", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#000000"),
  #           order = c("D1", "D2", "D3", "D4"),
  #           ylab = "Expressão", xlab = "Dia")
  
  
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






####### MÉTODOS GRAFICOS
# 
# # Precisa corrigir o parâmetro recebido
# grafico <- function(???) {
#   ###  Monta o gráfico
#   # Primeiro, obtemos o nome do gene.
#   nomeT <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i]) %>% 
#     select(description)
#   nome <- sprintf("%s", nomeT[1][1,])
#   
#   # Agora, listamos as comparações que serão realizadas.
#   my_comparisons <- list( c("D1", "D2"), c("D2", "D3"), 
#                           c("D3", "D4"), c("D1", "D3"), c("D1", "D4"),
#                           c("D2", "D4"))
#   
#   # Calculamos a altura máxima do boxplot para que o gráfico fique com proporções
#   # decentes.
#   altura = max(groupExpression$expression) + 0.1
#   
#   # Por fim, fazemos o gráfico.
#   setwd("~/IC Alexandre/Imgs/")
#   png(sprintf("%s.png",nome))
#   img <- ggboxplot(groupExpression, main = nome, x = "Dia", y = "expression",
#                    color = "Dia", palette = "jco")+
#     stat_compare_means(comparisons = my_comparisons, label="p.signif", hide.ns = TRUE)+
#     stat_compare_means(method = "anova", label.y = altura)+
#     stat_compare_means(label = "p.signif", method = "t.test",
#                        ref.group = controlExpression)   
#   ### TODO -- pode ser t.test acima? Não precisa ser welch? Não checamos variância...
#   plot(img)
#   dev.off()
#   
  # Nos links abaixo, parece que descobriram como esconder a comparação
  # quando não é significante. Mas tá dando trabalho pra caralho para entender
  # # o que fizeram
  # https://stackoverflow.com/questions/45476950/r-ggplot2-boxplots-ggpubr-stat-compare-means-not-working-properly
  # https://stackoverflow.com/questions/45552715/r-ggplot2-perform-pairwise-tests-per-pair-in-a-facet-and-show-the-p-values-wit
  # https://stackoverflow.com/questions/46446392/r-ggplot2-boxplots-with-significance-level-more-than-2-groups-kruskal-test-an    
  
  
  # DEBUG
  
  
}

# 
# stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   stat_compare_means(label.y = 50)     # Add global p-value
# 
# 
# 
# ggboxplot(ToothGrowth, x = "dose", y = "len",
#           color = "dose", palette = "jco")+ 
#   stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+
#   stat_compare_means(label.y = 45)
# 
# 
# ggboxplot(groupExpression, x = "Dia", y = "expression",
#           color = "Dia", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#000000"))+ 
#   stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+
#   stat_compare_means(label.y = 45)
# 
# # ggboxplot(dado, x = "Dia", y = "expression",
# #           color = "Dia", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#000000"),
# #           order = c("D1", "D2", "D3", "D4"),
# #           ylab = "Expressão", xlab = "Dia")
# 
# #
# # Se não der certo, rode um Kruskall Wallis
# # kruskal.test(weight ~ group, data = my_data)
# 
# 
# 
# ################
# #h<-hist(groupExpression$expression, breaks = 7, col="white")
# ################
# 
# 
# # 
# #   ###############  
# #   media_log <- mean(dado$log_expression)
# #   
# #   # DEBUG
# #   print(sprintf('%8.3f  %5.1f%%', media_log, i/numberTranscripts*100), quote = FALSE)
# #   
# #   log_averages [nrow(log_averages) + 1,] = c(listUniqueTranscripts[i], media_log)
# # }
# # 
# # return (log_averages)
# # }