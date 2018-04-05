library(readr)
library(dplyr)
library(car)
library(ggplot2)
#library(ggsignif)
library(ggpubr)


# Vittorio
setwd("~/IC Alexandre")
#changes were made
# Tiago
#setwd("~/mu-ray_data/")

# Abre os dados
clean_data <- read_csv("clean_data_without_duplicates.csv")

clean_data$expression <- log(clean_data$expression)

listUniqueTranscripts <- unique(clean_data$transcript_cluster_id)
numberTranscripts <- length(listUniqueTranscripts)

# 
# for (i in seq(numberTranscripts)) {
#   dado <- filter(data, transcript_cluster_id == listUniqueTranscripts[i])
#            %>% select(data, Dia, expression)

for (i in seq(numberTranscripts)) {
  groupExpression <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i], Dia != 'D0') %>% 
    select(Dia, expression)
  groupExpression$Dia <- ordered(groupExpression$Dia, levels = c("D1", "D2", "D3", "D4"))
  
  # Cria boxplot das espressões
  # ggboxplot(dado, x = "Dia", y = "expression",
  #           color = "Dia", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#000000"),
  #           order = c("D1", "D2", "D3", "D4"),
  #           ylab = "Expressão", xlab = "Dia")
  
  
  # Antes de computar o ANOVA, verificamos a normalidade dos dados e a homogeneidade das variâncias.
  s = shapiro.test(groupExpression$expression)
  l = leveneTest(expression ~ Dia, data = groupExpression)
  
  if (s$p.value > 0.05 && l$`Pr(>F)` > 0.05) {
    # Abaixo, computamos o ANOVA.
    res.aov <- aov(expression ~ Dia, data = groupExpression)
    pvalue <- summary(res.aov)[[1]][["Pr(>F)"]][[1]]
    
    if (pvalue<.05) {
      # Computamos o TukeyHSD para fazer várias comparações entre pares e localizar aonde
      # está a diferença (ou mais de uma, se houver).
      #TukeyHSD(res.aov)
      
      # Agora, verificamos a normalidade dos resíduos -- outro pressuposto necessário à
      # validade do ANOVA.
      aov_residuals <- residuals(object = res.aov )
      shapiro.test(x = aov_residuals )
      
      
      ###  Monta o gráfico
      # Primeiro, obtemos o nome do gene.
      nomeT <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i]) %>% 
        select(description)
      nome <- sprintf("%s", nomeT[1][1,])
      
      # Agora, listamos as comparações que serão realizadas.
      my_comparisons <- list( c("D1", "D2"), c("D2", "D3"), 
                              c("D3", "D4"), c("D1", "D3"), c("D1", "D4"),
                              c("D2", "D4"))
      
      # Calculamos a altura máxima do boxplot para que o gráfico fique com proporções
      # decentes.
      altura = max(groupExpression$expression) + 0.1
      
      # Por fim, fazemos o gráfico.
      setwd("~/IC Alexandre/Imgs/")
      png(sprintf("%s.png",nome))
      img <- ggboxplot(groupExpression, main = nome, x = "Dia", y = "expression",
                       color = "Dia", palette = "jco")+
        stat_compare_means(comparisons = my_comparisons)+
        stat_compare_means(method = "anova", label.y = altura)
      plot(img)
      dev.off()
      
      
      # DEBUG
      print(i)
    }
  }
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