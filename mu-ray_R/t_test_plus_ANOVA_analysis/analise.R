library(readr)
library(dplyr)



# Vittorio
setwd("~/IC Alexandre")

# Tiago
#setwd("~/mu-ray_data/")

# Abre os dados
clean_data <- read_csv("clean_data_without_duplicates.csv")
clean_data$expression <- log(clean_data$expression)



#### Genes estudando
dipep <- results %>% filter(transcript_cluster_id == 8002181)

tnfr <- results %>% filter(transcript_cluster_id == 8149733)

difCtrl2 <- results %>% filter(p_ctrlSepsis < .10^(-8)) #0.05 * 10^(-3)


#### Gráfico

T <- 8103769
dadosGene <- clean_data %>% filter(transcript_cluster_id == T)


X <- hist(dadosGene$expression)
plot(X)

img2 <- boxplot(expression ~ Dia, data = dadosGene)

plot(img2)




####### MÉTODOS GRAFICOS
# Primeiro, obtemos o nome do gene.
nomeT <- clean_data %>% filter(transcript_cluster_id == T) %>% 
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

#png(sprintf("%s.png",nome))

img2 

img <- ggboxplot(groupExpression, main = nome, x = "Dia", y = "expression",
                   color = "Dia", palette = "jco")+
    stat_compare_means(comparisons = my_comparisons, label="p.signif", hide.ns = TRUE)+
    stat_compare_means(method = "anova", label.y = altura)+
    stat_compare_means(label = "p.signif", method = "t.test",
                       ref.group = controlExpression)

### TODO -- pode ser t.test acima? Não precisa ser welch? Não checamos variância...
plot(img)
dev.off()
  
  
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