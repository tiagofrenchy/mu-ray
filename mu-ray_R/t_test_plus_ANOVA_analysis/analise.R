library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Vittorio
# setwd("~/IC Alexandre")

# Tiago
#setwd("~/mu-ray_data/")

# Abre os dados.
clean_data <- read_csv("clean_data_without_duplicates.csv")
clean_data$expression <- log(clean_data$expression)

results <- read_csv("Resultado dos testes estatísticos.csv")

p_value = .05

#### Selecionando genes que respeitam os pressupostos do ANOVA
listUniqueTranscripts <- unique(results$transcript_cluster_id)
numberTranscripts <- length(listUniqueTranscripts)

for (i in seq(numberTranscripts)) {
   
  results$"segue todos os pressupostos ANOVA"[i]  <- (results$normD1[i] < p_value &&
                                                      results$normD2[i] < p_value &&
                                                      results$normD3[i] < p_value &&
                                                      results$normD4[i] < p_value &&
                                                      results$normSepsis[i] < p_value &&
                                                      results$bartSepsis[i] < p_value &&
                                                      results$normResidual[[i]] < p_value)
  

  results$"segue todos os pressupostos T-test"[i] <- (results$normD1[i] < p_value &&
                                                      results$normD2[i] < p_value &&
                                                      results$normD3[i] < p_value &&
                                                      results$normD4[i] < p_value &&
                                                      results$normSepsis[i] < p_value &&
                                                      results$normControl[i] < p_value) 
}

results <- results[c("transcript_cluster_id", "description", "normD1", "normD2", "normD3",
                     "normD4", "normSepsis", "bartSepsis", "normControl", "normResidual",
                     "segue pressupostos ANOVA", "p_ANOVA", "p_tuk21", "p_tuk31", "p_tuk41",
                     "p_tuk32", "p_tuk42", "p_tuk43", "segue pressupostos T-test",
                     "p_ctrlDia1", "p_ctrlDia2", "p_ctrlDia3", "p_ctrlDia4", "p_ctrlSepsis")]

#### Genes estudando
dipep <- results %>% filter(transcript_cluster_id == 8002181)
tnfr <- results %>% filter(transcript_cluster_id == 8149733)

difCtrl2 <- results %>% filter(p_ctrlSepsis < .005 / 27000)
difANOVA <- difCtrl2 %>% filter(p_ANOVA < p_value / 27000)


#### Gráfico

gene <- 7899534
dadosGene <- clean_data %>% filter(transcript_cluster_id == gene)


##### FUNÇÃO GRÁFICA 02
boxplot(expression ~ Dia, data = dadosGene)

# Cria boxplot das espressões
# ggboxplot(dado, x = "Dia", y = "expression",
#           color = "Dia", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#000000"),
#           order = c("D1", "D2", "D3", "D4"),
#           ylab = "Expressão", xlab = "Dia")




###### FUNÇÃO GRÁFICA 01, IMPORTADA DO ARQUIVO ANTIG
# Obtemos nome e descrição do gene.
nomeT <- clean_data %>% filter(transcript_cluster_id == gene) %>% select(description)
nome <- sprintf("%s", nomeT[1][1,])


# Agora, listamos as comparações que serão realizadas.
my_comparisons <- list( c("D1", "D2"), c("D2", "D3"),
                        c("D3", "D4"), c("D1", "D3"), c("D1", "D4"),
                        c("D2", "D4"))

# Por fim, fazemos o gráfico.
setwd("~/IC Alexandre/Imgs/")

png(sprintf("%s.png",nome))

groupExpression <- clean_data %>% filter(transcript_cluster_id == gene, Dia != 'D0') %>% 
  select(Dia, expression)
groupExpression$Dia <- ordered(groupExpression$Dia, levels = c("D1", "D2", "D3", "D4"))

controlExpression <- clean_data %>% filter(transcript_cluster_id == gene, Dia == 'D0') %>% 
  select(expression)

# Calculamos a altura máxima do boxplot para que o gráfico fique com proporções
# decentes.
altura = max(groupExpression$expression) + 0.1

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