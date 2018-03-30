library(readr)
library(dplyr)
library(ggpubr)
library(Rcmdr)
library(nortest)

# Vittorio
setwd("~/IC Alexandre")

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

  
groupExpression <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[1], Dia != 'D0') %>% 
  select(Dia, expression)
groupExpression$Dia <- ordered(dado$Dia, levels = c("D1", "D2", "D3", "D4"))

# Cria boxplot das espressões
# ggboxplot(dado, x = "Dia", y = "expression",
#           color = "Dia", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#000000"),
#           order = c("D1", "D2", "D3", "D4"),
#           ylab = "Expressão", xlab = "Dia")


# Antes de computar o ANOVA, verificamo a normalidade dos dados e a homogeneidade das variâncias.


shapiro.test(groupExpression$expression)
leveneTest(expression ~ Dia, data = groupExpression)


################
#h<-hist(groupExpression$expression, breaks = 7, col="white")
###########################





# Abaixo, computamos o ANOVA.
res.aov <- aov(expression ~ Dia, data = data)
pvalue <- summary(res.aov)[[1]][["Pr(>F)"]][[1]]

if (pvalue<.05) {
  # Computamos o TukeyHSD para fazer várias comparações entre pares e localizar aonde
  # está a diferença (ou mais de uma, se houver).
  TukeyHSD(res.aov)
}


# 
#   The ANOVA test assumes that, the data are normally distributed
# and the variance across groups are homogeneous. We can check that with some diagnostic plots.
# 
# !!! Cheque a normalidade
#
# E a normalidade dos resíduos?
# Extract the residuals
# aov_residuals <- residuals(object = res.aov )
# # Run Shapiro-Wilk test
# shapiro.test(x = aov_residuals )
#
# Se não der certo, rode um Kruskall Wallis
# kruskal.test(weight ~ group, data = my_data)



# 
#   ###############  
#   media_log <- mean(dado$log_expression)
#   
#   # DEBUG
#   print(sprintf('%8.3f  %5.1f%%', media_log, i/numberTranscripts*100), quote = FALSE)
#   
#   log_averages [nrow(log_averages) + 1,] = c(listUniqueTranscripts[i], media_log)
# }
# 
# return (log_averages)
# }