library(readr)
library(dplyr)
library(ggpubr)

# Vittorio
setwd("~/IC Alexandre")

# Tiago
# setwd("~/Downloads/")

# Abre os dados
data <- read_csv("clean_data_without_duplicates.csv")

data$expression <- log(data$expression)

listUniqueTranscripts <- unique(data$transcript_cluster_id)
numberTranscripts <- length(listUniqueTranscripts)

# 
# for (i in seq(numberTranscripts)) {
#   dado <- filter(data, transcript_cluster_id == listUniqueTranscripts[i])
#            %>% select(data, Dia, expression)

  
dado <- data %>% filter(transcript_cluster_id == listUniqueTranscripts[1], Dia != 'D0') %>% 
  select(Dia, expression)
dado$Dia <- ordered(dado$Dia, levels = c("D1", "D2", "D3", "D4"))

# Cria boxplot das espressões
# ggboxplot(dado, x = "Dia", y = "expression",
#           color = "Dia", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#000000"),
#           order = c("D1", "D2", "D3", "D4"),
#           ylab = "Expressão", xlab = "Dia")

# Compute the analysis of variance
res.aov <- aov(expression ~ Dia, data = data)
# Summary of the analysis
summary(res.aov)


  
  ###############  
  media_log <- mean(dado$log_expression)
  
  # DEBUG
  print(sprintf('%8.3f  %5.1f%%', media_log, i/numberTranscripts*100), quote = FALSE)
  
  log_averages [nrow(log_averages) + 1,] = c(listUniqueTranscripts[i], media_log)
}

return (log_averages)
}