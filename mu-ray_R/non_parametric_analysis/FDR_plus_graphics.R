library(readr)
library(sgof)
library(dplyr)
library(ggplot2)

# Vittorio
setwd("~/IC Alexandre")

# Tiago
# setwd("~/mu-ray_data/")

# Abre os dados
clean_data <- read_csv("clean_data_without_duplicates.csv")
clean_data$expression <- log(clean_data$expression)

results <- read_csv("~/IC Alexandre/Resultado dos testes estatísticos NON-PARAMETRIC.csv")

fdr_smack  <- BH(results$p_SMACK, alpha = 0.001)
results <- results[order(results$p_SMACK),]
id_fdr <- head(results$transcript_cluster_id, n = fdr_smack$Rejections)


fdr_ctrl   <- BH(results$p_ctrlSepsis, alpha = 0.001)
results <- results[order(results$p_ctrlSepsis),]
id_ctrl <- head(results$transcript_cluster_id, n = fdr_ctrl$Rejections)


coolGenesId <- intersect(id_fdr, id_ctrl)

coolGenes <- clean_data %>% filter(transcript_cluster_id %in% coolGenesId)


listUniqueTranscripts <- unique(coolGenes$transcript_cluster_id)
numberTranscripts <- length(listUniqueTranscripts)

setwd("~/IC Alexandre/Imgs/")
for (i in seq(numberTranscripts)) {
  reallyCoolGenes <- coolGenes %>% filter(transcript_cluster_id == listUniqueTranscripts[i])
  
  nome <- reallyCoolGenes$description[1]
  id <- factor(listUniqueTranscripts[i])
  png(sprintf("[%s] - %s .png",id, nome))
  
  img <- ggplot(reallyCoolGenes, aes(x=Dia, y=expression, fill=Dia)) +
    geom_boxplot(alpha=0.3) +
    theme(legend.position="none") +
    scale_fill_brewer(palette="Dark2")
  
  print(img)
  dev.off()
}


#######3


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