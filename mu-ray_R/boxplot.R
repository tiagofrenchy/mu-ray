# Monta o boxplot da expressão do gene para o CTRL e o D1
# a partir do transcript_id (você precisa inserir manualmente
# e roar o código). 

### ATENÇÃO: *não* utiliza o logaritmo da expressão!

library(ggplot2)
library(readr)


# Vittorio
setwd("~/IC Alexandre")

# Tiago
#setwd("~/mu-ray_data/")

# Abre os dados
clean_data <- read_csv("clean_data_without_duplicates.csv")

# INSIRA ABAIXO:
your_gene <- 8103769

sepse <- clean_data %>% filter (transcript_cluster_id == your_gene & Dia == 'D1')
ctrl  <- clean_data %>% filter (transcript_cluster_id == your_gene & Dia == 'D0')
ctrl$Dia <- gsub("D0", "Controle", ctrl$Dia)

nome <- sepse$description[1]

ggplot() +
  geom_boxplot(data = ctrl, aes(x=Dia, y=expression, fill=Dia), alpha=0.5) +
  geom_boxplot(data = sepse, aes(x=Dia, y=expression, fill=Dia), alpha=0.5) +
  xlab("Grupos") +
  ylab("Expressão") +
  ggtitle(nome) +
  theme(legend.position="none") +
  scale_fill_brewer(palette = "Set1")

