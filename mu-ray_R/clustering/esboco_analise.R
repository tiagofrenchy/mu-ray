    #
    ###
######### Aqui não tem nada ainda mas pretendo trabalhar com fold-change, DEG,
######### coexpression network e up and down regulation
    ###
    #
  
# 
# Primeiro precisamos comparar a mudanca da media em todos os dias da sepse (e total)
# com a media dos controles para encontrar o fold change e determinar se a mudanca
# é positiva ou negativa (up and down-regulation). 
# 
# 
# 
# Depois, é necessário verificar se a mudanca é significativa (t test ou equivalente).
# Constatando diferenca significativa estabelecemos o grau de mudanca por um fator
# multiplicativo (!!! cuidado com o fato de estar usando log !!!)
# 
# aqui nós listamos todos os genes que são diferencialmente expressos (DEG)
# 
# Determinado o fold-change, podemos analizar covarianca de expressao por clustering 
# (talvez mesmo promotor ou fator estimulante). 
# 
# 
# 

library(ggplot2)
library(reshape2)
library(dplyr)

# Vittorio
# setwd("~/IC Alexandre")

# Tiago
setwd("~/mu-ray_data/")

regulation  <- read_csv("Regulação dos genes.csv")
regulation2 <- select(regulation, "media_controle", "media_sepse", "diferenca_medias")
regulation2$diferenca_medias <- abs(regulation2$diferenca_medias)
regulation2 <- melt(regulation2, id.vars="diferenca_medias")

ggplot(regulation2, aes(value, diferenca_medias, col=variable)) +
  geom_jitter(size=1.5, shape=20) +
  xlab("Change in Expression") +
  ylab("Expression of Genes") +
  stat_smooth()

  # ggtitle(nomeTit) +
  # #theme(legend.position="none", plot.title = element_text(family="Times", colour="black", size=6)) +
  # theme(legend.position="none",
  #       title = element_text(size = tamTit),
  #       axis.text = element_text(size = tamResto),
  #       axis.title = element_text(size = tamResto)) +
  # scale_fill_brewer(palette = "Set1") +
  # ggsave(sprintf("[%s] - %s.png",id, nome), path = paste(getwd(), "/boxplots", sep=""))
