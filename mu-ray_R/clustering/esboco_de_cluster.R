    #
    ###
######### Aqui não tem nada ainda mas pretendo trabalhar com fold-change, DEG,
######### coexpression network e up and down regulation
    ###
    #
library(readr)
library(dplyr)
library(tidyr)


# Vittorio
#setwd("~/IC Alexandre")

# Tiago
setwd("~/mu-ray_data/")

# Abre os dados
clean_data <- read_csv("clean_data_without_duplicates.csv")
clean_data$expression <- log(clean_data$expression)

listUniqueTranscripts <- unique(clean_data$transcript_cluster_id)
numberTranscripts <- length(listUniqueTranscripts)

# 
# Primeiro precisamos comparar a mudanca da media em todos os dias da sepse (e total)
# com a media dos controles para encontrar o fold change e determinar se a mudanca
# é positiva ou negativa (up and down-regulation. 
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


# Iniciamos o data frame da regulacao
regulated <- data.frame(transcript_cluster_id = double(), 
                      description   = character(),
                      p_SMACK = double(),
                      stringsAsFactors = FALSE)


for (i in seq(numberTranscripts)) {
  
  # Para o gene atual, selecionamos a sua expressão em todos os dias do grupo com sepse.
  sepsisExpression <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                            Dia != 'D0') %>% select(Dia, expression, Paciente)
  
  # Para o gene atual, selecionamos a sua expressão em todos os dias do grupo controle.
  controlExpression <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                            Dia == 'D0') %>% select(Dia, expression, Paciente)

  
  # 
  
    
}