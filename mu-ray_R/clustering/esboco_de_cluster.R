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

# Iniciamos o data frame da regulacao
regulated <- data.frame(transcript_cluster_id = double(), 
                      description   = character(),
                      p_SMACK = double(),
                      p_dun12 = double(),
                      p_dun13 = double(),
                      p_dun23 = double(),
                      p_dun14 = double(),
                      p_dun24 = double(),
                      p_dun34 = double(),
                      p_ctrlDia1 = double(),
                      p_ctrlDia2 = double(),
                      p_ctrlDia3 = double(),
                      p_ctrlDia4 = double(),
                      p_ctrlSepsis = double(),
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