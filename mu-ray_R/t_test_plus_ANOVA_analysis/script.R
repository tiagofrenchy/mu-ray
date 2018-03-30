library(readr)
library(dplyr)

# Vittorio
setwd("~/IC Alexandre")

# Tiago
# setwd("~/Downloads/")

# Abre os dados
clean_data <- read_csv("clean_data_without_duplicates.csv")

# Só estamos interessados nos controles
controle <- filter(clean_data, Dia == 'D0')
dia_1 <- filter(clean_data, Dia == 'D1')
dia_2 <- filter(clean_data, Dia == 'D2')
dia_3 <- filter(clean_data, Dia == 'D3')
dia_4 <- filter(clean_data, Dia == 'D4')


retornaMedia <- function(grupo){
  # Calculamos o log de todas as expressões.
  grupo["log_expression"] <- log(grupo$expression)
  
  # Cria um novo DF com 2 colunas. Para cada um dos genes (transcript_cluster), teremos
  # uma média do log da expressão (tal qual herdada do data frame antigo).
  log_averages <- data.frame(transcript_cluster_id=double(), 
                             average_expression=double(), 
                             stringsAsFactors = FALSE)
  
  # Calculamos o número de transcriptos que existem para serem analisados.
  listUniqueTranscripts <- unique(grupo$transcript_cluster_id)
  numberTranscripts <- length(listUniqueTranscripts)
  
  # No loop abaixo, para cada um de todos os genes do grupo (controle ou dia x da sepse)
  # calculamos a média do log de sua expressão.
  for (i in seq(numberTranscripts)) {
    # Filtramos nosso data frame de modo que só existirão observações do gene que estamos observando
    # no momento e de suas mensurações entre os controles (lembre-se: para cada controle, há um único
    # momento de observação: D0)
    dado <- filter(grupo, transcript_cluster_id == listUniqueTranscripts[i])
    
    media_log <- mean(dado$log_expression)
    
    # DEBUG
    print(sprintf('%8.3f  %5.1f%%', media_log, i/numberTranscripts*100), quote = FALSE)
    
    log_averages [nrow(log_averages) + 1,] = c(listUniqueTranscripts[i], media_log)
  }
  
  return (log_averages)
}

media_controle <- retornaMedia(controle)
media_dia_1 <- retornaMedia(dia_1)
media_dia_2 <- retornaMedia(dia_2)
media_dia_3 <- retornaMedia(dia_3)
media_dia_4 <- retornaMedia(dia_4)



############
# TODO: sprintf ? > criar nome diferenciado para cada grupo diferente na hora de salvar csv
# # Ninguém merece calcular isso de novo. 
# write.csv(log_averages, "media do LOG da expressao dentre os controles.csv", row.names = FALSE)
############



#t.test(tip ~ sex, data = tips)$p.value
