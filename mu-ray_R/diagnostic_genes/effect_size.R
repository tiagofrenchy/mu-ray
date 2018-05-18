# O valor de p só diz a significância -- aponta a diferença.
# Como selecionar os genes com maior potencial para focalizarmos a pesquisa? Pelo tamanho do efeito.
# Aqui, calculado como a distância entre boxplots dividida pelo tamanho médio dos boxplots.


#### NÃO ESTÁ SENDO FEITO PARA LOGARITMO DAS EXPRESSÕES, MAS PARA VALORES ABSOLUTOS.
#### Isso pode ser modificado com poucas linhas.
library(readr)

# Vittorio
setwd("~/IC Alexandre")

# Tiago
#setwd("~/mu-ray_data/")



# Abre os dados que serão utilizados
clean_data <- read_csv("~/IC Alexandre/clean_data_without_duplicates.csv")
regulation_fdr <- read_csv("~/IC Alexandre/Regulação dos genes - fdr valores absolutos.csv")


# Criamos uma coluna para o tamanho do efeito -- calculado como distânia entre os boxplots
# (descontados outliers) dividido pela média do "tamanho" dos boxplots.
regulation_fdr["effect_size"] <- 0

listUniqueTranscripts <- unique(regulation_fdr$transcript_cluster_id)
numberTranscripts <- length(uniqueTranscripts)

for (i in seq(numberTranscripts)) {
  control <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                   Dia == 'D0') %>% select(Dia, expression, Paciente)
  sepse  <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                   Dia == 'D1') %>% select(Dia, expression, Paciente)
  
  # Todas as medidas são baseadas em quartis (primeiro, terceiro) e a distância interquartil.
  IQ_sepse <- IQR(sepse$expression )
  IQ_ctrl  <- IQR(control$expression)
  
  Q1_sepse <- quantile(sepse$expression, 0.25)
  Q3_sepse <- quantile(sepse$expression, 0.75)
  
  Q1_ctrl  <- quantile(control$expression, 0.25)
  Q3_ctrl  <- quantile(control$expression, 0.75)
  
  
  #   Iremos dividir a distância entre os boxplots pela média do tamanho entre eles. Isso é independente
  # da posição relativa entre ambos. Calculamos abaixo.
  size_sepse <- (Q3_sepse + 1.5 * IQ_sepse) - (Q1_sepse - 1.5 * IQ_sepse)
  size_ctrl  <- (Q3_ctrl  + 1.5 * IQ_ctrl ) - (Q1_ctrl  - 1.5 * IQ_ctrl)
   
  avg_size <- (size_sepse + size_ctrl) / 2
  
  #   Para determinarmos a "distânia" entre os boxplots, precisamos saber qual está "em cima"
  # e qual está "em baixo" -- pois ela é definida pela distânia entre máximos/mínimos. É 
  # isso que faz o if/else abaixo. 
  effectSize <- 0
  
  if (regulation_fdr$media_sepse[i] > regulation_fdr$media_controle[i]){
    min_sepse <- Q1_sepse - 1.5*IQ_sepse
    max_ctrl  <- Q3_ctrl  + 1.5*IQ_ctrl
    
    effectSize <- (min_sepse - max_ctrl) / avg_size
  }
  # Se o controle está acima da sepse...
  else {
    min_ctrl  <- Q1_ctrl  - 1.5*IQ_ctrl
    max_sepse <- Q3_sepse + 1.5*IQ_sepse
    
    effectSize <- (min_ctrl - max_sepse) / avg_size
  }
  
  regulation_fdr$effect_size[i] <- effectSize
  
  # DEBUG
  print(sprintf('%5.1f%%', i/length(listUniqueTranscripts)*100), quote = FALSE)
}

write.csv(regulation_fdr, "Effect size dos diagnostios (com valores absolutos).csv", row.names = FALSE)