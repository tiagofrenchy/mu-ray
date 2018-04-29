# O valor de p só diz a significância -- aponta a diferença.
# Como selecionar os genes com maior potencial para focalizarmos a pesquisa? Pelo tamanho do efeito.
# Aqui, calculado como a distância entre boxplots dividida pelo tamanho médio dos boxplots.
# (supõe que o results_fdr_log está carregado na memória)

results_fdr_log["diff_max_min"] <- 0

uniqueTranscripts <- unique(results_fdr_log$transcript_cluster_id)
numberTranscripts <- length(uniqueTranscripts)

for (i in seq(numberTranscripts)) {
  control <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                   Dia == 'D0') %>% select(Dia, expression, Paciente)
  sepsis  <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                   Dia == 'D1') %>% select(Dia, expression, Paciente)
  
  if (results_fdr_log$media_log_sepse[i] > results_fdr_log$media_log_controle[i]){
    min_sepse <- min( log(sepsis$expression ) )
    max_ctrl  <- max( log(control$expression) )
    diff <- min_sepse - max_ctrl
  }
  else {
    min_ctrl  <- min( log(control$expression))
    max_sepse <- max( log(sepsis$expression )) 
    diff <- min_ctrl - max_sepse
  }
  
  results_fdr_log$diff_max_min[i] <- diff
  print(i)
}

results_fdr_log["diff_max_min_adjusted"] <- 0
for (i in seq(numberTranscripts)) {
  control <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                   Dia == 'D0') %>% select(Dia, expression, Paciente)
  sepsis  <- clean_data %>% filter(transcript_cluster_id == listUniqueTranscripts[i],
                                   Dia == 'D1') %>% select(Dia, expression, Paciente)
  
  sizeCtrl <- max( log(control$expression) ) - min( log(control$expression) )
  sizeSeps <- max( log(sepsis$expression) ) - min( log(sepsis$expression) )
  
  media = (sizeCtrl + sizeSeps) / 2
  
  results_fdr_log$diff_max_min_adjusted[i] <- (results_fdr_log$diff_max_min[i] / media)
  
  print(i)
}
