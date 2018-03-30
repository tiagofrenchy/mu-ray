library(readr)
library(dplyr)
library(gridExtra)
library(ggplot2)

# Vittorio
# setwd("~/IC Alexandre")

#Tiago
#setwd("~/mu-ray_data/")


# Abre os dados
clean_data <- read_csv("clean_data.csv")

# Só estamos interessados nos controles
controle <- filter(clean_data, Dia == 'D0')

# Calculamos o log de todas as expressões.
controle["log_expression"] <- log(controle$expression)



# Cria um novo DF com 2 colunas. Para cada um dos genes (transcript_cluster), teremos
# uma média da expressão (tal qual herdada do data frame antigo).
#
# Faz o mesmo para média do logaritmo da expressão. 
averages <- data.frame(transcript_cluster_id=double(), 
                       average_expression=double(),
                       stringsAsFactors = FALSE)

log_averages <- data.frame(transcript_cluster_id=double(), 
                           average_expression=double(), 
                           stringsAsFactors = FALSE)


# Calculamos o número de transcriptos (únicos) que existem para serem analisados.
listUniqueTranscripts <- unique(controle$transcript_cluster_id)
numberTranscripts <- length(listUniqueTranscripts)



# No loop abaixo, para cada um de todos os genes de nosso transcriptoma, calculamos a média
# de sua expressão entre os controles (e a média do LOGARITMO de suas expressões).
for (i in seq(numberTranscripts)) {
  # Filtramos nosso data frame de modo que só existirão observações do gene que estamos observando
  # no momento e de suas mensurações entre os controles (lembre-se: para cada controle, há um único
  # momento de observação: D0)
  dado <- filter(controle, transcript_cluster_id == listUniqueTranscripts[i])
  
  media     <- mean(dado$expression)
  media_log <- mean(dado$log_expression)
  
  # DEBUG
  print(sprintf('%8.1f  %8.3f  %5.1f%%', media, media_log, i/numberTranscripts*100), quote = FALSE)
  
  averages     [nrow(averages) + 1,]     = c(listUniqueTranscripts[i], media)
  log_averages [nrow(log_averages) + 1,] = c(listUniqueTranscripts[i], media_log)
}


############
# Ninguém merece calcular isso de novo.
write.csv(averages, "media da expressao dentre os controles.csv", row.names = FALSE)
write.csv(log_averages, "media do LOG da expressao dentre os controles.csv", row.names = FALSE)
############


###########################
# Faz histograma da distribuição da média das expressões e da méia dos 
# logaritmos das expressões. As indentações são necessárias por razões estúpidas.

{
  name <- sprintf("Distribuição da média das expressões")
  png(name)
  h <- hist(averages$average_expression, xlim = c(0,1000), breaks = 1000)#, xlim = c(0,3000), breaks = 1000)
  dev.off()
}

h<-hist(x, breaks = 7, col="white", xlab="Expression Intensity",
        ylab = "# of Observations", main=sprintf("%s-%s",dado$description[1],transcript_id)) 



{
  plot.new()
  name2 <- sprintf("Distribuição da média do logaritmo das expressões")
  png(name2)
  h <- hist(log_averages$average_expression) #, xlim = c(0,5), breaks = 1000)
  x <- log_averages$average_expression
  xfit<-seq(min(x),max(x),length=50) 
  yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
  yfit <- yfit*diff(h$mids[1:2])*length(x) 
  lines(xfit, yfit, col="blue", lwd=2)
  dev.off()
}
###########################


#######################################
# Denominamos localAverages para distinguir o objeto 'averages' global -- este
# localAverages receberá em um primeiro momento 'averages' de fato e,
# em um segundo momento, 'log_averages'.
analiseExpressao <- function(localAverages) {
  # Dividir em tercil não é uma boa idéia. Embora seja muito fácil e rápido, não teríamos
  # acesso rápido às observações no "topo" de cada tercil. Assim, melhor ordenar os dados em função
  # da expressão média e agir de acordo.
  #     -> Observação: '-' sinaliza que a ordenação seja decrescente
  #
  #     -> Observação2: orderedAverages é igual averages (também possui description,
  #         transcript_cluster_id, etc) -- mas está ordenada pelas expressões.
  orderedAverages <- localAverages[order(-localAverages$average_expression),]
  
  
  # Agora, vamos criar uma nova tabela com os poucos genes que selecionamos para averiguar a 
  # normalidade (só estamos interessados nos 05 primeiros genes de cada tercil).
  selectedGenes <- data.frame(transcript_id=double(), 
                              description=character(), 
                              average_expression=double(), 
                              shapiro_p_value=double(),
                              significant=logical(),
                              stringsAsFactors = FALSE)
  
  size <- length(orderedAverages$transcript_cluster_id)
  start <- ceiling(size/3)
  start2 <- start*2
  
  selectedGenes <- analiseTopTercil(1,      orderedAverages, selectedGenes)
  selectedGenes <- analiseTopTercil(start,  orderedAverages, selectedGenes)
  selectedGenes <- analiseTopTercil(start2, orderedAverages, selectedGenes)
 

  # Uses library(gridExtra)
  png("selectedGenes.png", height = 50*nrow(selectedGenes), width = 400*ncol(selectedGenes))
  grid.table(selectedGenes)
  dev.off()
}



###########################
# TODO O que fazemos com isso?
# faz histograma das medias e log das medias, para auxiliar na divisão dos grupos
# hist(averages$average_expression, xlim = c(0,3000), breaks = 1000)
# hist(averages$log_average_expression, xlim = c(0,5), breaks = 1000)
###########################



# Começamos a preencher a tabela com os 5 primeiros genes do 1º tercil.
# Para cada um deles,
#   - Copiamos o transcript id
#   - Obtemos um vetor com todas as expressões gênicas mensuradas no grupo controle
#   - Calcular o shapiro-wilks sobre o vetor acima
#   - Geramos o histograma e salvamos
#   - Inseririmos os dados acima em uma tabela, *junto com a descrição*, retornamos o dado

analiseTopTercil <- function(inicio, orderedAverages, selectedGenes) {
  for (i in inicio:(inicio+4) ) {
    transcript_id <- orderedAverages$transcript_cluster_id[i]
    dado <- filter(controle, transcript_cluster_id == transcript_id)
    test <- shapiro.test(dado$expression)
    
    sig = FALSE
    sig[test$p.value < 0.001] <- TRUE
    
    selectedGenes[nrow(selectedGenes) + 1,] = c(transcript_id, 
                                                dado$description[1], 
                                                orderedAverages$average_expression[i],
                                                test$p.value,
                                                sig)
    name <- sprintf("Gene_número_%d_(%d)--%s", i, transcript_id, dado$description)
    png(name)
    x <- dado$expression
    h<-hist(x, breaks = 7, col="white", xlab="Expression Intensity",
            ylab = "# of Observations", main=sprintf("%s-%s",dado$description[1],transcript_id)) 
    diff <- (max(x)-min(x))/8
    xfit<-seq((min(x)-diff),(max(x)+diff),length=100) 
    yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
    yfit <- yfit*diff(h$mids[1:2])*length(x) 
    lines(xfit, yfit, col="blue", lwd=2)
    d <- density(x)
    lines(x = d$x, y = d$y * length(x) * diff(h$breaks)[1], lwd = 2, col="red")
    dev.off()
  }
  return (selectedGenes)
}

analiseExpressao(averages)
analiseExpressao(log_averages)
