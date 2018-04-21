library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(gsubfn)

# Vittorio
# setwd("~/IC Alexandre")

#Tiago
#setwd("~/mu-ray_data/")


transcript_data <- read_csv("HuGene-1_0-st-v1.na36.hg19.transcript.csv",
                            comment = "#")
# Open trancript csv from affymetrix


separated_assignment <- transcript_data %>%
  separate(mrna_assignment, c("gene_ref", "ref_source", "description"),
           sep = "//", remove = "TRUE", extra = "drop", fill = "right")
#   A partir da coluna mrna_assignment, que possui as informações atualizadas sobre os genes, retiramos
# três tipos de informações que desejamos, distribuindo-as em colunas novas: 
# gene_ref (número de referência do gene), ref_source (de onde saiu a informação) e description (um breve
# detalhamento sobre a função do gene)


summarized_table <- select(separated_assignment, "transcript_cluster_id",
                           "ref_source", "gene_ref", "description", "category")
# Importamos as três colunas que criamos acima para um dataset novo. Além disso, importamos
# category que nos diz se o gene é um controle (iremos eliminar estes depois)


nao_main <- summarized_table[ grep("main", summarized_table$category, invert = TRUE) , ]
# seleciona todas as linhas com "category" != main

summarized_table <- summarized_table[ grep("main", summarized_table$category, invert = FALSE) , ]
# seleciona todas as linhas com "category" = main

microarray_sepsis <- read_excel("110914 RMA-Gene-Linear with annotation g.xls")
# write.csv(microarray_sepsis_provisorio, "110914 RMA-Gene-Linear with annotation g.csv")
# microarray_sepsis <- read_csv("110914 RMA-Gene-Linear with annotation g.csv")
# names(microarray_sepsis) <- gsub(x = names(microarray_sepsis)(1), pattern = ", P", replacement = ", D0P")
# # Abrimos os dados de microarray do projeto. Transformamos xls em csv
# # Renomeamos os pacientes do controle, inserindo D0 no começo, antes de "P".
# names(microarray_sepsis) <- gsub(x = names(microarray_sepsis), pattern = ",P", replacement = ",P0")  


names(microarray_sepsis) <- gsub(x=names(microarray_sepsis), pattern = "^P", replacement = "D0P")
names(microarray_sepsis) <- gsub(x=names(microarray_sepsis), pattern = "^D0Probe", replacement = "Probe")
#TODO: juntar duas linhas acima em uma regex


microarray_sepsis_clean <- microarray_sepsis %>%
  gather(observation, expression, -"Probe Set ID", -"mRna - Description") %>%
  separate(observation, c("Dia", "Paciente"), sep = 2, convert = TRUE)
# O comando gather pega múltiplas colunas e as colapsa em pares key-value, duplicando as demais
# como necessário.
#   Os primeiros parâmetros são as colunas que queremos criar: 
#       - observation (que terá o número do paciente e o dia)[isso atualmente é o header]
#       - expression (que terá a expressão do gene)[isso atualmente está sob o header]
#   Os últimos parâmetros são as colunas que queremos retirar do data set (não utilizar).
#
#
# O comando separate separa as observações em "Dia" e "Paciente". O argumento sep=2 indica que 
# a separação deve ser feita após dois caracteres. O argumento convert=TRUE indica que queremos
# converter as strings em variáveis numéricas quando possível.
#
# Nos utilizamos a saída do separate para alimentar o elemento "observation" do gather.
#
# (final result: put observations in rows and variables in columns in the microarray data)


final_sepsis_table <- merge(summarized_table, microarray_sepsis_clean,
                            by.x = "transcript_cluster_id", by.y = "Probe Set ID")
# Junta nossos dois data sets (o que retiramos dos dados da AFFX e o que modificamos dos dados
# do projeto) em um único novo elemento.


# Deletamos todos os genes que são controle?
# 
# length(unique(microarray_sepsis$"Probe Set ID"))
# length(unique(nao_main$transcript_cluster_id))
# length(unique(summarized_table$transcript_cluster_id))
# length(unique(final_sepsis_table$transcript_cluster_id))
# 
# a <- microarray_sepsis$`Probe Set ID`[(microarray_sepsis$`Probe Set ID` %in% nao_main$transcript_cluster_id)]
# #View(a)
# hello <- merge(nao_main, a, by.x = "transcript_cluster_id", by.y = 1)
# View(hello)
# # comando para achar os valores em comum nas colunas correspondentes de DF diferentes



final_data <- subset(final_sepsis_table, select = -c(category,`mRna - Description`))
# deleta todas informações desnecessárias da tabela final

write.csv(final_data, "clean_data.csv", row.names = FALSE)
# escreve final_data em um arquivo csv

rm(list = ls())
# limpa o environment
