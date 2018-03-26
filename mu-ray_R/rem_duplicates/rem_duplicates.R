library(readr)
library(dplyr)

# Abre os cois arrquivos que contém os dados limpos, suas expreções e log destas.
cl_data <- read_csv("clean_data.csv")
log_avg <- read_csv("media do LOG da expressao dentre os controles.csv")

# Junta os dados limpos com descrição e os logs da expressão conforme o número de identificação.
data_log <- merge(cl_data, log_avg, by = "transcript_cluster_id")

# Procura todas as entradas distintas segundo a descreição, o dia, o paciente, a expressão média
data_log_no_dupl <- distinct(data_log, description, Dia, Paciente, average_expression, .keep_all = TRUE)
#View(data_log_no_dupl)

# Seleciona todas as colunas exceto o log
clean_data_no_dupl <- select(data_log_no_dupl, 1:7)
#View(data_log_no_dupl)
write.csv(clean_data_no_dupl, "clean_data_without_duplicates.csv", row.names = FALSE)

# Encontra as diferenças entre os dados iniciais e os dados sem duplicados, resultando numa lista dos
# duplicados (tudo está computanto, nenhum dado foi perdido)
duplicated_data <- setdiff(cl_data, data_log_no_dupl)
#View(duplicated_data)
write.csv(duplicated_data, "duplicated_data.csv", row.names = FALSE)
