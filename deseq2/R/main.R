# INSTALACIÓN DE PAQUETES (ejecutar solo una vez si es necesario)
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }
# BiocManager::install(c("DESeq2", "SummarizedExperiment"))
# install.packages(c("tidyverse", "pheatmap"))


suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
})

# Cargar los módulos con las funciones 
source("preprocessing.R")
source("analysis_dea.R")

# Rutas a los archivos de entrada y directorio de salida
gct_file_gz <- "../../example/datasets/acc_tcga_gdc/data_mrna_seq_read_counts.txt"
metadata_file <- "../../example/datasets/acc_tcga_gdc/metadata.tsv"
output_dir <- "../../example/results_r/acc_tcga_gdc"

# Parámetros para la comparación de grupos
columna_de_grupo <- "OS_STATUS"  # Nombre de la columna en 'metadata_file' que define los grupos.
grupo1 <- "0:LIVING"              # Nombre del primer grupo (usado como control en la comparación).
grupo2 <- "1:DECEASED"           # Nombre del segundo grupo.
  

# Crear el directorio de salida si no existe
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cat("INFO - Directorio de resultados:", output_dir, "\n")

# Preprocesamiento de datos:
# Carga, valida y filtra la matriz de conteos y los metadatos.
# Retorna una lista con la matriz de conteos y los metadatos listos para DESeq2.
cat("INFO - Iniciando preprocesamiento de datos...\n")
datos_preprocesados <- preparar_datos(
  gct_path = gct_file_gz,
  metadata_path = metadata_file,
  col_grupo = columna_de_grupo,
  nombre_grupo1 = grupo1,
  nombre_grupo2 = grupo2
)

# -- PASO 2: Análisis de Expresión Diferencial --
# Ejecuta el pipeline de DESeq2.
# Retorna una lista que contiene el objeto 'dds' y la tabla de resultados.
cat("INFO - Ejecutando análisis con DESeq2...\n")
resultados_deseq <- ejecutar_analisis_deseq(
  matriz_conteos = datos_preprocesados$matriz_conteos,
  metadatos = datos_preprocesados$metadatos,
  col_grupo = columna_de_grupo,
  nombre_grupo1 = grupo1,
  nombre_grupo2 = grupo2
)

# -- PASO 3: Guardar resultados principales --
# Guardar la tabla completa de resultados del análisis en un archivo CSV.
cat("INFO - Guardando tabla de resultados...\n")
ruta_resultados_csv <- file.path(output_dir, sprintf("resultados_completos_%s_vs_%s.csv", grupo2, grupo1))
utils::write.csv(as.data.frame(resultados_deseq$resultados), file = ruta_resultados_csv)

cat("Análisis completado! Resultados en la carpeta:", output_dir, "\n")