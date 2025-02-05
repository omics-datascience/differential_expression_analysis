# Función para verificar e instalar paquetes
Sys.getenv("R_LIBS_USER")
Sys.getenv("R_LIBS")

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install()

verificar_e_instalar <- function(paquete) {
  if (!require(paquete, quietly = TRUE)) {
    BiocManager::install(paquete)
  } 
}

# Instalar el paquete matrixStats si no lo tienes
if (!require("matrixStats")) install.packages("matrixStats")
if(!require(ggrepel)){install.packages("ggrepel")}

# Verificar e instalar limma
verificar_e_instalar("limma")
# Verificar e instalar edgeR (solo para graficas de volcan)
verificar_e_instalar("edgeR")

# Cargar las bibliotecas (ya deberían estar instaladas con el paso anterior)
library(limma)
library(edgeR)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(ggrepel)

#### Cargo y proceso metadatos clinicos ####
cat("Cargando dataset de metadatos...", "\n")
# Asegur que el dataset tenga los IDs de muestra en la primera columna y los grupos en la segunda.
clinical_data <- read.table("/home/mauricio/Desktop/prueba_expresion/datasets/luad_tcga/data_clinical_patient.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
sample_data <- read.table("/home/mauricio/Desktop/prueba_expresion/datasets/luad_tcga/data_clinical_sample.txt", header = TRUE, sep = "\t", check.names = FALSE)

##### Proceso clinical data #####
clinical_data <- clinical_data[, c("PATIENT_ID", "SEX")]
clinical_data$SEX <- toupper(as.character(clinical_data$SEX))
# unique(clinical_data$SEX) # chequeo que solo haya dos formas de expresar el sexo biologico
# Elimino duplicados
clinical_data <- clinical_data[!duplicated(clinical_data), ]

##### Proceso sample data #####
sample_data <- sample_data[, c("SAMPLE_ID", "PATIENT_ID")]

##### Mergeo clinical & sample data #####
# Realizar el merge entre sample_data y clinical_data usando 'PATIENT_ID'
metadatos <- merge(sample_data, clinical_data, by = "PATIENT_ID")
# Seleccionar solo las columnas 'SAMPLE_ID' y 'SEX'
metadatos <- metadatos[, c("SAMPLE_ID", "SEX")]
# Reemplazar los guiones medios por puntos en la columna 'SAMPLE_ID' (Porque asi se llamaran en el dataset de RNASeq por defecto segun R)
metadatos$SAMPLE_ID <- gsub("-", ".", metadatos$SAMPLE_ID)


# Número de muestras (filas) antes de eliminar duplicados
n_before <- nrow(metadatos)
metadatos <- metadatos[!duplicated(metadatos), ]
# Número de muestras después de eliminar duplicados
n_after <- nrow(metadatos)

cat("Dataset de metadatos: Resumen", "\n")
cat("Dataset de metadatos: Filas (muestras) eliminadas por duplicacion:", n_before - n_after, "\n")
cat("Dataset de metadatos: Numero de muestras:", n_after, "\n")
metadatos %>%
  dplyr::count(SEX) %>%
  dplyr::mutate(Proporcion = 100 * n / sum(n))
# SAMPLE_ID tiene los muestras 
# SEX tiene los grupos (FEMALE Y MALE)

#### Cargo y proceso datos de RNASeq (TPM data) ####
cat("Cargando dataset de RNASeq normalizados a TPM...", "\n")
rnaseq_tpm_data <- read.table("/home/mauricio/Desktop/prueba_expresion/datasets/luad_tcga/data_mrna_seq_v2_rsem.txt", header = TRUE, sep = "\t")

# Identificar y eliminar los Hugo_Symbols duplicados
duplicados <- rnaseq_tpm_data %>% 
  filter(duplicated(Hugo_Symbol)) %>% 
  pull(Hugo_Symbol)
num_duplicados <- length(duplicados)

rnaseq_tpm_data <- rnaseq_tpm_data %>%
  distinct(Hugo_Symbol, .keep_all = TRUE) %>%   # Elimina duplicados por 'Hugo_Symbol'
  select(-Entrez_Gene_Id)                       # Elimina la columna 'Entrez_Gene_Id'

cat("Dataset de RNASeq: Número de genes (filas) eliminados:", num_duplicados, "\n")
cat("Dataset de RNASeq: Hugo_Symbols duplicados:", unique(duplicados), "\n")

# Obtener los valores de SAMPLE_ID del dataframe metadatos
# Obtengo muestras en comun entre los dataset de datos y metadatos
valid_columns <- intersect(colnames(rnaseq_tpm_data), metadatos$SAMPLE_ID)
# Filtrar las columnas de rnaseq_tpm_data que coincidan con los SAMPLE_ID
total_muestras_antes <- ncol(rnaseq_tpm_data) - 1  # Restamos 1 porque "Hugo_Symbol" no es una muestra
rnaseq_tpm_data <- rnaseq_tpm_data[, c("Hugo_Symbol", valid_columns)]
total_muestras_despues <- ncol(rnaseq_tpm_data) - 1 

cat("Dataset de RNASeq: Número de muestras antes de la filtrar:", total_muestras_antes, "\n")
cat("Dataset de RNASeq:Número de muestras después de la filtrar:", total_muestras_despues, "\n")
cat("Dataset de RNASeq:Número de muestras eliminadas:", total_muestras_antes - total_muestras_despues, "\n")

cat("Dataset de RNASeq: Resumen Final", "\n")
cat("Dataset de RNASeq: Numero de genes:", nrow(rnaseq_tpm_data), "\n")
cat("Dataset de RNASeq: Numero de muestras:", ncol(rnaseq_tpm_data) - 1, "\n")

filas_con_na <- is.na(rnaseq_tpm_data[, 1])
# Eliminar las filas con NA 
rnaseq_tpm_data <- rnaseq_tpm_data[!filas_con_na, ]


rownames(rnaseq_tpm_data) <- rnaseq_tpm_data[, 1]
rnaseq_tpm_data <- rnaseq_tpm_data[, -1] # Elimina la primera columna de genes (ahora es el índice)

# Elimino las muestras en metadatos que finalmente no estuvieron en el dataset de RNASeq
metadatos <- metadatos %>% dplyr::filter(SAMPLE_ID %in% valid_columns)


#### Análisis de expresión diferencial con limma ####




# Transformar los datos a escala logarítmica (Asi lo exige lima cuando usamos datos en TPM)
# Separar los símbolos de genes Hugo y los datos de expresión
log_data <- log2(rnaseq_tpm_data + 1)
sum(is.na(log_data)) 

# Filtrar genes con varianza cero o cercana a cero
varianzas <- matrixStats::rowVars(as.matrix(log_data))
summary(varianzas) # Muestra un resumen de las varianzas
hist(log(varianzas)) # Histograma de las varianzas (en escala logarítmica)
umbral_varianza <- 1e-4  # umbral pequeño para varianza "cercana a cero"
genes_a_eliminar <- names(varianzas[varianzas <= umbral_varianza])
log_data <- log_data[!(rownames(log_data) %in% genes_a_eliminar), ]

# armar una matriz de diseño
group <- factor(metadatos$SEX) # Extraer las etiquetas de grupo como factor
design <- model.matrix(~ group)

# Realizar el análisis de expresión diferencial con limma
fit <- lmFit(log_data, design)
ebayes_fit <- eBayes(fit)

# Obtener los resultados
all_results <- topTable(ebayes_fit, coef = 2, adjust.method = "BH", number = Inf)


write.table(all_results, file = "resultados_limma_tpm.csv", sep = ",", row.names = TRUE)

# Ordenar por valor p ajustado (de menor a mayor)
ordered_results_pvalue <- all_results[order(all_results$adj.P.Val), ]

# Ordenar por log fold change (de mayor a menor)
ordered_results_logFC <- all_results[order(all_results$logFC, decreasing = TRUE), ]

# Ordenar por valor p ajustado y luego por log fold change
ordered_results_both <- all_results[order(all_results$adj.P.Val, decreasing = c(FALSE)), ]
top_50_pvalue <- head(ordered_results_pvalue, 50)


# Paso 7: Explorar y guardar los resultados
head(top_50_pvalue)
summary(top_50_pvalue)



# Crear el gráfico de volcán con ggplot2
ggplot(ordered_results_pvalue, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Significant", "Not Significant")), size=2) + # Colorear los puntos según el valor p y logFC
  geom_text_repel(data = subset(ordered_results_pvalue, adj.P.Val < 0.05 & abs(logFC) > 1), # Etiquetas para genes significativos
                  aes(label = rownames(subset(ordered_results_pvalue, adj.P.Val < 0.05 & abs(logFC) > 1))), # Etiquetas con nombres de genes
                  size = 3, 
                  nudge_x = 0.2, 
                  nudge_y = 0.2,
                  segment.size = 0.2,
                  max.overlaps = Inf) + # Manejo de solapamiento de etiquetas
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) + # Colores personalizados
  geom_hline(yintercept = -log10(0.05), linetype="dashed") + # Línea horizontal para el umbral de p-valor
  geom_vline(xintercept = c(-1, 1), linetype="dashed") + # Líneas verticales para el umbral de logFC
  labs(title = "Gráfico de Volcán", x = "Log Fold Change", y = "-log10(p-value ajustado)") + # Etiquetas de los ejes y título
  theme_bw() + # Estilo de tema en blanco y negro
  theme(legend.position = "none") # Ocultar la leyenda
