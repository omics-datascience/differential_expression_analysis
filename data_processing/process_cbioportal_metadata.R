process_cbioportal_metadata <- function(clinical_path, sample_path, clinical_attribute) {
  library(rlang)  # Import rlang for sym function
  library(dplyr)  # Import dplyr for %>% and other functions
  cat("Loading cBioPortal metadata dataset...", "\n")

  # Load files
  clinical_data <- read.table(clinical_path, header = TRUE, skip = 4, sep = "\t", check.names = FALSE, quote = "", fill = TRUE)
  sample_data <- read.table(sample_path, header = TRUE, skip = 4, sep = "\t", check.names = FALSE, quote = "", fill = TRUE)
  ##### Process clinical data #####
  clinical_data <- clinical_data[, c("PATIENT_ID", clinical_attribute)]
  # Eliminar filas con NA en el atributo clínico
  clinical_data <- clinical_data[!is.na(clinical_data[[clinical_attribute]]), ]

  # Si la columna es de tipo texto (character o factor), limpiar y convertir a mayúsculas
  column_values <- clinical_data[[clinical_attribute]]
  if (is.character(column_values) || is.factor(column_values)) {
    # Eliminar espacios al inicio y al final, luego convertir a mayúsculas
    cleaned_values <- trimws(as.character(column_values))
    clinical_data[[clinical_attribute]] <- toupper(cleaned_values)
  }

  # Eliminar filas duplicadas
  clinical_data <- clinical_data[!duplicated(clinical_data), ]

  ##### Process sample data #####
  sample_data <- sample_data[, c("SAMPLE_ID", "PATIENT_ID")]

  ##### Merge clinical & sample data #####
  metadata <- merge(sample_data, clinical_data, by = "PATIENT_ID")
  metadata <- metadata[, c("SAMPLE_ID", clinical_attribute)]
  metadata$SAMPLE_ID <- gsub("-", ".", metadata$SAMPLE_ID)  # Replace hyphens with dots

  # Count duplicates
  n_before <- nrow(metadata)
  metadata <- metadata[!duplicated(metadata), ]
  n_after <- nrow(metadata)

  # Summary messages
  cat("Metadata dataset: Summary", "\n")
  cat("Metadata dataset: Rows (samples) removed due to duplication:", n_before - n_after, "\n")
  cat("Metadata dataset: Number of samples:", n_after, "\n")

  # Summary by clinical_attribute
 summary <- metadata %>%
  dplyr::count(!!sym(clinical_attribute)) %>%
  dplyr::mutate(Proportion = 100 * n / sum(n))

  return(list(metadata = metadata, summary = summary))
}
