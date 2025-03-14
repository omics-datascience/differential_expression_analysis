process_cbioportal_metadata <- function(clinical_path, sample_path) {
  cat("Loading cBioPortal metadata dataset...", "\n")

  # Load files
  clinical_data <- read.table(clinical_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  sample_data <- read.table(sample_path, header = TRUE, sep = "\t", check.names = FALSE)

  ##### Process clinical data #####
  clinical_data <- clinical_data[, c("PATIENT_ID", "SEX")]
  clinical_data$SEX <- toupper(as.character(clinical_data$SEX))
  clinical_data <- clinical_data[!duplicated(clinical_data), ]  # Remove duplicates

  ##### Process sample data #####
  sample_data <- sample_data[, c("SAMPLE_ID", "PATIENT_ID")]

  ##### Merge clinical & sample data #####
  metadata <- merge(sample_data, clinical_data, by = "PATIENT_ID")
  metadata <- metadata[, c("SAMPLE_ID", "SEX")]
  metadata$SAMPLE_ID <- gsub("-", ".", metadata$SAMPLE_ID)  # Replace hyphens with dots

  # Count duplicates
  n_before <- nrow(metadata)
  metadata <- metadata[!duplicated(metadata), ]
  n_after <- nrow(metadata)

  # Summary messages
  cat("Metadata dataset: Summary", "\n")
  cat("Metadata dataset: Rows (samples) removed due to duplication:", n_before - n_after, "\n")
  cat("Metadata dataset: Number of samples:", n_after, "\n")

  # Summary by sex
  summary <- metadata %>%
    dplyr::count(SEX) %>%
    dplyr::mutate(Proportion = 100 * n / sum(n))

  return(list(metadata = metadata, summary = summary))
}
