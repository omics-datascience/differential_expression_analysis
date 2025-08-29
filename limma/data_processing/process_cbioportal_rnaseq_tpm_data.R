process_rnaseq_tpm_data <- function(rnaseq_path, metadata) {
  cat("Loading RNA-Seq TPM dataset...", "\n")
  
  # Load RNASeq data
  rnaseq_tpm_data <- read.table(rnaseq_path, header = TRUE, sep = "\t")

  # Identify and remove duplicate Hugo_Symbols
  duplicates <- rnaseq_tpm_data %>% 
    filter(duplicated(Hugo_Symbol)) %>% 
    pull(Hugo_Symbol)
  num_duplicates <- length(duplicates)

  rnaseq_tpm_data <- rnaseq_tpm_data %>%
    distinct(Hugo_Symbol, .keep_all = TRUE) %>%  # Remove duplicates by 'Hugo_Symbol'
    select(-Entrez_Gene_Id)                      # Remove 'Entrez_Gene_Id' column

  cat("RNA-Seq Dataset: Number of genes (rows) removed:", num_duplicates, "\n")
  cat("RNA-Seq Dataset: Duplicate Hugo_Symbols:", unique(duplicates), "\n")

  # Get SAMPLE_ID values from metadata dataframe
  valid_columns <- intersect(colnames(rnaseq_tpm_data), metadata$SAMPLE_ID)
  
  # Filter columns of RNA-Seq data that match SAMPLE_ID
  total_samples_before <- ncol(rnaseq_tpm_data) - 1  # Subtract 1 because "Hugo_Symbol" is not a sample
  rnaseq_tpm_data <- rnaseq_tpm_data[, c("Hugo_Symbol", valid_columns)]
  total_samples_after <- ncol(rnaseq_tpm_data) - 1 

  cat("RNA-Seq Dataset: Number of samples before filtering:", total_samples_before, "\n")
  cat("RNA-Seq Dataset: Number of samples after filtering:", total_samples_after, "\n")
  cat("RNA-Seq Dataset: Number of samples removed:", total_samples_before - total_samples_after, "\n")

  cat("RNA-Seq Dataset: Final Summary", "\n")
  cat("RNA-Seq Dataset: Number of genes:", nrow(rnaseq_tpm_data), "\n")
  cat("RNA-Seq Dataset: Number of samples:", ncol(rnaseq_tpm_data) - 1, "\n")

  # Remove rows with NA in the gene expression data
  rows_with_na <- is.na(rnaseq_tpm_data[, 1])
  rnaseq_tpm_data <- rnaseq_tpm_data[!rows_with_na, ]

  # Set rownames to Hugo_Symbol and remove the first column
  rownames(rnaseq_tpm_data) <- rnaseq_tpm_data[, 1]
  rnaseq_tpm_data <- rnaseq_tpm_data[, -1]  # Remove the first column (now index)

  # Remove the samples from metadata that are not in the RNA-Seq dataset
  metadata <- metadata %>% dplyr::filter(SAMPLE_ID %in% valid_columns)

  # Print the summary of the processed RNA-Seq data
  cat("Processed RNA-Seq Data Summary:\n")
  cat("Number of genes:", nrow(rnaseq_tpm_data), "\n")
  cat("Number of samples:", ncol(rnaseq_tpm_data), "\n")

  return(list(rnaseq_tpm_data = rnaseq_tpm_data, metadata = metadata))
}
