library(httr)
library(jsonlite)

# Function to retrieve a list of studies from cBioPortal
get_cbioportal_studies <- function() {
  # Construct API URL to get all studies with a summary projection
  url_get_studies <- "https://www.cbioportal.org/api/studies?projection=SUMMARY&pageSize=10000000&pageNumber=0&direction=ASC" # nolint

  # Request to the API
  response <- GET(url_get_studies, add_headers(Accept = "application/json"))
  stop_for_status(response)  # Stop execution if request failed

  # Parse response as text and then convert from JSON
  data <- content(response, as = "text", encoding = "UTF-8")
  studies <- fromJSON(data)

  # Create a data frame with only selected fields from the response
  data.frame(
    studyId = studies$studyId,
    name = studies$name,
    description = studies$description,
    stringsAsFactors = FALSE
  )
}

# Function to get RNA-seq molecular profiles for a given study ID
get_molecular_profiles_rna_seq <- function(study_id) {
  # Construct API URL to get molecular profiles for the given study
  url_molecular_profiles <- paste0("https://www.cbioportal.org/api/studies/",
    study_id,
    "/molecular-profiles?projection=SUMMARY&pageSize=10000000&pageNumber=0"
  )

  # Request to the API
  response <- GET(url_molecular_profiles, add_headers(Accept = "application/json"))
  stop_for_status(response)  # Stop if request failed

  # Parse the JSON response
  data <- content(response, as = "text", encoding = "UTF-8")
  profiles <- fromJSON(data)

  # Filter profiles to include only specific RNA-seq expression types
  filtered_profiles <- profiles[
    profiles$molecularAlterationType == "MRNA_EXPRESSION" &
      profiles$name %in% c(
        "mRNA expression (RNA Seq V2 RSEM)",
        "mRNA expression (TPM)",
        "mRNA expression (FPKM)"
      ),
  ]

  # Create a data frame with relevant fields
  data.frame(
    molecularProfileId = filtered_profiles$molecularProfileId,
    molecularAlterationType = filtered_profiles$molecularAlterationType,
    datatype = filtered_profiles$datatype,
    name = filtered_profiles$name,
    description = filtered_profiles$description,
    studyId = filtered_profiles$studyId,
    stringsAsFactors = FALSE
  )
}

# Function to find studies that have RNA-seq data and save their profiles to a file
get_studies_with_rna_seq <- function() {
  studies <- get_cbioportal_studies()  # Retrieve all studies
  studies_with_rna_seq <- c()  # Initialize empty vector to store matching study IDs

  # Create an empty TSV file with column headers
  column_names <- c("molecularProfileId", "molecularAlterationType", "datatype", "name", "description", "studyId")
  empty_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
  colnames(empty_df) <- column_names
  write.table(empty_df, file = "datasets_with_rna_seq.tsv", sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)

  # Loop through each study to check for RNA-seq profiles
  for (study_id in studies$studyId) {
    message("Checking study: ", study_id)
    tryCatch({
      profiles <- get_molecular_profiles_rna_seq(study_id)
      if (nrow(profiles) > 0) {
        # Append RNA-seq profiles to TSV file
        write.table(profiles, file = "cbioportal_datasets_with_rna_seq.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE) # nolint
        studies_with_rna_seq <- c(studies_with_rna_seq, study_id)  # Add correct profile to list
      }
    },
    error = function(e) {
      message("Error for study ", study_id, ": ", e$message)
    })
  }

  studies_with_rna_seq
}

# Function to download and extract dataset files for a study
download_and_extract_dataset <- function(study_id, download_dir = "datasets") {
  if (!dir.exists(download_dir)) {
    dir.create(download_dir)  # Create directory if it doesn't exist
  }

  # Construct URL and file paths
  url <- paste0("https://cbioportal-datahub.s3.amazonaws.com/", study_id, ".tar.gz")
  file_name <- paste0(study_id, ".tar.gz")
  file_path <- file.path(download_dir, file_name)

  download_yn = FALSE
  # Attempt to download the tar.gz file
  tryCatch({
    message(paste("Downloading study:", study_id))
    download.file(url, file_path, mode = "wb", quiet = TRUE)
    download_yn = TRUE
  }, error = function(e) {
    message(paste("Error downloading file:", e$message))
  })

  # If download succeeded, extract and clean up
  if (download_yn) {
    untar(file_path, exdir = download_dir)  # Extract contents
    file.remove(file_path)  # Delete the tar.gz file

    extracted_folder <- file.path(download_dir, study_id)

    if (dir.exists(extracted_folder)) {
      # Delete unnecessary subdirectories if present
      case_list_path <- file.path(extracted_folder, "case_lists")
      normals_path <- file.path(extracted_folder, "normals")
      validation_reports_path <- file.path(extracted_folder, "validation_reports")
      if (dir.exists(case_list_path)) {
        unlink(case_list_path, recursive = TRUE)
      }
      if (dir.exists(normals_path)) {
        unlink(normals_path, recursive = TRUE)
      }
      if (dir.exists(validation_reports_path)) {
        unlink(validation_reports_path, recursive = TRUE)
      }

      # Keep only a selected set of files, delete others
      files_in_folder <- list.files(extracted_folder)
      files_to_keep <- c("data_clinical_patient.txt",
        "data_clinical_sample.txt",
        "data_mrna_seq_v2_rsem.txt", "data_mrna_seq_fpkm.txt",
        "data_mrna_seq_tpm.txt", "data_mrna_seq_fpkm.txt",
        "data_mrna_seq_read_counts.txt",
        "data_mrna_seq_tpm_zscores_ref_all_samples.txt",
        "data_mrna_seq_fpkm_zscores_ref_all_samples.txt",
        "data_mrna_seq_read_counts_zscores_ref_all_samples.txt"
      )

      for (file in files_in_folder) {
        if (!(file %in% files_to_keep)) {
          file_to_delete <- file.path(extracted_folder, file)
          file.remove(file_to_delete)  # Remove unwanted file
        }
      }
    }
  }

  file_path  # Return the path to the downloaded file (even if it was deleted later)
}

# Main execution: get all RNA-seq studies and download their datasets
rna_seq_studies <- get_studies_with_rna_seq()
for (study_id in rna_seq_studies) {
  downloaded_data_dir <- download_and_extract_dataset(study_id, download_dir = "datasets")
  if (!is.null(downloaded_data_dir)) {
    message(paste("Dataset", study_id, "downloaded and extracted to", downloaded_data_dir))
  } else {
    message(paste("Download or extraction failed for study", study_id))
  }
}
