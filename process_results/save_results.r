# Function to generate a timestamp string for filenames
generate_timestamp <- function() {
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  return(timestamp)
}

# write results to a CSV file
write_results_to_csv <- function(results, file_name_prefix) {
  timestamp <- generate_timestamp()  # Get the current timestamp
  file_name <- paste0(file_name_prefix, "_", timestamp, ".csv")  # Create the full filename
  write.table(results, file = file_name, sep = ",", row.names = TRUE)
  cat("Results written to", file_name, "\n")
}

