# transform the RNA-Seq TPM data to a log2 scale
transform_to_log2 <- function(rnaseq_tpm_data) {
  log_data <- log2(rnaseq_tpm_data + 1)
  cat("Number of NA values in log-transformed data:", sum(is.na(log_data)), "\n")
  return(log_data)
}


# Create Distribution and Boxplot Graphs
create_plots_for_expresion_data <- function(rnaseq_tpm_data, log_data, log_data_filtered) {  
  pdf("plots.pdf", width = 11, height = 8.5)
  # Melt the data for plotting
  tpm_melted <- melt(rnaseq_tpm_data, variable.name = "Sample", value.name = "Expression")
  log_melted <- melt(log_data, variable.name = "Sample", value.name = "Expression")
  log_filtered_melted <- melt(log_data_filtered, variable.name = "Sample", value.name = "Expression")
  
  # Histogram before and after log transformation
  p1 <- ggplot(tpm_melted, aes(x = Expression)) +
    geom_histogram(bins = 50, fill = "skyblue", color = "black") +
    ggtitle("Distribución TPM") +
    theme_minimal()
  
  p2 <- ggplot(log_melted, aes(x = Expression)) +
    geom_histogram(bins = 50, fill = "salmon", color = "black") +
    ggtitle("Distribución log2(TPM + 1)") +
    theme_minimal()

   p3 <- ggplot(log_filtered_melted, aes(x = Expression)) +
    geom_histogram(bins = 50, fill = "salmon", color = "black") +
    ggtitle("Distribución log2(TPM + 1) Filtered") +
    theme_minimal()
  
  # Show histograms side by side
  grid.arrange(p1, p2, p3, ncol = 3)
  
  # Boxplots for comparing distributions
  p4 <- ggplot(tpm_melted, aes(x = Sample, y = Expression)) +
    geom_boxplot(fill = "skyblue") +
    ggtitle("Boxplot TPM") +
    theme_minimal()
  
  p5 <- ggplot(log_filtered_melted, aes(x = Sample, y = Expression)) +
    geom_boxplot(fill = "salmon") +
    ggtitle("Boxplot log2(TPM + 1) filtered") +
    theme_minimal()
  
  # Show boxplots side by side
  grid.arrange(p4, p5, ncol = 2)
}


# Filter Genes Based on Variance
filter_genes_by_variance <- function(log_data, threshold = 1e-4) {
  # library(matrixStats)
  
  # Calculate variance for each gene (row)
  variances <- rowVars(as.matrix(log_data))
  cat("Summary of variances:\n")
  summary(variances) # Display a summary of variances
  
  # Filter out genes with variance below the threshold
  genes_to_remove <- names(variances[variances <= threshold])
  log_data_filtered <- log_data[!(rownames(log_data) %in% genes_to_remove), ]
  
  cat("Number of genes removed due to low variance:", length(genes_to_remove), "\n")
  return(log_data_filtered)
}


# Remove low expression genes
filter_lowly_expressed_genes <- function(log_data, threshold_percentile = 0.15) {
  # Calculate the average expression for each gene
  avg_expression <- rowMeans(log_data, na.rm = TRUE)
  
  # Calculate the threshold based on the percentile
  threshold <- quantile(avg_expression, threshold_percentile)
  
  # Filter genes with average expression below the threshold
  genes_to_keep <- names(avg_expression[avg_expression >= threshold])
  log_data_filtered <- log_data[genes_to_keep, ]
  
  cat("Number of genes removed due to low expression:", nrow(log_data) - length(genes_to_keep), "\n")
  return(log_data_filtered)
}

