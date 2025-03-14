# Perform Differential Expression Analysis with limma
perform_differential_expression <- function(log_data_filtered, metadata) {
  library(limma)
  
  # Create a design matrix based on the metadata group (SEX in this case)
  group <- factor(metadata$SEX)
  design <- model.matrix(~ group)
  
  # Perform differential expression analysis with limma
  fit <- lmFit(log_data_filtered, design)
  ebayes_fit <- eBayes(fit)
  
  # Get the results (adjusted p-values using Benjamini-Hochberg correction)
  results <- topTable(ebayes_fit, coef = 2, adjust.method = "BH", number = Inf)
  
  return(results)
}
