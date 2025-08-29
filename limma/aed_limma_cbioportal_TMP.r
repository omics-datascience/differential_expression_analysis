library(limma)
library(edgeR)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]         # ej: "luad_tcga_pub"
atributo <- args[2]        # ej: "SEX"
if (length(args) < 2) {
  stop("ERROR: Faltan argumentos. Se necesitan dos argumentos: 'dataset' y 'atributo'.", call. = FALSE)
}

if (!dir.exists("output")) {
  dir.create("output")
} 

# Construir paths a partir del dataset
clinical_path <- file.path("../datasets", dataset, "data_clinical_patient.txt")
sample_path <- file.path("../datasets", dataset, "data_clinical_sample.txt")
tpm_rnaseq_path <- file.path("../datasets", dataset, "data_mrna_seq_v2_rsem.txt")


#### Clinical metadata processing and loading ####
source("data_processing/process_cbioportal_metadata.R")
result <- process_cbioportal_metadata(clinical_path, sample_path, atributo)
# Store processed data in variables
metadata <- result$metadata
summary <- result$summary


#### RNASeq data processing and loading ####
source("data_processing/process_cbioportal_rnaseq_tpm_data.R")
result <- process_rnaseq_tpm_data(tpm_rnaseq_path, metadata)
# Store processed data in variables
rnaseq_tpm_data <- result$rnaseq_tpm_data
metadata <- result$metadata


#### Differential expression analysis with limma ####

# Perform the transformations and analysis
source("data_processing/datasets_processing.R")
source("differential_expression_analysis/dea_with_limma.r")

# Transform data to log2 scale
log_data <- transform_to_log2(rnaseq_tpm_data)

# Filter out lowly expressed genes
log_data_filtered <- filter_lowly_expressed_genes(log_data)

# Filter genes based on variance
log_data_filtered <- filter_genes_by_variance(log_data_filtered)

# Create plots after filtering
create_plots_for_expresion_data(rnaseq_tpm_data, log_data, log_data_filtered, paste0("output/", dataset, "_", atributo, "_plots.pdf"))

# Perform differential expression analysis
all_results <- perform_differential_expression(log_data_filtered, metadata, atributo)


##### Process results #####
source("process_results/process_limma_results.r")

ordered_results_pvalue <- order_results_by_pvalue(all_results)
# ordered_results_logFC <- order_results_by_logFC(all_results)
# ordered_results_both <- order_results_by_both(all_results)

# Get the top 50 results by adjusted p-value
top_50_pvalue <- get_top_50_by_pvalue(ordered_results_pvalue)

# Write the results to CSV files
write.csv(ordered_results_pvalue, file = paste0("output/", dataset, "_", atributo, "_results.csv"), row.names = TRUE)
write.csv(top_50_pvalue, file = paste0("output/", dataset, "_", atributo, "_top50_results.csv"), row.names = TRUE)


# Creo y guardo grafico de volcano
results_plot <- create_volcano_plot(ordered_results_pvalue,
       p_value_threshold = 0.05,
       logFC_threshold = 1,
       title = paste0(dataset, ":", dataset, "Diff_Expression_By", atributo))

ggsave(filename = paste0("output/", dataset, "_", atributo, "_volcano_plot.png"),
       plot = results_plot, width = 10, height = 7)
