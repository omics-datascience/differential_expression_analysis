# differential_expression_analysis

This project is based on differential expression analysis with the following characteristics:

- It uses the [lung adenocarcinoma (LUAD)](https://www.cbioportal.org/study/summary?id=luad_tcga) dataset from The Cancer Genome Atlas Program (TCGA), obtained from cBioportal.
- It uses the SEX variable to compare genes and obtain the differentially expressed ones.  
- Use the [Limma library](https://bioconductor-org.translate.goog/packages/release/bioc/html/limma.html?_x_tr_sl=en&_x_tr_tl=es&_x_tr_hl=es&_x_tr_pto=tc/) of R to perform the DEA.

## Requerimientos

- R version 4.4.2
  
## Run DEA

Use the following command to perform differential expression analysis  

```R
Rscript aed_limma_cbioportal_TMP.r
```

The script downloads all necessary libraries before running the analysis.  

## Operating details

Of the dataset downloaded from cBioportal, only the files "datasets/luad_tcga/data_clinical_patient.txt", "datasets/luad_tcga/data_clinical_sample.txt", and "datasets/luad_tcga/data_mrna_seq_v2_rsem.txt" are used for the analysis. These correspond to the patient clinical data, the sample data, and the gene expression data, respectively.  

Expression data: The dataset has RNASeq data expressed in Transcripts per million (TPM) values ​​using [RSEM](https://github.com/deweylab/RSEM).  

### Processing of clinical and sample data

- Filter only the PATIENT_ID and SEX columns in clinical data.
- Convert SEX values ​​to uppercase in clinical data.
- Remove duplicate records in clinical data.
- Filter only the SAMPLE_ID and PATIENT_ID columns in sample data.
- Perform a merge between sample_data and clinical_data using PATIENT_ID as the key.
- Keep only the SAMPLE_ID and SEX columns.
- Replace hyphens (-) with periods (.) in SAMPLE_ID.
- Ensure there are no duplicates by removing them from the dataset (there should be no duplicates at this point).
- Count the number of samples per sex, and calculate the percentage of each sex in the sample.
- Store the entire processed dataset in a dataframe called 'metadata'

### RNASeq dataset processing

- Genes with duplicate names (Hugo_Symbol) are identified.
- Duplicate genes are removed, keeping only the first occurrence of each Hugo_Symbol.
- The 'Entrez_Gene_Id' column is removed.
- Rows where Hugo_Symbol is NA are removed.
- Hugo_Symbol is assigned as row names.
- Samples in metadata that are not present in rnaseq_tpm_data are removed.
- The result is returned in a data frame named 'rnaseq_tpm_data'.

### DEA with limma

- Create boxplots to compare expression distributions between samples.
- Convert TPM values ​​to a log2 scale (TPM + 1) to stabilize variance.
- Generate histograms of the expression distribution before and after the log2 transformation.
- Save the plots to a plots.pdf file.
- Create a design matrix with the variable SEX as a factor.
- Calculate the variance of each gene in the log-transformed data.
- Eliminate genes with a variance less than the threshold (1e-4 by default).
- Calculate the mean expression of each gene, set a percentile-based threshold (15% by default), and then eliminate genes with mean expression below this threshold.
- Use the limma package to identify genes differentially expressed by sex (SEX), using the dataset modified with the two points above.
- Fit a linear model (lmFit) and apply the eBayes adjustment.
- Extracts the results with corrected p-values ​​(Benjamini-Hochberg correction).
- Returns the table of differentially expressed genes.
- Sorts the results of the differential analysis by the adjusted p-value (adj.P.Val), in ascending order (lowest to highest).
- Sorts the results by the fold change (logFC), in descending order (highest to lowest).
- Sorts the results in two steps: first by the adjusted p-value in ascending order and then by logFC in descending order.
- Extracts the 50 genes with the lowest adjusted p-values ​​(the most significant).
- Generates a "Volcano Plot" to visualize the relationship between logFC and -log10 (adjusted p-value). Distinguishes significant genes (in red) from non-significant ones (in black).
- Saves the results in two CSV files (one with all results and one with the top 50)