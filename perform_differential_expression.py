# Importaciones necesarias
import pandas as pd
import numpy as np
import warnings
from itertools import combinations
import json

# Suprimir warning específico
warnings.filterwarnings('ignore',
                       message='Environment variable "XPC_SERVICE_NAME" redefined by R',
                       category=UserWarning)

import rpy2.robjects as robjects
from rpy2.robjects.conversion import get_conversion, localconverter
from rpy2.robjects.pandas2ri import converter
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# Activar conversión automática entre pandas y R
# pandas2ri.activate() # Desactivado para evitar problemas con la conversión automática


def perform_differential_expression(log_data_filtered, metadata, clinical_attribute):
    """
    Perform Differential Expression Analysis using limma (R package) via rpy2.
    """

    # Import R packages
    limma = importr('limma')  # limma for differential expression
    stats = importr('stats')  # stats for model matrix
    base = importr('base')    # base R functions

    # 1. Validation: Ensure the clinical attribute has at least two categories
    # This is necessary because differential expression requires at least two groups to compare.
    unique_values = metadata[clinical_attribute].dropna().unique()
    if len(unique_values) < 2:
        raise ValueError(f"At least two categories are required in '{clinical_attribute}'.")

    # 2. Group vector creation
    # Converts the clinical attribute to a categorical variable (factor in R).
    # This tells the model to treat the values as groups, not as numeric values.
    group_values = metadata[clinical_attribute].astype(str).values
    group_levels = sorted(np.unique(group_values))  # Get all unique group names sorted
    group_dict = {k: i for i, k in enumerate(group_levels)}  # Map group names to indices (not strictly needed)
    group_factor = pd.Categorical(group_values, categories=group_levels)

    # 3. Conversion to R objects
    # Converts the filtered log-expression data (Pandas DataFrame) to an R matrix.
    with localconverter(get_conversion() + converter):
        r_log_data = pandas2ri.py2rpy(log_data_filtered)
    # Convert the group factor to an R factor vector
    r_group = robjects.FactorVector(group_factor)

    # 4. Design matrix construction (no intercept)
    # The design matrix encodes the group structure for the linear model.
    # Using '~ 0 + group' means no intercept: each group gets its own column.
    formula = robjects.Formula('~ 0 + group')
    env = robjects.Environment()
    env['group'] = r_group
    r_design = stats.model_matrix(formula, env)
    r_design.colnames = robjects.StrVector(group_levels)  # Set column names to group names

    # 5. Linear model fitting with limma
    # Fit the linear model to estimate mean expression for each gene in each group.
    fit = limma.lmFit(r_log_data, r_design)

    # 6. Create all possible pairwise contrasts (all-vs-all)
    # For each pair of groups, create a contrast expression like 'groupB - groupA'.
    pares = list(combinations(group_levels, 2))
    contrastes = [f"{b} - {a}" for a, b in pares]
    contrast_matrix = limma.makeContrasts(
        contrasts=robjects.StrVector(contrastes),
        levels=r_design
    )

    # 7. Apply contrasts and empirical Bayes moderation
    # Apply the contrasts to the fitted model, then use eBayes to stabilize variance estimates.
    fit2 = limma.contrasts_fit(fit, contrast_matrix)
    fit2 = limma.eBayes(fit2)

    # 8. Extract results for the first contrast
    # Get the table of differential expression results for the first contrast (logFC, p-value, adjusted p-value, etc.).
    results = limma.topTable(
        fit2,
        coef=1,  # First contrast
        number=robjects.r('Inf'),  # All genes
        adjust_method="BH"  # Benjamini-Hochberg adjustment
    )
    # Convert the R data frame to a Pandas DataFrame for further analysis in Python.
    results_df = pandas2ri.rpy2py(results)
    return results_df


def load_data(rnaseq_path, clinical_path):
    df_RNAseq = pd.read_csv(rnaseq_path)
    df_clinical = pd.read_csv(clinical_path)
    return df_RNAseq, df_clinical


def prepare_data(df_RNAseq, df_clinical, sample_column='SAMPLE_ID'):
    # If the first column contains gene names, use as index
    if 'Unnamed: 0' in df_RNAseq.columns:
        df_RNAseq = df_RNAseq.set_index('Unnamed: 0')
    # Filter to keep only common samples between RNAseq and clinical data
    if sample_column in df_clinical.columns:
        common_samples = sorted(set(df_RNAseq.columns) & set(df_clinical[sample_column]))
        df_RNAseq_filtered = df_RNAseq[common_samples]
        df_clinical_filtered = df_clinical[df_clinical[sample_column].isin(common_samples)]
        df_clinical_filtered = df_clinical_filtered.set_index(sample_column).reindex(common_samples).reset_index()
    else:
        df_RNAseq_filtered = df_RNAseq
        df_clinical_filtered = df_clinical
    return df_RNAseq_filtered, df_clinical_filtered


# Export to HTML with numbers in scientific notation
def export_html(results_fmt, output_path):
    html = results_fmt.to_html(classes='table table-striped', index=True)
    with open(output_path, 'w') as f:
        f.write("""
        <html>
        <head>
          <link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
        </head>
        <body>
        <div class=\"container\">
        <h2>Results Report</h2>
        """)
        f.write(html)
        f.write("""
        </div>
        </body>
        </html>
        """)
    print(f"Resultados exportados a {output_path}")


def main():
    rnaseq_path = "examples/RNAseq_log_TMP_acc_tcga.csv"
    clinical_path = "examples/clinical_data_acc_tcga.csv"
    clinical_attribute = "SEX"
    output_html = "examples/informe_resultados_completo.html"

    df_RNAseq, df_clinical = load_data(rnaseq_path, clinical_path)
    print("Datos cargados:")
    print(f"RNA-seq shape: {df_RNAseq.shape}")
    print(f"Clinical data shape: {df_clinical.shape}")
    print(f"Unique values in {clinical_attribute}: {df_clinical[clinical_attribute].unique()}")

    df_RNAseq_filtered, df_clinical_filtered = prepare_data(df_RNAseq, df_clinical)
    print(f"Datos filtrados - RNA-seq: {df_RNAseq_filtered.shape}, Clinical: {df_clinical_filtered.shape}")

    try:
        results = perform_differential_expression(
            df_RNAseq_filtered,
            df_clinical_filtered,
            clinical_attribute
        )
        print("Análisis completado exitosamente!")
        print(f"Resultados shape: {results.shape}")
        print("\nPrimeras filas de los resultados:")
        print(results.head())
        top_genes = results.nsmallest(10, 'adj.P.Val')
        top_genes_reset = top_genes.reset_index()
        response = {
            "status": "success",
            "message": "Differential expression analysis completed successfully.",
            "results": top_genes_reset.to_dict(orient='records')
        }
        print("\nResponse JSON:" + json.dumps(response, indent=2))
        results_fmt = format_scientific(results)
        export_html(results_fmt, output_html)
    except Exception as e:
        print(f"Error en el análisis: {e}")
        print("Verificar que R y el paquete limma estén instalados correctamente")


# Formatea todas las columnas numéricas a notación científica con 2 decimales
def format_scientific(df, decimals=4):
    df_fmt = df.copy()
    for col in df_fmt.select_dtypes(include=[np.number]).columns:
        df_fmt[col] = df_fmt[col].apply(lambda x: f"{x:.{decimals}e}")
    return df_fmt


# Función auxiliar para verificar instalación de R y limma
def check_r_installation():
    """Verificar que R y limma estén disponibles"""
    try:
        # Verificar R
        r_version = robjects.r('R.version.string')[0]
        print(f"R version: {r_version}")

        # Verificar limma
        limma = importr('limma')
        print("Limma package loaded successfully")

        return True
    except Exception as e:
        print(f"Error checking R installation: {e}")
        print("Asegúrate de tener R instalado y el paquete limma disponible")
        print("Para instalar limma en R: install.packages('BiocManager'); BiocManager::install('limma')")
        return False


# Verificar instalación
check_r_installation()

if __name__ == "__main__":
    main()
