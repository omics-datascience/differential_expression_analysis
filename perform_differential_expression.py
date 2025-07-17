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
from matplotlib import pyplot as plt

# Activar conversión automática entre pandas y R
# pandas2ri.activate() # Desactivado para evitar problemas con la conversión automática

# Nueva funcionalidad procesamiento de datos clinicos
def process_cbioportal_metadata(clinical_path, sample_path, clinical_attribute):
    """
    Procesar metadatos de cBioPortal - versión Python de la función R.

    Args:
        clinical_path: Ruta al archivo de datos clínicos
        sample_path: Ruta al archivo de datos de muestra
        clinical_attribute: Atributo clínico a procesar

    Returns:
        dict: Diccionario con metadata procesada y resumen
    """
    print("Cargando dataset de metadatos de cBioPortal...")

    # Cargar archivos (skip=4 para saltar las primeras 4 líneas)
    clinical_data = pd.read_csv(clinical_path, sep='\t', skiprows=4)
    sample_data = pd.read_csv(sample_path, sep='\t', skiprows=4)

    ##### Procesar datos clínicos #####
    # Filtrar solo las columnas PATIENT_ID y el atributo clínico
    clinical_data = clinical_data[['PATIENT_ID', clinical_attribute]]

    # Eliminar filas con NA en el atributo clínico
    clinical_data = clinical_data.dropna(subset=[clinical_attribute])

    # Si la columna es de tipo texto, limpiar y convertir a mayúsculas
    if clinical_data[clinical_attribute].dtype == 'object':
        # Eliminar espacios al inicio y al final, luego convertir a mayúsculas
        clinical_data[clinical_attribute] = clinical_data[clinical_attribute].astype(str).str.strip().str.upper()

    # Eliminar filas duplicadas
    clinical_data = clinical_data.drop_duplicates()

    ##### Procesar datos de muestra #####
    sample_data = sample_data[['SAMPLE_ID', 'PATIENT_ID']]

    ##### Fusionar datos clínicos y de muestra #####
    metadata = pd.merge(sample_data, clinical_data, on='PATIENT_ID', how='inner')
    metadata = metadata[['SAMPLE_ID', clinical_attribute]]

    # Reemplazar guiones con puntos en SAMPLE_ID
    metadata['SAMPLE_ID'] = metadata['SAMPLE_ID'].str.replace('-', '.')

    # Contar duplicados
    n_before = len(metadata)
    metadata = metadata.drop_duplicates()
    n_after = len(metadata)

    # Mensajes de resumen
    print("Dataset de metadatos: Resumen")
    print(f"Dataset de metadatos: Filas (muestras) eliminadas por duplicación: {n_before - n_after}")
    print(f"Dataset de metadatos: Número de muestras: {n_after}")

    # Resumen por atributo clínico - CORRECCIÓN AQUÍ
    summary = (metadata[clinical_attribute]
               .value_counts()
               .reset_index())
    # Corregir los nombres de las columnas después de value_counts()
    summary.columns = [clinical_attribute, 'n']
    # Asegurarse de que 'n' sea numérico
    summary['n'] = pd.to_numeric(summary['n'], errors='coerce')
    summary['Proportion'] = 100 * summary['n'] / summary['n'].sum()

    return {
        'metadata': metadata,
        'summary': summary
    }


# Nueva funcionalidad para procesar los datos de RNA-seq
def process_rnaseq_tpm_data(rnaseq_path, metadata):
    """
    Procesar datos de RNA-Seq TPM - versión Python de la función R.

    Args:
        rnaseq_path: Ruta al archivo de datos de RNA-Seq TPM
        metadata: DataFrame con metadatos que contiene columna SAMPLE_ID

    Returns:
        dict: Diccionario con rnaseq_tpm_data procesado y metadata filtrado
    """
    print("Cargando dataset de RNA-Seq TPM...")

    # Cargar datos de RNA-Seq
    rnaseq_tpm_data = pd.read_csv(rnaseq_path, sep='\t')

    # SOLUCIÓN: Reemplazar guiones por puntos en los nombres de las columnas sino no puede machear con SAMPLE_ID
    columns_to_rename = {}
    for col in rnaseq_tpm_data.columns:
        if col not in ['Hugo_Symbol', 'Entrez_Gene_Id']:
            columns_to_rename[col] = col.replace('-', '.')

    rnaseq_tpm_data = rnaseq_tpm_data.rename(columns=columns_to_rename)

    # Identificar y eliminar Hugo_Symbols duplicados
    duplicates = rnaseq_tpm_data[rnaseq_tpm_data['Hugo_Symbol'].duplicated()]['Hugo_Symbol']
    num_duplicates = len(duplicates)


    # Eliminar duplicados manteniendo la primera ocurrencia
    rnaseq_tpm_data = rnaseq_tpm_data.drop_duplicates(subset=['Hugo_Symbol'], keep='first')

    # Eliminar la columna Entrez_Gene_Id si existe
    if 'Entrez_Gene_Id' in rnaseq_tpm_data.columns:
        rnaseq_tpm_data = rnaseq_tpm_data.drop('Entrez_Gene_Id', axis=1)

    print(f"Dataset de RNA-Seq: Número de genes (filas) eliminados: {num_duplicates}")
    if num_duplicates > 0:
        print(f"Dataset de RNA-Seq: Hugo_Symbols duplicados: {duplicates.unique()}")

    # Obtener valores de SAMPLE_ID del DataFrame de metadatos
    valid_columns = list(set(rnaseq_tpm_data.columns) & set(metadata['SAMPLE_ID']))

    # Filtrar columnas de datos de RNA-Seq que coinciden con SAMPLE_ID
    total_samples_before = len(rnaseq_tpm_data.columns) - 1  # Restar 1 porque "Hugo_Symbol" no es una muestra
    rnaseq_tpm_data = rnaseq_tpm_data[['Hugo_Symbol'] + valid_columns]
    total_samples_after = len(rnaseq_tpm_data.columns) - 1

    print(f"Dataset de RNA-Seq: Número de muestras antes del filtrado: {total_samples_before}")
    print(f"Dataset de RNA-Seq: Número de muestras después del filtrado: {total_samples_after}")
    print(f"Dataset de RNA-Seq: Número de muestras eliminadas: {total_samples_before - total_samples_after}")

    print("Dataset de RNA-Seq: Resumen Final")
    print(f"Dataset de RNA-Seq: Número de genes: {len(rnaseq_tpm_data)}")
    print(f"Dataset de RNA-Seq: Número de muestras: {len(rnaseq_tpm_data.columns) - 1}")

    # Eliminar filas con NA en los datos de expresión génica
    # Eliminar filas donde Hugo_Symbol es NA
    rnaseq_tpm_data = rnaseq_tpm_data.dropna(subset=['Hugo_Symbol'])

    # Establecer Hugo_Symbol como índice y eliminar la columna
    rnaseq_tpm_data = rnaseq_tpm_data.set_index('Hugo_Symbol')

    # Eliminar las muestras de metadata que no están en el dataset de RNA-Seq
    metadata_filtered = metadata[metadata['SAMPLE_ID'].isin(valid_columns)]

    # Imprimir el resumen de los datos procesados de RNA-Seq
    print("Resumen de Datos de RNA-Seq Procesados:")
    print(f"Número de genes: {len(rnaseq_tpm_data)}")
    print(f"Número de muestras: {len(rnaseq_tpm_data.columns)}")

    return {
        'rnaseq_tpm_data': rnaseq_tpm_data,
        'metadata': metadata_filtered
    }



# Vamos a poner todas las funciones de datasets_processing
def transform_to_log2(rnaseq_tpm_data):
    """
    Transform RNA-Seq TPM data to log2 scale.

    Parameters:
    rnaseq_tpm_data: pandas DataFrame or numpy array
        RNA-Seq TPM data to be transformed

    Returns:
    pandas DataFrame or numpy array (same type as input)
        Log2-transformed data
    """
    # Check for negative values
    if isinstance(rnaseq_tpm_data, pd.DataFrame):
        has_negative = (rnaseq_tpm_data < 0).any().any()
    else:
        has_negative = np.any(rnaseq_tpm_data < 0)

    if has_negative:
        warnings.warn("Input data contains negative values. Log transformation may not be appropriate.")

    # Apply log2 transformation
    log_data = np.log2(rnaseq_tpm_data + 1)

    # Count NA values
    if isinstance(log_data, pd.DataFrame):
        na_count = log_data.isna().sum().sum()
    else:
        na_count = np.isnan(log_data).sum()

    print(f"Number of NA values in log-transformed data: {na_count}")

    return log_data

def filter_lowly_expressed_genes(log_data, threshold_percentile=0.15):
    """
    Filter out lowly expressed genes based on mean expression percentile.

    Parameters:
    log_data: pandas DataFrame
        Log-transformed gene expression data (genes in rows, samples in columns)
    threshold_percentile: float, default=0.15
        Percentile threshold for filtering (0.15 = 15th percentile)

    Returns:
    pandas DataFrame
        Filtered gene expression data with lowly expressed genes removed
    """
    # Calculate the average expression for each gene (across all samples)
    avg_expression = log_data.mean(axis=1)

    # Calculate the threshold based on the percentile
    threshold = np.percentile(avg_expression, threshold_percentile * 100)

    # Filter genes with average expression above or equal to the threshold
    genes_to_keep = avg_expression >= threshold
    log_data_filtered = log_data[genes_to_keep]

    # Count genes removed
    genes_removed = len(log_data) - len(log_data_filtered)

    print(f"Number of genes removed due to low expression: {genes_removed}")

    return log_data_filtered

def filter_genes_by_variance(log_data, threshold=1e-4):
    """
    Filter genes based on variance - remove genes with low variance.

    Parameters:
    log_data: pandas DataFrame
        Log-transformed gene expression data (genes in rows, samples in columns)
    threshold: float, default=1e-4
        Variance threshold - genes with variance below this will be removed

    Returns:
    pandas DataFrame
        Filtered gene expression data with low-variance genes removed
    """
    # Calculate variance for each gene (row)
    variances = log_data.var(axis=1)

    # Display summary of variances
    print("Summary of variances:")
    print(f"Min: {variances.min():.6e}")
    print(f"1st Qu.: {variances.quantile(0.25):.6e}")
    print(f"Median: {variances.median():.6e}")
    print(f"Mean: {variances.mean():.6e}")
    print(f"3rd Qu.: {variances.quantile(0.75):.6e}")
    print(f"Max: {variances.max():.6e}")
    print(f"Variance threshold: {threshold:.6e}")

    # Filter out genes with variance below the threshold
    genes_to_keep = variances > threshold
    log_data_filtered = log_data[genes_to_keep]

    # Count genes removed
    genes_removed = len(log_data) - len(log_data_filtered)

    print(f"Number of genes removed due to low variance: {genes_removed}")

    return log_data_filtered

def create_plots_for_expression_data(rnaseq_tpm_data, log_data, log_data_filtered, file_path, max_samples=20):
    """
    Crear gráficos para análisis de datos de expresión y guardar en PDF.

    Parámetros:
    rnaseq_tpm_data: pandas DataFrame
        Datos originales RNA-seq TPM (genes en filas, muestras en columnas)
    log_data: pandas DataFrame
        Datos de expresión transformados a log2
    log_data_filtered: pandas DataFrame
        Datos de expresión log2 filtrados
    file_path: str
        Ruta para guardar el archivo PDF con los gráficos
    max_samples: int, default=20
        Número máximo de muestras a incluir en los gráficos para mejorar rendimiento
    """
    # Configurar estilo de gráficos
    plt.style.use('default')
    sns.set_palette("husl")

    # Limitar número de muestras para los gráficos (mejora rendimiento)
    sample_cols_tpm = rnaseq_tpm_data.columns[:max_samples] if len(rnaseq_tpm_data.columns) > max_samples else rnaseq_tpm_data.columns
    sample_cols_log = log_data.columns[:max_samples] if len(log_data.columns) > max_samples else log_data.columns
    sample_cols_filtered = log_data_filtered.columns[:max_samples] if len(log_data_filtered.columns) > max_samples else log_data_filtered.columns

    # Crear PDF para guardar todos los gráficos
    with PdfPages(file_path) as pdf:
        # Gráfico 1: Boxplot de datos TPM originales
        print("Creando gráfico 1: Boxplot de datos TPM originales")
        fig1 = plt.figure(figsize=(12, 6))
        sns.boxplot(data=rnaseq_tpm_data[sample_cols_tpm].T)
        plt.title('Distribución de valores TPM originales (por muestra)')
        plt.xlabel('Muestras')
        plt.ylabel('Expresión TPM')
        plt.xticks(rotation=45)
        plt.tight_layout()
        pdf.savefig(fig1)
        plt.close(fig1)

        # Gráfico 2: Histograma de datos TPM originales
        print("Creando gráfico 2: Histograma de datos TPM originales")
        fig2 = plt.figure(figsize=(10, 6))
        # Tomar solo una muestra de los datos para el histograma
        flat_data = rnaseq_tpm_data[sample_cols_tpm].values.flatten()
        flat_data = flat_data[flat_data > 0.1]  # Filtrar valores muy bajos para mejor visualización
        if len(flat_data) > 100000:  # Limitar número de puntos
            flat_data = np.random.choice(flat_data, size=100000, replace=False)
        plt.hist(flat_data, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        plt.title('Distribución de valores TPM originales')
        plt.xlabel('Expresión TPM')
        plt.ylabel('Frecuencia')
        plt.yscale('log')
        plt.tight_layout()
        pdf.savefig(fig2)
        plt.close(fig2)

        # Gráfico 3: Boxplot de datos transformados log2
        print("Creando gráfico 3: Boxplot de datos transformados log2")
        fig3 = plt.figure(figsize=(12, 6))
        sns.boxplot(data=log_data[sample_cols_log].T)
        plt.title('Distribución de valores de expresión transformados a Log2 (por muestra)')
        plt.xlabel('Muestras')
        plt.ylabel('Log2(TPM + 1)')
        plt.xticks(rotation=45)
        plt.tight_layout()
        pdf.savefig(fig3)
        plt.close(fig3)

        # Gráfico 4: Histograma de datos log2
        print("Creando gráfico 4: Histograma de datos log2")
        fig4 = plt.figure(figsize=(10, 6))
        flat_log_data = log_data[sample_cols_log].values.flatten()
        if len(flat_log_data) > 100000:  # Limitar número de puntos
            flat_log_data = np.random.choice(flat_log_data, size=100000, replace=False)
        plt.hist(flat_log_data, bins=50, alpha=0.7, color='lightcoral', edgecolor='black')
        plt.title('Distribución de valores transformados a Log2')
        plt.xlabel('Log2(TPM + 1)')
        plt.ylabel('Frecuencia')
        plt.tight_layout()
        pdf.savefig(fig4)
        plt.close(fig4)

        # Gráfico 5: Boxplot de datos log2 filtrados
        print("Creando gráfico 5: Boxplot de datos log2 filtrados")
        fig5 = plt.figure(figsize=(12, 6))
        sns.boxplot(data=log_data_filtered[sample_cols_filtered].T)
        plt.title('Distribución de valores Log2 filtrados (por muestra)')
        plt.xlabel('Muestras')
        plt.ylabel('Log2(TPM + 1)')
        plt.xticks(rotation=45)
        plt.tight_layout()
        pdf.savefig(fig5)
        plt.close(fig5)

        # Gráfico 6: Histograma de datos log2 filtrados
        print("Creando gráfico 6: Histograma de datos log2 filtrados")
        fig6 = plt.figure(figsize=(10, 6))
        flat_filtered_data = log_data_filtered[sample_cols_filtered].values.flatten()
        if len(flat_filtered_data) > 100000:  # Limitar número de puntos
            flat_filtered_data = np.random.choice(flat_filtered_data, size=100000, replace=False)
        plt.hist(flat_filtered_data, bins=50, alpha=0.7, color='lightgreen', edgecolor='black')
        plt.title('Distribución de valores Log2 filtrados')
        plt.xlabel('Log2(TPM + 1)')
        plt.ylabel('Frecuencia')
        plt.tight_layout()
        pdf.savefig(fig6)
        plt.close(fig6)

        # Gráfico 7: Comparación en histogramas
        print("Creando gráfico 7: Comparación de datos en histogramas")
        fig7, axes = plt.subplots(1, 3, figsize=(15, 5))

        # TPM original (escala log para mejor visualización)
        sample_tpm = rnaseq_tpm_data[sample_cols_tpm].values.flatten()
        sample_tpm = sample_tpm[sample_tpm > 0]
        if len(sample_tpm) > 50000:
            sample_tpm = np.random.choice(sample_tpm, size=50000, replace=False)
        axes[0].hist(sample_tpm, bins=50, alpha=0.7, color='skyblue')
        axes[0].set_title('TPM Original')
        axes[0].set_xlabel('TPM')
        axes[0].set_ylabel('Frecuencia')
        axes[0].set_yscale('log')

        # Log2 transformado
        sample_log = log_data[sample_cols_log].values.flatten()
        if len(sample_log) > 50000:
            sample_log = np.random.choice(sample_log, size=50000, replace=False)
        axes[1].hist(sample_log, bins=50, alpha=0.7, color='lightcoral')
        axes[1].set_title('Log2 Transformado')
        axes[1].set_xlabel('Log2(TPM + 1)')
        axes[1].set_ylabel('Frecuencia')

        # Log2 filtrado
        sample_filtered = log_data_filtered[sample_cols_filtered].values.flatten()
        if len(sample_filtered) > 50000:
            sample_filtered = np.random.choice(sample_filtered, size=50000, replace=False)
        axes[2].hist(sample_filtered, bins=50, alpha=0.7, color='lightgreen')
        axes[2].set_title('Log2 Filtrado')
        axes[2].set_xlabel('Log2(TPM + 1)')
        axes[2].set_ylabel('Frecuencia')

        plt.tight_layout()
        pdf.savefig(fig7)
        plt.close(fig7)

        # Gráfico 8: Estadísticas de resumen
        print("Creando gráfico 8: Estadísticas de resumen")
        fig8, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Calcular estadísticas de resumen
        summary_stats = pd.DataFrame({
            'Dataset': ['TPM Original', 'Log2 Transformado', 'Log2 Filtrado'],
            'Número de Genes': [len(rnaseq_tpm_data), len(log_data), len(log_data_filtered)],
            'Número de Muestras': [len(rnaseq_tpm_data.columns), len(log_data.columns), len(log_data_filtered.columns)]
        })

        # Crear gráfico de barras
        ax1.bar(summary_stats['Dataset'], summary_stats['Número de Genes'],
                color=['skyblue', 'lightcoral', 'lightgreen'])
        ax1.set_title('Número de Genes en cada paso de procesamiento')
        ax1.set_ylabel('Número de Genes')
        ax1.tick_params(axis='x', rotation=45)

        ax2.bar(summary_stats['Dataset'], summary_stats['Número de Muestras'],
                color=['skyblue', 'lightcoral', 'lightgreen'])
        ax2.set_title('Número de Muestras en cada paso de procesamiento')
        ax2.set_ylabel('Número de Muestras')
        ax2.tick_params(axis='x', rotation=45)

        plt.tight_layout()
        pdf.savefig(fig8)
        plt.close(fig8)

    print(f"Gráficos guardados en: {file_path}")


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

# Formatea todas las columnas numéricas a notación científica con 2 decimales
def format_scientific(df, decimals=4):
    df_fmt = df.copy()
    for col in df_fmt.select_dtypes(include=[np.number]).columns:
        df_fmt[col] = df_fmt[col].apply(lambda x: f"{x:.{decimals}e}")
    return df_fmt

def main():
    clinical_path = "datasets/acc_tcga/data_clinical_patient.txt"
    sample_path = "datasets/acc_tcga/data_clinical_sample.txt"
    clinical_attribute = "SEX"

    process = process_cbioportal_metadata(
        clinical_path,
        sample_path,
        clinical_attribute
    )
    print("\nvamos a ver que sale......")
    print("\nMetadatos procesados:")
    print(process['metadata'].head())
    print("\nResumen de metadatos:")
    print(process['summary'])

    rnaseq_path = "datasets/acc_tcga/data_mrna_seq_v2_rsem.txt"
    metadata = process['metadata']  # Usando metadata del ejemplo anterior

    resultado = process_rnaseq_tpm_data(rnaseq_path, metadata)
    print("\nDatos de RNA-Seq procesados:")
    print(f"Shape: {resultado['rnaseq_tpm_data'].shape}")
    print(f"Primeras 5 filas y 5 columnas:")
    print(resultado['rnaseq_tpm_data'].iloc[:5, :5])
    print(f"\nMetadata filtrado: {len(resultado['metadata'])} muestras")

    # Transformar los datos de RNA-Seq a escala log2
    print("\nTransformando datos de RNA-Seq a escala log2...")
    log_data_filtered = transform_to_log2(resultado['rnaseq_tpm_data'])

    # Filtrar genes con baja expresión
    print("\nFiltrando genes con baja expresión...")
    filtered_log_data = filter_lowly_expressed_genes(log_data_filtered)

    # Filtrar genes con baja varianza
    print("\nFiltrando genes con baja varianza...")
    filtered_log_data = filter_genes_by_variance(filtered_log_data)

    # Plots
    print("\nCreando plots para los datos de expresión...")
    print("Saltando reportes...")
    # NO EJECUTAR
    # create_plots_for_expression_data(
    #     resultado['rnaseq_tpm_data'],
    #     log_data_filtered,
    #     filtered_log_data,
    #     "expression_data_analysis_2.pdf"
    # )

    df_RNAseq = filtered_log_data
    df_clinical = resultado['metadata']

    # Usar SEX como parámetro clinical_attribute
    clinical_attribute = "SEX"

    print("Datos cargados:")
    print(f"RNA-seq shape: {df_RNAseq.shape}")
    print(f"Clinical data shape: {df_clinical.shape}")
    print(f"Unique values in {clinical_attribute}: {df_clinical[clinical_attribute].unique()}")

    # Preparar datos para el análisis
    # Asegurarse de que los datos de RNA-seq estén en el formato correcto
    # (genes en filas, muestras en columnas)

    # Si la primera columna contiene nombres de genes, usar como índice
    if 'Unnamed: 0' in df_RNAseq.columns:
        df_RNAseq = df_RNAseq.set_index('Unnamed: 0')

    sample_column = 'SAMPLE_ID'
    if sample_column in df_clinical.columns:
        common_samples = sorted(set(df_RNAseq.columns) & set(df_clinical[sample_column]))
        df_RNAseq_filtered = df_RNAseq[common_samples]
        df_clinical_filtered = df_clinical[df_clinical[sample_column].isin(common_samples)]
        df_clinical_filtered = df_clinical_filtered.set_index(sample_column).reindex(common_samples).reset_index()
    else:
        df_RNAseq_filtered = df_RNAseq
        df_clinical_filtered = df_clinical

    print(f"Datos filtrados - RNA-seq: {df_RNAseq_filtered.shape}, Clinical: {df_clinical_filtered.shape}")

    # Ejecutar análisis de expresión diferencial
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

        # Mostrar genes más significativos
        # print("\nGenes más significativos (top 10):")
        top_genes = results.nsmallest(10, 'adj.P.Val')
        # print(top_genes[['logFC', 'P.Value', 'adj.P.Val']])

        top_genes_reset = top_genes.reset_index()  # El índice (nombre del gen) pasa a columna
        # API response
        response = {
            "status": "success",
            "message": "Differential expression analysis completed successfully.",
            "results": top_genes_reset.to_dict(orient='records')
        }
        print("\nResponse JSON:" + json.dumps(response, indent=2))

    except Exception as e:
        print(f"Error en el análisis: {e}")
        print("Verificar que R y el paquete limma estén instalados correctamente")

    results_fmt = format_scientific(results)

    # Exporta a HTML con los números en notación científica
    html = results_fmt.to_html(classes='table table-striped', index=True)
    with open('examples/informe_resultados_completo2.html', 'w') as f:
        f.write("""
        <html>
        <head>
          <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
        </head>
        <body>
        <div class="container">
        <h2>Informe de Resultados</h2>
        """)
        f.write(html)
        f.write("""
        </div>
        </body>
        </html>
        """)

    print("Resultados exportados a examples/informe_resultados_completo2.html")
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
