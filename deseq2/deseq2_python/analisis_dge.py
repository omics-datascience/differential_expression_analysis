# analisis_dge.py

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def run_deseq_analysis(counts, metadata, group_col, group1, group2):
    """
    Ejecuta el pipeline completo de PyDESeq2 sobre los datos preprocesados.

    Returns:
        tuple: Una tupla con el objeto DeseqDataSet (dds) y el DataFrame 
               de resultados ordenado.
    """
    # Pre-filtrado: Eliminar genes con conteos totales muy bajos
    counts_to_keep = counts.sum(axis=1) >= 10
    counts_filtered = counts[counts_to_keep]

    # Crear el objeto DESeqDataSet
    dds = DeseqDataSet(
        counts=counts_filtered.T,
        metadata=metadata,
        design_factors=group_col
    )

    # Ejecutar el pipeline de DESeq2
    print("INFO - Ajustando el modelo DESeq2...")
    dds.deseq2()

    # Obtener los resultados del análisis estadístico
    print(f"INFO - Realizando contraste: {group2} vs {group1}")
    stat_res = DeseqStats(dds, contrast=(group_col, group2, group1))

    # Ejecuta el cálculo y muestra el resumen en la consola
    stat_res.summary()

    # Accede al DataFrame de resultados y lo ordena por 'padj'
    results_df = stat_res.results_df.sort_values('padj')

    return dds, results_df
