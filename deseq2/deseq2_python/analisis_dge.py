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
    # Pre-filtrado: Eliminar genes con conteos totales muy bajos.
    # Suma los conteos de cada gen a lo largo de todas las muestras y luego se eliminan
    # los genes con muy baja expresión (<10 en todas las muestras), ya que no aportan poder 
    # estadístico al análisis y pueden introducir ruido.
    counts_to_keep = counts.sum(axis=1) >= 10
    counts_filtered = counts[counts_to_keep]

    # Crear el objeto DESeqDataSet.
    # Este objeto contiene todos los datos y la configuración del experimento.
    # counts=counts_filtered.T: Pasa la matriz de conteos filtrada. Importante el uso de .T,
    # que transpone la matriz. PyDESeq2 espera que las muestras estén en las filas y los 
    # genes en las columnas.
    # metadata=metadata: Proporciona el DataFrame de metadatos, que describe cada muestra.
    # design_factors=group_col: Especifica cuál columna de los metadatos se usará para el 
    # modelo estadístico. Es la "fórmula del diseño" que le dice a DESeq2 qué grupos comparar.
    dds = DeseqDataSet(
        counts=counts_filtered.T,
        metadata=metadata,
        design_factors=group_col
    )

    # Ejecutar el pipeline de DESeq2
    print("INFO - Ajustando el modelo DESeq2...")
    # Ejecuta la secuencia de pasos principal de DESeq2, que internamente realiza tres calculos:
    # Size Factor Estimation: Normaliza los conteos para corregir las diferencias en la 
    # profundidad de secuenciación entre muestras.
    # Dispersion Estimation: Calcula la variabilidad de los conteos para cada gen, un paso 
    # importante para la estadística del modelo.
    # Model Fitting and Testing: Ajusta un modelo lineal generalizado binomial negativo a 
    # los datos y realiza las pruebas estadísticas.
    dds.deseq2()

    print(f"INFO - Realizando contraste: {group2} vs {group1}")
    # Obtener los resultados del análisis estadístico.add(stat_res = DeseqStats(...): Crea un objeto para calcular y almacenar los resultados de una comparación específica (contraste).
    # dds: Se le pasa el objeto dds que ya contiene el modelo ajustado.
    # contrast=('columna', 'grupo2', 'grupo1'): Este es el argumento más importante. Define el
    # contraste exacto: dentro de la columna group_col, compara el nivel group2 contra el 
    # nivel group1. Esto determina la dirección de los resultados (un log2FoldChange positivo
    # significará que el gen está más expresado en group2 que en group1)
    stat_res = DeseqStats(dds, contrast=(group_col, group2, group1))

    # Ejecuta el cálculo y muestra el resumen en la consola
    stat_res.summary()

    # Accede al DataFrame de resultados y lo ordena por 'padj'
    results_df = stat_res.results_df.sort_values('padj')

    return dds, results_df
