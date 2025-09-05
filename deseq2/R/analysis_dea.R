#' @title Ejecutar Análisis DESeq2
#' @description Realiza el análisis de expresión diferencial usando DESeq2.
#' @param matriz_conteos La matriz de conteos de genes (genes en filas, muestras en columnas).
#' @param metadatos El data frame de metadatos, con filas ordenadas para coincidir con las columnas de 'matriz_conteos'.
#' @param col_grupo El nombre de la columna en 'metadatos' que define los grupos experimentales.
#' @param nombre_grupo1 El nombre del grupo de referencia (control).
#' @param nombre_grupo2 El nombre del grupo a comparar contra la referencia.
#' @return Una lista con dos elementos: 'dds' (el objeto DESeqDataSet después de ejecutar DESeq)
#'         y 'resultados' (un data frame con los resultados ordenados por p-valor ajustado).

ejecutar_analisis_deseq <- function(matriz_conteos, metadatos, col_grupo, nombre_grupo1, nombre_grupo2) {

  # Crear el objeto DESeqDataSet a partir de la matriz de conteos y los metadatos.
  # La fórmula de diseño `~ col_grupo` indica a DESeq2 que modele los conteos en función de la variable de grupo.
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = matriz_conteos,
                                        colData = metadatos,
                                        design = stats::as.formula(paste("~", col_grupo)))

  # Pre-filtrado: Eliminar genes con conteos muy bajos (p. ej., menos de 10 conteos totales en todas las muestras).
  # Esto reduce el ruido y mejora la potencia estadística.
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep,]

  # Ejecutar el pipeline de DESeq2. Esta única función se encarga de:
  # 1. Estimar los factores de normalización (size factors).
  # 2. Estimar la dispersión de los datos.
  # 3. Ajustar un modelo lineal generalizado y realizar las pruebas de hipótesis.
  dds <- DESeq2::DESeq(dds)

  # Obtener los resultados del análisis.
  # El contraste define la comparación: `c(columna, grupo_a_comparar, grupo_base)`.
  # Se calculará el log2FoldChange como (grupo2 / grupo1).

  res <- DESeq2::results(dds, contrast = c(col_grupo, gsub(":", "_", nombre_grupo2), gsub(":", "_", nombre_grupo1)))

  # Ordenar la tabla de resultados por el p-valor ajustado (padj) para ver los genes más significativos primero.
  res_ordered <- res[order(res$padj),]

  # Devolver el objeto dds y los resultados
  return(list(dds = dds, resultados = res_ordered))
}