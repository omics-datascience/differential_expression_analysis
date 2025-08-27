#' @title Preparar Datos para DESeq2
#' @description Carga, filtra y formatea la matriz de conteos y los metadatos.
#' @param gct_path Ruta al archivo de conteos en formato GCT (comprimido en .gz).
#'               - Formato esperado: Archivo de texto con genes en filas y muestras en columnas.
#'                 Las primeras 2 líneas son cabeceras GCT. La primera columna debe ser 'Name' (IDs de genes).
#' @param metadata_path Ruta al archivo de metadatos.
#'                - Formato esperado: Archivo de texto delimitado por tabulaciones.
#'                  La primera columna debe contener los IDs de las muestras (y ser el nombre de la columna).
#' @param col_grupo Nombre de la columna en los metadatos que contiene la información de los grupos a comparar.
#' @param nombre_grupo1 String con el nombre del primer grupo de comparación. Es el grupo de referencia (el control)
#' @param nombre_grupo2 String con el nombre del segundo grupo de comparación.
#' @return Una lista con dos elementos: 'matriz_conteos' (la matriz de conteos filtrada y ordenada)
#'         y 'metadatos' (el data frame de metadatos filtrado y ordenado).

preparar_datos <- function(gct_path, metadata_path, col_grupo, nombre_grupo1, nombre_grupo2) {

  # --- Carga y filtrado de metadatos ---

  # Cargar el archivo completo de metadatos. Se asume que la primera columna contiene los IDs de muestra.
  all_metadata <- utils::read.delim(metadata_path, row.names = 1)

  # Validar que la columna y los grupos especificados existan.
  if (!col_grupo %in% colnames(all_metadata)) {
    stop(sprintf("ERROR: La columna '%s' no existe en el archivo de metadatos.", col_grupo))
  }
  niveles_grupo <- unique(all_metadata[[col_grupo]])
  if (!nombre_grupo1 %in% niveles_grupo) {
    stop(sprintf("ERROR: El grupo '%s' no se encuentra en la columna '%s'.", nombre_grupo1, col_grupo))
  }
  if (!nombre_grupo2 %in% niveles_grupo) {
    stop(sprintf("ERROR: El grupo '%s' no se encuentra en la columna '%s'.", nombre_grupo2, col_grupo))
  }

  # Filtrar metadatos para conservar únicamente las muestras de los grupos de interés.
  metadata_filtrado <- all_metadata %>%
    dplyr::filter(.data[[col_grupo]] %in% c(nombre_grupo1, nombre_grupo2))

  cat(sprintf("INFO - Muestras encontradas: %d para '%s' y %d para '%s'.\n",
              sum(metadata_filtrado[[col_grupo]] == nombre_grupo1), nombre_grupo1,
              sum(metadata_filtrado[[col_grupo]] == nombre_grupo2), nombre_grupo2))


  # --- Carga y filtrado de la matriz de conteos ---

  # Cargar el archivo GCT de conteos. Se saltan las 2 primeras líneas de metadatos del formato.
  counts_tibble <- readr::read_tsv(gct_path, skip = 2, show_col_types = FALSE)

  # Convertir el tibble a una matriz numérica.
  # DESeq2 requiere una matriz con genes en filas y muestras en columnas.
  count_matrix <- counts_tibble %>%
    dplyr::select(-Description) %>%                # La columna 'Description' no es necesaria.
    tibble::column_to_rownames("Name") %>%         # La columna 'Name' se convierte en los nombres de las filas (genes).
    as.matrix()

  # Filtrar la matriz de conteos para quedarnos solo con las columnas (muestras)
  # que están presentes en nuestros metadatos filtrados.
  muestras_a_mantener <- rownames(metadata_filtrado)
  count_matrix_filtrada <- count_matrix[, intersect(muestras_a_mantener, colnames(count_matrix))]

  # Sincronizar el orden: Reordenar los metadatos para que el orden de sus filas (muestras)
  # coincida exactamente con el orden de las columnas (muestras) en la matriz de conteos.
  # Esto es un requisito estricto de DESeq2.
  metadata_final <- metadata_filtrado[colnames(count_matrix_filtrada), ]

  # Convierte la columna a factor
  metadata_final[[col_grupo]] <- factor(metadata_final[[col_grupo]])

  # Establece explícitamente cuál es el grupo de referencia (el control)
  # 'grupo1' debería ser tu grupo control o base.
  metadata_final[[col_grupo]] <- relevel(metadata_final[[col_grupo]], ref = grupo1)

  # Devolver los datos listos para el análisis
  return(list(matriz_conteos = count_matrix_filtrada, metadatos = metadata_final))
}