# Perform Differential Expression Analysis with limma
perform_differential_expression <- function(log_data_filtered, metadata, clinical_attribute) {
  library(limma)
  
  variables_count <- length(unique(metadata[[clinical_attribute]]))
  if (variables_count < 2) {
    message("El atributo clínico '", clinical_attribute, "' tiene solo una categoría: '", unique(metadata[[clinical_attribute]]), "'.")
    stop("Se necesitan al menos dos categorias para realizar análisis diferencial.")
  } 
  
  # Convertir la columna a factor le indica a R que sus valores son categorías (como grupos
  # de edad) y no números con significado matemático. Esto es clave para que el modelo trate
  # cada grupo por separado al hacer comparaciones.
  group <- factor(metadata[[clinical_attribute]])

  # Construir matriz de diseño sin intercepto (para permitir comparaciones flexibles)
  # No usamos intercepto para que R cree una columna en la matriz de diseño por cada grupo,
  # lo que permite comparar todos contra todos con contrastes explícitos. Con intercepto,
  # uno de los grupos se usa como referencia y no se modela directamente.
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  # Con intercepto seria asi: design <- model.matrix(~ group)

  # Perform differential expression analysis with limma:
  # Ajuste del modelo lineal para estimar la expresión media de cada gen en cada grupo y 
  # calcula las diferencias entre grupos. Esto nos permite identificar genes cuya expresión 
  # varía significativamente según la variable clínica de interés.
  fit <- lmFit(log_data_filtered, design)

  # Crear todos los contrastes posibles: todos contra todos
  # Para comparar cada par de grupos de forma directa (por ejemplo, grupo 2 vs grupo 1, grupo 3 vs grupo 2, etc.),
  # generamos todas las combinaciones posibles de dos grupos usando combn().
  # Luego armamos las expresiones de contraste en formato "grupoB - grupoA",
  # que le indica a limma que debe calcular la diferencia de expresión promedio entre esos dos grupos.
  grupos <- levels(group)
  pares <- combn(grupos, 2, simplify = FALSE)
  contrastes <- sapply(pares, function(par) paste0(par[2], " - ", par[1]))
  contrast_matrix <- makeContrasts(contrasts = contrastes, levels = design)

  # Ajuste de contrastes y eBayes
  # Aplicamos los contrastes definidos a los modelos ajustados para estimar, por cada gen,
  # las diferencias de expresión entre pares de grupos.
  # Luego usamos eBayes() para aplicar un ajuste bayesiano que estabiliza las estimaciones de varianza,
  # lo cual mejora la detección de diferencias significativas, especialmente en estudios con pocas muestras.
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  # Resultados por contrastes
  # Extraer la tabla de resultados para cada contraste, incluyendo estadísticas 
  # de expresión diferencial (logFC, p-valor, p-valor ajustado por método BH, etc.), 
  # ordenadas por significancia, para facilitar la interpretación y selección de 
  # genes diferencialmente expresados.
  results <- lapply(seq_len(ncol(contrast_matrix)), function(i) {
    topTable(fit2, coef = i, number = Inf, adjust.method = "BH")
  })
  names(results) <- colnames(contrast_matrix)

  return(results[[1]])
}
