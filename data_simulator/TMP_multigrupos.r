simular_tpm_multigrupos <- function(
  n_genes = 8000,
  n_muestras = 200,
  n_grupos = 4,
  n_deg = 80,
  cambios_esperados = NULL, # Vector con fold changes para cada grupo
  ruido_sd = 0,
  semilla = 42
) {
  #' Simula un dataset sintético de expresión génica TPM con múltiples grupos.
  #'
  #' Esta función genera una matriz de expresión en TPM simulada para un número arbitrario de grupos,
  #' donde un conjunto de genes (n_deg) está diferencialmente expresado de forma exclusiva en cada grupo.
  #' Para cada gen diferencialmente expresado, se aplica un fold change multiplicativo alto en el grupo asignado
  #' y un fold change bajo en los demás grupos, simulando un patrón de expresión distintivo.
  #'
  #' @param n_genes Número total de genes a simular.
  #' @param n_muestras Número total de muestras (debe ser múltiplo de n_grupos para balancear grupos).
  #' @param n_grupos Número de grupos a simular (e.g. 3 para "latino", "oriental", "americano").
  #' @param n_deg Número total de genes diferencialmente expresados (repartidos equitativamente entre grupos).
  #' @param cambios_esperados Vector numérico con fold changes esperados para cada grupo en sus genes DEA.
  #'                        Si NULL, se asume fold change 2 para todos.
  #' @param ruido_sd Desviación estándar del ruido técnico aditivo.
  #' @param semilla Semilla para reproducibilidad de la simulación.
  #'
  #' @return Lista con tres elementos:
  #' \item{tpm_df}{Data frame con expresión TPM (genes x muestras).}
  #' \item{metadata}{Data frame con muestra y grupo asignado.}
  #' \item{info_deg}{Data frame indicando para cada gen si es DEG y a qué grupo está asignado.}
  
  # Validaciones básicas de parámetros de entrada
  if (n_deg > n_genes) stop("n_deg no puede ser mayor que n_genes.")
  if (n_grupos < 2) stop("Debe haber al menos 2 grupos.")
  if (is.null(cambios_esperados)) cambios_esperados <- rep(2, n_grupos)  # Si no se pasa fold change, asignar 2 a todos
  if (length(cambios_esperados) != n_grupos) stop("Cambios esperados debe tener un valor por grupo.")
  if (n_muestras %% n_grupos != 0) stop("n_muestras debe ser múltiplo de n_grupos para balancear grupos.")

  # Fijar semilla para reproducibilidad
  set.seed(semilla)

  # Crear nombres para genes y muestras
  genes <- paste0("gen_", seq_len(n_genes))
  muestras <- paste0("muestra_", seq_len(n_muestras))

  # Asignar grupos balanceados a las muestras
  n_por_grupo <- n_muestras / n_grupos
  grupos <- rep(paste0("grupo", seq_len(n_grupos)), each = n_por_grupo)

  # Simular expresión base TPM: distribución exponencial con media ~30 TPM
  rate <- 1/30
  tpm_base <- matrix(rexp(n_genes * n_muestras, rate = rate), nrow = n_genes, ncol = n_muestras)

  # Añadir efecto biológico por muestra (variabilidad natural entre individuos)
  efecto_bio <- rnorm(n_muestras, mean = 0, sd = 0)  # ruido biológico normal, media 0, sd 3
  tpm_base <- sweep(tpm_base, 2, efecto_bio, "+")    # se suma efecto_bio a cada columna (muestra)

  # Distribuir genes DEA entre los grupos lo más equitativamente posible
  n_deg_por_grupo <- rep(floor(n_deg / n_grupos), n_grupos)
  sobrante <- n_deg - sum(n_deg_por_grupo)            # genes restantes por repartir
  if (sobrante > 0) {
    n_deg_por_grupo[1:sobrante] <- n_deg_por_grupo[1:sobrante] + 1
  }

  # Elegir aleatoriamente qué genes serán DEA
  indices_deg <- sample(n_genes, n_deg, replace = FALSE)

  # Vector para asignar a qué grupo pertenece cada gen DEA (NA si no es DEG)
  grupo_deg_vector <- rep(NA_character_, n_genes)
  start <- 1
  for (i in seq_len(n_grupos)) {
    inds <- indices_deg[start:(start + n_deg_por_grupo[i] - 1)]
    grupo_deg_vector[inds] <- paste0("grupo", i)      # asigna grupo i a genes DEA correspondientes
    start <- start + n_deg_por_grupo[i]
  }

  # Para cada grupo, aplicar fold change multiplicativo a sus genes DEA en muestras del grupo
  # Y reducir la expresión de esos genes en los otros grupos a 0.3 veces el valor base (expresión baja)
  for (i in seq_len(n_grupos)) {
    grupo_name <- paste0("grupo", i)
    inds_genes <- which(grupo_deg_vector == grupo_name)    # genes DEA del grupo i
    inds_muestras <- which(grupos == grupo_name)           # muestras del grupo i

    # Expresión alta en su grupo: multiplicar fold change
    tpm_base[inds_genes, inds_muestras] <- tpm_base[inds_genes, inds_muestras] * cambios_esperados[i]

    # Expresión baja en los otros grupos: multiplicar por 0.3
    otros_grupos <- setdiff(seq_len(n_grupos), i)
    for (g2 in otros_grupos) {
      inds_muestras_g2 <- which(grupos == paste0("grupo", g2))
      tpm_base[inds_genes, inds_muestras_g2] <- tpm_base[inds_genes, inds_muestras_g2] * 0.3
    }
  }

  # Añadir ruido técnico aditivo normal a toda la matriz TPM simulada
  ruido <- matrix(rnorm(n_genes * n_muestras, mean = 0, sd = ruido_sd), nrow = n_genes, ncol = n_muestras)
  tpm_base <- tpm_base + ruido

  # Forzar a que no haya valores negativos (no tiene sentido expresión < 0)
  tpm_base[tpm_base < 0] <- 0

  # Convertir matriz TPM en data frame y asignar nombres a filas y columnas
  tpm_df <- as.data.frame(tpm_base)
  rownames(tpm_df) <- genes
  colnames(tpm_df) <- muestras

  # Data frame con metadata de las muestras y su grupo
  metadata <- data.frame(
    sample = muestras,
    group = grupos,
    stringsAsFactors = FALSE
  )

  # Data frame con info de genes: si es DEG y grupo asignado
  info_deg <- data.frame(
    gene_id = genes,
    grupo_deg = grupo_deg_vector,
    es_DEG = !is.na(grupo_deg_vector)
  )

  # Retornar lista con los tres data frames generados
  return(list(tpm_df = tpm_df, metadata = metadata, info_deg = info_deg))
}
