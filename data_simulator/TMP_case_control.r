simular_tpm_case_control<- function(
  n_genes = 8000,
  n_muestras = 60,
  n_deg = 50,
  cambio_esperado = 2.5,
  ruido_sd = 3,
  semilla = 42
) {
  #' Simula un dataset sintético de expresión génica en TPM para análisis diferencial.
  #'
  #' Simula una matriz de expresión TPM con dos grupos de muestras ('control' y 'caso'),
  #' donde un número específico de genes (n_deg) están diferencialmente expresados en el grupo 'caso'.
  #' Se añade ruido técnico para representar variabilidad realista.
  #'
  #' @param n_genes Número total de genes a simular.
  #' @param n_muestras Número total de muestras (se divide en dos grupos iguales o casi iguales).
  #' @param n_deg Número exacto de genes diferencialmente expresados (DEGs).
  #' @param cambio_esperado Fold change aplicado a los genes DEGs en el grupo 'caso'. El fold change (cambio 
  #' en la expresión) aplicado a los genes DE en el grupo caso respecto al grupo control. Por 
  #' ejemplo, 4 implica que los DEGs tendrán 4 veces más expresión en caso.
  #' @param ruido_sd Desviación estándar del ruido técnico aditivo aplicado a todos los valores TPM.
  #' @param semilla Semilla para la generación de números aleatorios (garantiza reproducibilidad).
  #'
  #' @return Una lista con tres elementos:
  #' \item{tpm_df}{Matriz de expresión TPM simulada (genes x muestras).}
  #' \item{metadata}{Data frame con el grupo de cada muestra ('control' o 'caso').}
  #' \item{info_deg}{Data frame indicando qué genes son diferencialmente expresados (TRUE/FALSE).}
  #'
  #' @examples
  #' resultado <- simular_tpm(n_genes=15000, n_muestras=10, n_deg=400, cambio_esperado=5, ruido_sd=2)
  #' tpm <- resultado$tpm_df
  #' metadata <- resultado$metadata
  #' info_deg <- resultado$info_deg

  if (n_deg > n_genes) {
    stop("El número de genes DE (n_deg) no puede ser mayor que el total de genes (n_genes).")
  }

  set.seed(semilla)

  # Crear nombres
  genes <- paste0("gen_", seq_len(n_genes))
  muestras <- paste0("muestra_", seq_len(n_muestras))

  # Definir grupos: mitad control, mitad caso (o casi)
  grupos <- c(rep("control", floor(n_muestras / 2)), rep("caso", ceiling(n_muestras / 2)))

  # Simulación base: distribución exponencial para TPM
  # Exponencial con media 30
  rate = 1/30
  tpm_base <- matrix(rexp(n_genes * n_muestras, rate = rate),
                    nrow = n_genes, ncol = n_muestras)

  # Variabilidad entre individuos (biología real)
  efecto_biologico <- rnorm(n_muestras, 0, 3)
  tpm_base <- sweep(tpm_base, 2, efecto_biologico, "+")

  # Seleccionar genes diferencialmente expresados
  indices_deg <- sample(n_genes, n_deg, replace = FALSE)

  # Índices de muestras caso
  indices_caso <- which(grupos == "caso")

  # Aplicar fold change a genes DE en grupo caso
  tpm_base[indices_deg, indices_caso] <- tpm_base[indices_deg, indices_caso] * cambio_esperado

  # Añadir ruido normal aditivo (Variabilidad de medición (secuenciación, RSEM))
  ruido <- matrix(rnorm(n_genes * n_muestras, mean = 0, sd = ruido_sd), nrow = n_genes, ncol = n_muestras)
  tpm_base <- tpm_base + ruido

  # No permitir valores negativos
  tpm_base[tpm_base < 0] <- 0

  # Construir data frames
  tpm_df <- as.data.frame(tpm_base)
  rownames(tpm_df) <- genes
  colnames(tpm_df) <- muestras

  metadata <- data.frame(
    sample = muestras,
    group = grupos,
    stringsAsFactors = FALSE
  )

  info_deg <- data.frame(
    gene_id = genes,
    es_DEG = seq_len(n_genes) %in% indices_deg
  )

  return(list(tpm_df = tpm_df, metadata = metadata, info_deg = info_deg))
}
