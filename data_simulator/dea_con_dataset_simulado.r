#!/usr/bin/env Rscript

# ------------------ Librerías -------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(pROC)
})

# ------------------ Cargar scripts --------------------------
source("data_simulator/TMP_multigrupos.r")
source("differential_expression_analysis/dea_with_limma.r")
source("process_results/process_limma_results.r")

# ------------------ Configuración ---------------------------
n_datasets       <- 10

umbral_significancia <- 0.05
# umbral_significancia: Este es el umbral de significancia que se usará para determinar si un gen es
# diferencialmente expresado.

# Parámetros de simulación
n_genes = 8000
n_muestras = 200
n_grupos = 4      # sólo dos grupos: grupo1 y grupo2 (ej. control y caso)
n_deg = 80        # Número total de genes que serán diferencialmente expresados
cambios_esperados = c(5, 5, 5, 5)  # fold change para cada grupo
ruido_sd = 0 # Mucha dispersión en los datos si es > 5
# n_deg: cuántos genes queremos que sean diferentes entre los dos grupos
# cambio_esperado: cuánto queremos que cambie la expresión de los genes en el grupo "caso". Si 
#   ponemos cambio_esperado = 4, eso significa que esos genes aparecerán 4 veces más 
#   expresados en los casos que en los controles.
# ruido_sd: Agrega variación aleatoria a todos los datos. Cuanto más alto sea ruido_sd, más
#   dispersos serán los resultados, aunque los promedios sean iguales. En este caso, no 
#   queremos ruido, así que lo dejamos en 0.

dir.create("data_simulator/datasets", showWarnings = FALSE)
dir.create("data_simulator/output", showWarnings = FALSE)

nombre_base      <- paste0("simu_",n_grupos,"_grupos")
# ------------------ Funciones -------------------------------
generar_dataset <- function(i, base_name) {
  nombre_dataset <- paste0(base_name, "_", i)
  dir.create(file.path("data_simulator/datasets", nombre_dataset), showWarnings = FALSE)

  set.seed(1000 + i)
  semilla_aleatoria <- sample(1:1000000, 1)

  resultado <- simular_tpm_multigrupos(
    n_genes = n_genes,
    n_muestras = n_muestras,
    n_grupos = n_grupos,
    n_deg = n_deg,
    cambios_esperados = cambios_esperados,
    ruido_sd = ruido_sd,
    semilla = semilla_aleatoria
  )

  list(
    nombre = nombre_dataset,
    semilla = semilla_aleatoria,
    tpm = resultado$tpm_df,
    metadata = resultado$metadata,
    info_deg = resultado$info_deg
  )
}

guardar_datos_simulados <- function(nombre, tpm, metadata, info_deg, semilla) {
  write.csv(tpm, file.path("data_simulator/datasets", nombre, "matriz_tpm_simulada.csv"), row.names = TRUE)
  write.csv(metadata, file.path("data_simulator/datasets", nombre, "metadata_muestras.csv"), row.names = FALSE)
  write.csv(info_deg, file.path("data_simulator/datasets", nombre, "info_genes_deg.csv"), row.names = FALSE)

  # Guardar semilla utilizada
  writeLines(as.character(semilla), con = file.path("data_simulator/datasets", nombre, "semilla_usada.txt"))
}

explorar_dataset_simulado <- function(tpm_df, nombre_dataset, es_log2 = FALSE) {
  nombre_archivo <- ifelse(es_log2, 
                           paste0(nombre_dataset, "_log2_tpm_distribution.png"), 
                           paste0(nombre_dataset, "_tpm_distribution.png"))
  titulo <- ifelse(es_log2,
                   paste("Distribución de log2(TPM+1) -", nombre_dataset),
                   paste("Distribución de TPM -", nombre_dataset))
  
  png(filename = file.path("data_simulator/datasets", nombre_dataset, nombre_archivo),
      width = 800, height = 600)
  
  hist(as.numeric(as.matrix(tpm_df)), breaks = 100,
       main = titulo,
       xlab = ifelse(es_log2, "log2(TPM + 1)", "TPM"),
       col = ifelse(es_log2, "tomato", "steelblue"))
  
  dev.off()
}

evaluar_resultados <- function(all_results, info_deg) {
  all_results$gene_id <- rownames(all_results)

  eval_tbl <- all_results %>%
    left_join(info_deg[, c("gene_id", "es_DEG")], by = "gene_id") %>%
    mutate(
      es_DEG = ifelse(is.na(es_DEG), FALSE, es_DEG),
      sig = adj.P.Val < umbral_significancia
    )

  tp <- sum(eval_tbl$es_DEG & eval_tbl$sig)
  fn <- sum(eval_tbl$es_DEG & !eval_tbl$sig)
  fp <- sum(!eval_tbl$es_DEG & eval_tbl$sig)
  tn <- sum(!eval_tbl$es_DEG & !eval_tbl$sig)

  sens <- tp / (tp + fn)
  prec <- tp / (tp + fp)
  f1 <- 2 * prec * sens / (prec + sens)

  roc_obj <- roc(eval_tbl$es_DEG, -log10(eval_tbl$adj.P.Val), quiet = TRUE)

  list(
    eval_tbl = eval_tbl,
    metrics = list(tp = tp, fn = fn, fp = fp, tn = tn, sens = sens, prec = prec, f1 = f1),
    roc = roc_obj
  )
}

guardar_resultados <- function(nombre, top_50, eval_tbl, roc_obj) {
  write.csv(top_50, file = file.path("data_simulator/output", paste0(nombre, "_group_top50_results.csv")), row.names = TRUE)
  write.table(eval_tbl, file = file.path("data_simulator/output", paste0(nombre, "_eval_limma.csv")), sep = "\t", row.names = FALSE)

  png(filename = file.path("data_simulator/output", paste0(nombre, "_roc_curve.png")), width = 800, height = 600)
  plot(roc_obj, main = paste("ROC -", nombre))
  dev.off()
}

# ------------------ Ejecución principal ----------------------
for (i in seq_len(n_datasets)) {

  ds <- generar_dataset(i, nombre_base)
  guardar_datos_simulados(ds$nombre, ds$tpm, ds$metadata, ds$info_deg, ds$semilla)

  explorar_dataset_simulado(ds$tpm, ds$nombre, FALSE)
  ds$tpm_log = log2(ds$tpm + 1)
  explorar_dataset_simulado(ds$tpm_log, ds$nombre, TRUE)

  all_results <- perform_differential_expression(ds$tpm_log, ds$metadata, "group")
  # all_results: lista de data frames con resultados de DEA por contraste
  # Cada data frame contiene columnas como logFC, AveExpr, t, P.Value, adj.P.Val, etc.
  # Cada data frame corresponde a un contraste entre grupos (ej. grupo2 vs grupo1)

  for (i in seq_along(all_results)) {
    nombre <- names(all_results)[i]
    df <- all_results[[i]]
    
    ordered_results <- order_results_by_pvalue(df)
    top_50 <- get_top_50_by_pvalue(ordered_results)

    eval_result <- evaluar_resultados(df, ds$info_deg)
    eval_tbl <- eval_result$eval_tbl
    met <- eval_result$metrics
    roc_obj <- eval_result$roc
    auc <- auc(roc_obj)

    cat("Comparación:", nombre, "\n")
    message(sprintf(
      "Dataset %s — TP=%d  FN=%d  FP=%d  TN=%d  Sens=%.2f  Prec=%.2f  F1=%.2f  AUC=%.2f",
      ds$nombre, met$tp, met$fn, met$fp, met$tn, met$sens, met$prec, met$f1, auc
    ))

    guardar_resultados(ds$nombre, top_50, eval_tbl, roc_obj)

    rm(df, ordered_results, top_50, eval_tbl, eval_result, roc_obj)
    gc()

  }
}