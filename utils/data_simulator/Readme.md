# Simulación de Experimentos Case-Control con RNA-seq

El script **dea_con_dataset_simulado.r** permite generar datasets simulados de expresión génica tipo **case-control**, realizar análisis de expresión diferencial con **limma**, y evaluar el rendimiento de detección usando métricas como **F1-score** y **AUC**. Es útil para validar el proceso de analisis de expresion diferencial desarrollado.

---

## ¿Qué hace el script?

Para cada dataset simulado, el script:

1. Genera una matriz TPM de expresión génica con diferencias entre grupos.
2. Guarda los datos simulados, metadatos y genes diferenciales esperados.
3. Realiza análisis de expresión diferencial con `limma`.
4. Evalúa sensibilidad, precisión, F1-score y AUC.
5. Guarda resultados tabulares y una curva ROC por dataset.

---

## Parámetros configurables

Los siguientes parámetros se configuran directamente en la sección de **configuración del script**:

| Parámetro              | Descripción                                                                 | Ejemplo |
|------------------------|-----------------------------------------------------------------------------|---------|
| `n_datasets`           | Número de datasets simulados que se generarán.                             | `5`     |
| `nombre_base`          | Nombre base para las carpetas y archivos de salida.                        | `"0_dataset_simulado_case_control"` |
| `n_genes`              | Número total de genes en cada simulación.                                  | `8000`  |
| `n_muestras`           | Número total de muestras por dataset.                                      | `150`   |
| `n_deg`                | Número de genes que serán diferencialmente expresados (DEGs).              | `5`     |
| `cambio_esperado`      | Factor de cambio de expresión en los casos respecto a controles.           | `5`     |
| `ruido_sd`             | Desviación estándar del ruido aleatorio agregado a la expresión.           | `0`     |
| `umbral_significancia`| Umbral para considerar significativa una diferencia de expresión (`adj.P.Val`). | `0.05`  |

---

## Generacion de valores de expresion

Resumen de la generación de valores en el dataset simulado

1. Fijación de la semilla (set.seed)
   - Permite que los valores generados sean reproducibles.
   - semilla = 42 por defecto para controlar toda la aleatoriedad.
2. Generación de expresión de base (TPM)
   - Se simulan los valores de expresión de todos los genes en todas las muestras usando una distribución exponencial:
        rexp(n_genes * n_muestras, rate = 1/20)
   - Esto genera muchos valores bajos (simulando genes con baja expresión) y unos pocos altos (genes mas expresados). Simula una distribución natural de expresión, algunos genes tienen muy baja expresión (TPM < 1), otros moderada (~10-100) y no hay ninguno con muy alta expresion (>1000).
   - Estos valores van a parar a una matriz de TPM base (tpm_base).
3. Asignación aleatoria de genes diferencialmente expresados (DEGs)
   - Se eligen al azar n_deg genes que serán los diferencialmente expresados:
       sample(n_genes, n_deg)
   - Estos genes van a tener mayor expresión en las muestras del grupo caso.
4. Aplicación de fold change:
   - A los genes DEGs en el grupo caso se les multiplica su valor de TPM por un cambio_esperado (por ejemplo, 4 veces más expresión que en control):
5. Agregado de ruido técnico:
   - Se genera un ruido aleatorio con distribución normal (media 0, desviación estándar ruido_sd) para simular variación experimental.
   - Este ruido se suma a todos los valores de TPM.
6. Corrección de valores negativos:
   - Como el ruido puede generar valores negativos, se reemplazan por cero.

## Estructura de carpetas y archivos

Después de la ejecución, se generan:

- Carpeta `datasets/` con subcarpetas por dataset, que contienen:
  - `matriz_tpm_simulada.csv`
  - `metadata_muestras.csv`
  - `info_genes_deg.csv`
  - `semilla_usada.txt`

- Carpeta `output/` con:
  - Resultados del análisis (`*_top50_results.csv`)
  - Evaluación del desempeño (`*_eval_limma.csv`)
  - Curva ROC (`*_roc_curve.png`)

---

## Ejecución desde línea de comandos

```bash
Rscript dea_con_dataset_simulado.R
```

---

## Dependencias

- `dplyr`
- `pROC`
- Scripts auxiliares:
  - `data_simulator/TMP_case_control.r`
  - `differential_expression_analysis/dea_with_limma.r`
  - `process_results/process_limma_results.r`

## Funcionamiento de creacion de datasets simulados  

El script generador de datasets **TPM** sintéticos para estudios *case‑control* de expresión génica es **TMP_case_control.r**

---

### ¿Qué hace?  

Crea una matriz de expresión génica (TPM) con dos grupos de muestras —`control` y `caso`— donde un conjunto definido de genes está diferencialmente expresado (DEG) en el grupo `caso`.  
Se añade ruido técnico para reflejar variabilidad experimental y se devuelven:

| Objeto | Descripción |
|--------|-------------|
| `tpm_df` | Data frame de dimensión **genes × muestras** con los valores TPM simulados. |
| `metadata` | Data frame con la asignación de cada muestra a su grupo (`control` / `caso`). |
| `info_deg` | Data frame que marca qué genes son DEG (`TRUE/FALSE`). |

---

### Argumentos

| Argumento | Tipo / Valor por defecto | Descripción |
|-----------|--------------------------|-------------|
| `n_genes` | `numeric` · `10000` | Número total de genes a simular. |
| `n_muestras` | `numeric` · `12` | Número total de muestras (se dividen entre grupos). |
| `n_deg` | `numeric` · `500` | Genes que serán diferencialmente expresados. |
| `cambio_esperado` | `numeric` · `4` | **Fold change** aplicado a los DEGs en el grupo `caso`. |
| `ruido_sd` | `numeric` · `5` | Desviación estándar del ruido normal aditivo. |
| `semilla` | `numeric` · `42` | Semilla para garantizar reproducibilidad. |

---

### Flujo de simulación

1. **Distribución base** de TPM: valores `rexp(rate = 1/20)`.  
2. Selección aleatoria de `n_deg` genes como DEG.  
3. Aplicación de `cambio_esperado` a los DEGs en las muestras `caso`.  
4. Adición de **ruido** normal con `sd = ruido_sd`.  
5. Reemplazo de valores negativos por 0.

---

### Uso  

```r
resultado <- simular_tpm_case_control(
  n_genes = 8000,
  n_muestras = 20,
  n_deg = 300,
  cambio_esperado = 3,
  ruido_sd = 2,
  semilla = 123
)

# Objetos resultantes
tpm      <- resultado$tpm_df      # matriz genes × muestras
metadata <- resultado$metadata    # tabla de grupos
info_deg <- resultado$info_deg    # genes diferenciales
```

---

### Reproducibilidad

La semilla (`semilla`) se fija internamente con `set.seed()`, garantizando que se obtenga exactamente el mismo dataset cuando se repiten los mismos parámetros.  
