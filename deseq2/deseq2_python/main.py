# main.py

import os
# Importa las funciones desde tus nuevos módulos
from preprocessing import prepare_data
from analisis_dge import run_deseq_analysis

# --- PARÁMETROS DE CONFIGURACIÓN ---
# Rutas a los archivos de entrada y directorio de salida
GCT_FILE_GZ = "../example/CCLE_RNAseq_genes_counts_20180929.gct.gz"
METADATA_FILE = "../example/Cell_lines_annotations_20181226.txt"
OUTPUT_DIR = "../example/results_python"

# Parámetros para la comparación de grupos
COLUMNA_DE_GRUPO = "Pathology"  # Nombre de la columna en 'METADATA_FILE' que define los grupos.
GRUPO1 = "primary"              # Nombre del primer grupo (usado como control/referencia).
GRUPO2 = "metastasis"           # Nombre del segundo grupo a comparar.


def main():
    """Función principal para ejecutar el pipeline de análisis."""

    # Crear el directorio de salida si no existe
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"INFO - Directorio de resultados: '{OUTPUT_DIR}'")

    try:
        # 1. Preprocesamiento de datos (usando la función importada)
        print("\n--- Paso 1: Iniciando preprocesamiento de datos ---")
        count_matrix, metadata = prepare_data(
            gct_path=GCT_FILE_GZ,
            metadata_path=METADATA_FILE,
            group_col=COLUMNA_DE_GRUPO,
            group1=GRUPO1,
            group2=GRUPO2
        )

        # 2. Análisis de Expresión Diferencial (usando la función importada)
        print("\n--- Paso 2: Ejecutando análisis con PyDESeq2 ---")
        dds, results_df = run_deseq_analysis(
            counts=count_matrix,
            metadata=metadata,
            group_col=COLUMNA_DE_GRUPO,
            group1=GRUPO1,
            group2=GRUPO2
        )

        # 3. Guardado de resultados
        print("\n--- Paso 3: Guardando tabla de resultados ---")
        results_path = os.path.join(OUTPUT_DIR, f"resultados_{GRUPO2}_vs_{GRUPO1}.csv")
        results_df.to_csv(results_path)
        print(f"INFO - Tabla de resultados guardada en '{results_path}'")

    except (FileNotFoundError, ValueError) as e:
        print(f"\nCRITICAL - {e}")
        return  # Termina la ejecución si hay un error

    print(f"\nAnálisis completado! Revisa los resultados en la carpeta: '{OUTPUT_DIR}'")


if __name__ == "__main__":
    main()