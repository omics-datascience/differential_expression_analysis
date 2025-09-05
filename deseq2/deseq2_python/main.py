# main.py

import os
# Importa las funciones desde tus nuevos módulos
from preprocessing import prepare_data
from analisis_dge import run_deseq_analysis

# --- PARÁMETROS DE CONFIGURACIÓN ---
# Rutas a los archivos de entrada y directorio de salida
GCT_FILE_GZ = "/home/mauri/Documentos/differential_expression_analysis/example/datasets/brain_cptac_gdc/data_mrna_seq_read_counts.txt"
METADATA_FILE = "/home/mauri/Documentos/differential_expression_analysis/example/datasets/brain_cptac_gdc/metadata.tsv"
OUTPUT_DIR = "/home/mauri/Documentos/differential_expression_analysis/example/results_python/brain_cptac_gdc"

# Parámetros para la comparación de grupos
COLUMNA_DE_GRUPO = "SEX"  # Nombre de la columna en 'METADATA_FILE' que define los grupos.
GRUPO1 = "Male"          # Nombre del primer grupo (usado como control/referencia).
GRUPO2 = "Female"           # Nombre del segundo grupo a comparar.


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

        # 2. Análisis de Expresión Diferencial
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
        return 

    print(f"\nAnálisis completado! Revisa los resultados en la carpeta: '{OUTPUT_DIR}'")


if __name__ == "__main__":
    main()