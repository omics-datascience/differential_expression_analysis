import pandas as pd
import argparse
import sys


def create_metadata_file(patient_file, sample_file, output_file="metadata.tsv"):
    """
    Une los datos de dos archivos TSV (paciente y muestra) para crear un archivo de metadatos unificado.
    Args:
        patient_file (str): La ruta al archivo data_clinical_patient.txt.
        sample_file (str): La ruta al archivo data_clinical_sample.txt.
        output_file (str): El nombre del archivo de salida.
    """
    try:
        df_patient = pd.read_csv(patient_file, sep='\t', skiprows=4)
        df_sample = pd.read_csv(sample_file, sep='\t', skiprows=4)

        if 'PATIENT_ID' not in df_patient.columns:
            print(f"Error: La columna 'PATIENT_ID' no se encuentra en el archivo {patient_file}", file=sys.stderr)
            return
        if 'PATIENT_ID' not in df_sample.columns or 'SAMPLE_ID' not in df_sample.columns:
            print(f"Error: Las columnas 'PATIENT_ID' y 'SAMPLE_ID' no se encuentran en el archivo {sample_file}", file=sys.stderr)
            return

        # Unir los dataframes basados en la columna 'PATIENT_ID'
        # Esto asocia los datos del paciente con sus muestras correspondientes
        df_merged = pd.merge(df_sample[['PATIENT_ID', 'SAMPLE_ID']],
                             df_patient,
                             on='PATIENT_ID',
                             how='left')
        # Elimino duplicados
        df_merged_deduplicated = df_merged.drop_duplicates()

        # Reordenar las columnas para que 'SAMPLE_ID' sea la primera
        cols = ['SAMPLE_ID'] + [col for col in df_merged_deduplicated.columns if col != 'SAMPLE_ID']
        df_final = df_merged_deduplicated[cols]

        # Guardar el resultado en el archivo de salida
        df_final.to_csv(output_file, sep='\t', index=False)
        # Reordenar las columnas para que 'SAMPLE_ID' sea la primera
        cols = ['SAMPLE_ID'] + [col for col in df_merged.columns if col != 'SAMPLE_ID']
        df_final = df_merged[cols]

        # Guardar el resultado en el archivo de salida
        df_final.to_csv(output_file, sep='\t', index=False)
        print(f"Archivo '{output_file}' creado exitosamente.")
    except Exception as e:
        print(f"Ocurrió un error inesperado: {e}", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Crea un archivo de metadatos TSV uniendo datos de paciente y muestra.')
    parser.add_argument('patient_file', help='Ruta al archivo data_clinical_patient.txt')
    parser.add_argument('sample_file', help='Ruta al archivo data_clinical_sample.txt')
    parser.add_argument('-o', '--output', default='metadata.tsv', help='Nombre del archivo de salida (por defecto: metadata.tsv)')

    args = parser.parse_args()

    create_metadata_file(args.patient_file, args.sample_file, args.output)
