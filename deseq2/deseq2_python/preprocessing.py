# preprocessing.py
import pandas as pd


def prepare_data(gct_path, metadata_path, group_col, group1, group2):
    """
    Carga, valida, filtra y sincroniza los datos de conteos y metadatos.

    Returns:
        tuple: Una tupla conteniendo el DataFrame de conteos filtrado y el 
               DataFrame de metadatos final y sincronizado.
    """
    # Carga y filtrado de metadatos
    try:
        all_metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
    except FileNotFoundError:
        raise FileNotFoundError(f"ERROR: No se encontró el archivo de metadatos en '{metadata_path}'")

    # Validar que la columna y los grupos especificados existan
    if group_col not in all_metadata.columns:
        raise ValueError(f"ERROR: La columna '{group_col}' no existe en el archivo de metadatos.")

    unique_groups = all_metadata[group_col].unique()
    if group1 not in unique_groups:
        raise ValueError(f"ERROR: El grupo '{group1}' no se encuentra en la columna '{group_col}'.")
    if group2 not in unique_groups:
        raise ValueError(f"ERROR: El grupo '{group2}' no se encuentra en la columna '{group_col}'.")

    # Filtrar metadatos para conservar únicamente las muestras de los grupos de interés
    metadata_filtered = all_metadata[all_metadata[group_col].isin([group1, group2])].copy()

    print(f"INFO - Muestras encontradas: {(metadata_filtered[group_col] == group1).sum()} para '{group1}' y "
          f"{(metadata_filtered[group_col] == group2).sum()} para '{group2}'.")

    # Carga y filtrado de la matriz de conteos
    matrix_type = None
    try:
        if gct_path.endswith("gct") or gct_path.endswith("gct.gz"):
            matrix_type = "GCT"
            counts_df = pd.read_csv(gct_path, sep='\t', skiprows=2, index_col=0)
        elif gct_path.endswith("csv") or gct_path.endswith("csv.gz"):
            matrix_type = "CSV"
            counts_df = pd.read_csv(gct_path, index_col=0)
        elif gct_path.endswith("tsv") or gct_path.endswith("tsv.gz"):
            matrix_type = "TSV"
            counts_df = pd.read_csv(gct_path, sep='\t', index_col=0)
        elif gct_path.endswith("txt") or gct_path.endswith("txt.gz"):
            matrix_type = "TSV"
            counts_df = pd.read_csv(gct_path, sep='\t', index_col=0)
        else:
            raise ValueError("ERROR: El archivo de conteos debe estar en formato .gct, .csv, .tsv o .txt (o sus versiones comprimidas).")
        # index_col=0 indica que la primera columna del archivo (índice 0) debe usarse como 
        # el índice de la tabla (DataFrame). En este contexto, esta columna contiene los 
        # nombres de los genes
    except FileNotFoundError:
        raise FileNotFoundError(f"ERROR: No se encontró el archivo GCT en '{gct_path}'")

    # Limpiar y asegurar que la matriz es numérica
    # Elimina la columna llamada "Description" del DataFrame. Esta columna contiene 
    # descripciones de los genes que no son necesarias para el análisis numérico
    if matrix_type == "GCT" and 'Description' in counts_df.columns:
        count_matrix = counts_df.drop(columns=['Description']).apply(pd.to_numeric)
    else:
        count_matrix = counts_df.apply(pd.to_numeric)

    # Sincronizar muestras entre la matriz de conteos y los metadatos
    # Encuentra los nombres de las muestras que están presentes en ambos DataFrames (en las
    # columnas de count_matrix y en el índice de metadata_filtered). Esto es crucial para 
    # asegurar que solo se analicen las muestras que tienen tanto datos de conteo como 
    # información descriptiva (grupos).
    common_samples = count_matrix.columns.intersection(metadata_filtered.index)

    count_matrix_filtered = count_matrix[common_samples]
    metadata_final = metadata_filtered.loc[common_samples]

    # Asegurar que el orden es idéntico (requisito de PyDESeq2)
    # Reordena las filas del DataFrame de metadatos (metadata_final) para que sigan 
    # exactamente el mismo orden que las columnas del DataFrame de conteos (count_matrix_filtered). 
    # Algunas herramientas de análisis, como PyDESeq2, requieren que este orden sea idéntico
    # para funcionar correctamente.
    metadata_final = metadata_final.reindex(count_matrix_filtered.columns)

    # Convertir la columna de grupo a categórica y establecer el nivel de referencia (control)
    metadata_final[group_col] = pd.Categorical(metadata_final[group_col], categories=[group1, group2])

    return count_matrix_filtered, metadata_final