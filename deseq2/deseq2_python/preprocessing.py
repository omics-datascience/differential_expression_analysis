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
    try:
        counts_df = pd.read_csv(gct_path, sep='\t', skiprows=2, index_col=0)
    except FileNotFoundError:
        raise FileNotFoundError(f"ERROR: No se encontró el archivo GCT en '{gct_path}'")

    # Limpiar y asegurar que la matriz es numérica
    count_matrix = counts_df.drop(columns=['Description']).apply(pd.to_numeric)

    # Sincronizar muestras entre la matriz de conteos y los metadatos
    common_samples = count_matrix.columns.intersection(metadata_filtered.index)

    count_matrix_filtered = count_matrix[common_samples]
    metadata_final = metadata_filtered.loc[common_samples]

    # Asegurar que el orden es idéntico (requisito de PyDESeq2)
    metadata_final = metadata_final.reindex(count_matrix_filtered.columns)

    # Convertir la columna de grupo a categórica y establecer el nivel de referencia (control)
    metadata_final[group_col] = pd.Categorical(metadata_final[group_col], categories=[group1, group2])

    return count_matrix_filtered, metadata_final