# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# from scipy.stats import pearsonr
# import numpy as np

# archivo_r = 'results_r/resultados_completos_metastasis_vs_primary.csv'
# archivo_python = 'results_python/resultados_metastasis_vs_primary.csv'
# umbral_padj = 0.01

# # Cargar los datos de los archivos CSV
# try:
#     df_r = pd.read_csv(archivo_r)
#     df_py = pd.read_csv(archivo_python)
# except FileNotFoundError as e:
#     print(f"Error: No se pudo encontrar el archivo {e.filename}.")
#     print("Por favor, asegúrate de que los nombres de archivo y la ubicación sean correctos.")
#     exit()

# print("--- Comparación de Genes Significativos ---")

# # Filtrar genes significativos en cada dataframe
# df_r_significativos = df_r[df_r['padj'] < umbral_padj]
# df_py_significativos = df_py[df_py['padj'] < umbral_padj]

# # Convertir las listas de nombres de genes a 'sets' para una comparación eficiente
# genes_r_set = set(df_r_significativos['Name'])
# genes_py_set = set(df_py_significativos['Name'])

# # Realizar las comparaciones
# genes_comunes = genes_r_set.intersection(genes_py_set)
# genes_solo_r = genes_r_set.difference(genes_py_set)
# genes_solo_py = genes_py_set.difference(genes_r_set)

# # Calcular los porcentajes de solapamiento
# porcentaje_en_r = (len(genes_comunes) / len(genes_py_set)) * 100 if len(genes_py_set) > 0 else 0
# porcentaje_en_py = (len(genes_comunes) / len(genes_r_set)) * 100 if len(genes_r_set) > 0 else 0

# # Imprimir los resultados en la consola
# print(f"Umbral de p-valor ajustado: {umbral_padj}")
# print(f"Total de genes significativos en R: {len(genes_r_set)}")
# print(f"Total de genes significativos en Python: {len(genes_py_set)}")
# print("-" * 38)
# print(f"Genes significativos en AMBOS análisis: {len(genes_comunes)} ")
# print(f"Genes significativos SÓLO en R: {len(genes_solo_r)}")
# print(f"Genes significativos SÓLO en Python: {len(genes_solo_py)}")
# print("-" * 38)
# print("Métricas de Superposición:")
# print(f"{porcentaje_en_py:.2f}% de los genes de R también fueron encontrados por Python.")
# print(f"{porcentaje_en_r:.2f}% de los genes de Python también fueron encontrados por R.")
# print("-" * 38)

# print(f"Nombres de genes SÓLO en R: {list(genes_solo_r)}")
# print(f"Nombres de genes SÓLO en Python: {list(genes_solo_py)}")

# # Crear el gráfico de dispersión comparativo
# print("\n--- Creando gráfico de dispersión comparativo ---")
# # Unir los dos dataframes usando la columna 'Name'
# df_comparativo = pd.merge(
#     df_r,
#     df_py,
#     on='Name',
#     suffixes=('_r', '_py')
# )
# df_comparativo.dropna(subset=['log2FoldChange_r', 'log2FoldChange_py'], inplace=True)

# # Calcular el coeficiente de correlación de Pearson
# corr, _ = pearsonr(df_comparativo['log2FoldChange_r'], df_comparativo['log2FoldChange_py'])

# # Crear el gráfico
# plt.style.use('seaborn-v0_8-whitegrid')
# fig, ax = plt.subplots(figsize=(8, 8))
# sns.scatterplot(
#     data=df_comparativo, x='log2FoldChange_r', y='log2FoldChange_py',
#     alpha=0.5, s=15, ax=ax
# )

# # Añadir línea de referencia
# lims = [
#     np.min([ax.get_xlim(), ax.get_ylim()]),
#     np.max([ax.get_xlim(), ax.get_ylim()]),
# ]
# ax.plot(lims, lims, 'r--', alpha=0.75, zorder=0, label='Correlación Perfecta (y=x)')

# # Personalizar el gráfico
# ax.set_title('Comparación de log2 Fold Change: R vs. Python', fontsize=16, weight='bold')
# ax.set_xlabel('log2FC (DESeq2 en R)', fontsize=12)
# ax.set_ylabel('log2FC (PyDESeq2 en Python)', fontsize=12)
# ax.set_aspect('equal', 'box')
# ax.legend()
# ax.text(0.05, 0.95, f'Correlación de Pearson (todos los genes): {corr:.4f}',
#         transform=ax.transAxes, fontsize=12,
#         verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

# # Guardar y mostrar el gráfico
# plt.savefig('comparacion_log2fc.png', dpi=300)

# # Imprimir resumen final en la consola
# print("Gráfico de dispersión guardado como 'comparacion_log2fc.png'")
# print(f"Se compararon los `log2FoldChange` de {len(df_comparativo)} genes presentes en ambos archivos.")
# print(f"El coeficiente de correlación de Pearson es: {corr:.4f}")

# print("""La correlación de Pearson entre los `log2FoldChange` de los dos métodos nos indica qué tan alineados están los resultados en cuanto a magnitud y dirección de los cambios de expresión génica. 
# Usamos **`log2FoldChange`** porque representa directamente el cambio en la expresión de los genes, y es comparable entre diferentes métodos:
# Un valor cercano a **1** indica alta concordancia, mientras que un valor cercano a **0** sugiere poca relación.
# """)


import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import numpy as np

# Directorios base
dir_r = "results_r"
dir_py = "results_python"
umbral_padj = 0.01


# Función para cargar CSV y asegurarse de que la primera columna sin nombre se llame "Name"
def cargar_csv(archivo):
    df = pd.read_csv(archivo)

    # Detectar si la primera columna no tiene nombre o se llama Unnamed
    primera_col = df.columns[0]
    if primera_col == "" or primera_col.startswith("Unnamed"):
        df = df.rename(columns={primera_col: "Name"})
    elif "Name" not in df.columns:
        # Si no hay ninguna columna llamada Name, forzar a que la primera sea Name
        df = df.rename(columns={primera_col: "Name"})
    
    return df


# Función para encontrar el único CSV dentro de cada subcarpeta
def buscar_archivos(directorio):
    archivos = []
    for carpeta in os.listdir(directorio):
        ruta_carpeta = os.path.join(directorio, carpeta)
        if os.path.isdir(ruta_carpeta):
            encontrados = glob.glob(os.path.join(ruta_carpeta, "*.csv"))
            if len(encontrados) == 1:
                archivos.append((carpeta, encontrados[0]))
            elif len(encontrados) > 1:
                print(f"Atención: en {ruta_carpeta} hay más de un CSV. Se usará el primero.")
                archivos.append((carpeta, encontrados[0]))
    return dict(archivos)

# Buscar archivos en ambos directorios
archivos_r = buscar_archivos(dir_r)
archivos_py = buscar_archivos(dir_py)

# Comparar solo pruebas que existan en ambos directorios
pruebas_comunes = set(archivos_r.keys()).intersection(archivos_py.keys())

if not pruebas_comunes:
    print("No se encontraron pruebas comunes entre results_r y results_python.")
    exit()

for prueba in pruebas_comunes:
    archivo_r = archivos_r[prueba]
    archivo_python = archivos_py[prueba]

    print(f"\n=== Comparando prueba: {prueba} ===")
    
    # Cargar los datos
    try:
        df_r = cargar_csv(archivo_r)
        df_py = cargar_csv(archivo_python)
    except FileNotFoundError as e:
        print(f"Error: No se pudo encontrar el archivo {e.filename}.")
        continue

    # --- Comparación de genes significativos ---
    df_r_significativos = df_r[df_r['padj'] < umbral_padj]
    df_py_significativos = df_py[df_py['padj'] < umbral_padj]

    genes_r_set = set(df_r_significativos['Name'])
    genes_py_set = set(df_py_significativos['Name'])

    genes_comunes = genes_r_set.intersection(genes_py_set)
    genes_solo_r = genes_r_set.difference(genes_py_set)
    genes_solo_py = genes_py_set.difference(genes_r_set)

    porcentaje_en_r = (len(genes_comunes) / len(genes_py_set)) * 100 if len(genes_py_set) > 0 else 0
    porcentaje_en_py = (len(genes_comunes) / len(genes_r_set)) * 100 if len(genes_r_set) > 0 else 0

    print(f"Umbral de p-valor ajustado: {umbral_padj}")
    print(f"Total de genes significativos en R: {len(genes_r_set)}")
    print(f"Total de genes significativos en Python: {len(genes_py_set)}")
    print("-" * 38)
    print(f"Genes significativos en AMBOS análisis: {len(genes_comunes)} ")
    print(f"Genes significativos SÓLO en R: {len(genes_solo_r)}")
    print(f"Genes significativos SÓLO en Python: {len(genes_solo_py)}")
    print("-" * 38)
    print("Métricas de Superposición:")
    print(f"{porcentaje_en_py:.2f}% de los genes de R también fueron encontrados por Python.")
    print(f"{porcentaje_en_r:.2f}% de los genes de Python también fueron encontrados por R.")
    print("-" * 38)

    # --- Gráfico comparativo ---
    df_comparativo = pd.merge(df_r, df_py, on="Name", suffixes=("_r", "_py"))
    df_comparativo.dropna(subset=["log2FoldChange_r", "log2FoldChange_py"], inplace=True)

    if df_comparativo.empty:
        print("No hay genes comunes para graficar.")
        continue

    corr, _ = pearsonr(df_comparativo["log2FoldChange_r"], df_comparativo["log2FoldChange_py"])
    print(f"Coeficiente de correlación de Pearson (log2FC R vs Python): {corr:.4f}")
    
    
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(8, 8))
    sns.scatterplot(
        data=df_comparativo, x="log2FoldChange_r", y="log2FoldChange_py",
        alpha=0.5, s=15, ax=ax
    )

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),
        np.max([ax.get_xlim(), ax.get_ylim()]),
    ]
    ax.plot(lims, lims, 'r--', alpha=0.75, zorder=0, label="Correlación Perfecta (y=x)")

    ax.set_title(f"Comparación de log2 Fold Change: R vs. Python\n({prueba})", fontsize=16, weight="bold")
    ax.set_xlabel("log2FC (DESeq2 en R)", fontsize=12)
    ax.set_ylabel("log2FC (PyDESeq2 en Python)", fontsize=12)
    ax.set_aspect("equal", "box")
    ax.legend()
    ax.text(0.05, 0.95, f"Correlación de Pearson: {corr:.4f}",
            transform=ax.transAxes, fontsize=12,
            verticalalignment="top", bbox=dict(boxstyle="round,pad=0.5", fc="wheat", alpha=0.5))

    # Guardar gráfico
    out_file = f"comparacion_log2fc_{prueba}.png"
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"Gráfico guardado como '{out_file}'")

