import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import numpy as np

# --- CONFIGURACIÓN ---
# Reemplaza estos nombres con los de tus archivos
archivo_r = '/home/mauri/Documentos/differential_expression_analysis/deseq2/example/results_r/resultados_completos_metastasis_vs_primary.csv'
archivo_python = '/home/mauri/Documentos/differential_expression_analysis/deseq2/example/results_python/resultados_metastasis_vs_primary.csv'
umbral_padj = 0.01
# ---------------------

# 1. Cargar los datos de los archivos CSV
try:
    df_r = pd.read_csv(archivo_r)
    df_py = pd.read_csv(archivo_python)
except FileNotFoundError as e:
    print(f"Error: No se pudo encontrar el archivo {e.filename}.")
    print("Por favor, asegúrate de que los nombres de archivo y la ubicación sean correctos.")
    exit()

print("--- Comparación de Genes Significativos ---")

# Filtrar genes significativos en cada dataframe
df_r_significativos = df_r[df_r['padj'] < umbral_padj]
df_py_significativos = df_py[df_py['padj'] < umbral_padj]

# Convertir las listas de nombres de genes a 'sets' para una comparación eficiente
genes_r_set = set(df_r_significativos['Name'])
genes_py_set = set(df_py_significativos['Name'])

# Realizar las comparaciones
genes_comunes = genes_r_set.intersection(genes_py_set)
genes_solo_r = genes_r_set.difference(genes_py_set)
genes_solo_py = genes_py_set.difference(genes_r_set)

# Calcular los porcentajes de solapamiento
porcentaje_en_r = (len(genes_comunes) / len(genes_py_set)) * 100 if len(genes_py_set) > 0 else 0
porcentaje_en_py = (len(genes_comunes) / len(genes_r_set)) * 100 if len(genes_r_set) > 0 else 0

# Imprimir los resultados en la consola
print(f"Umbral de p-valor ajustado: {umbral_padj}")
print(f"Total de genes significativos en R: {len(genes_r_set)}")
print(f"Total de genes significativos en Python: {len(genes_py_set)}")
print("-" * 35)
print(f"Genes significativos en AMBOS análisis: {len(genes_comunes)} ")
print(f"Genes significativos SÓLO en R: {len(genes_solo_r)}")
print(f"Genes significativos SÓLO en Python: {len(genes_solo_py)}")
print("-" * 35)
print("Métricas de Superposición:")
print(f"{porcentaje_en_py:.2f}% de los genes de R también fueron encontrados por Python.")
print(f"{porcentaje_en_r:.2f}% de los genes de Python también fueron encontrados por R.")
print("-" * 35)

# print(f"Nombres de genes SÓLO en R: {list(genes_solo_r)}")
# print(f"Nombres de genes SÓLO en Python: {list(genes_solo_py)}")

# 2. Crear el gráfico de dispersión comparativo
print("\n--- Creando gráfico de dispersión comparativo ---")
# Unir los dos dataframes usando la columna 'Name'
df_comparativo = pd.merge(
    df_r,
    df_py,
    on='Name',
    suffixes=('_r', '_py')
)
df_comparativo.dropna(subset=['log2FoldChange_r', 'log2FoldChange_py'], inplace=True)

# Calcular el coeficiente de correlación de Pearson
corr, _ = pearsonr(df_comparativo['log2FoldChange_r'], df_comparativo['log2FoldChange_py'])

# Crear el gráfico
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(8, 8))
sns.scatterplot(
    data=df_comparativo, x='log2FoldChange_r', y='log2FoldChange_py',
    alpha=0.5, s=15, ax=ax
)

# Añadir línea de referencia
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),
    np.max([ax.get_xlim(), ax.get_ylim()]),
]
ax.plot(lims, lims, 'r--', alpha=0.75, zorder=0, label='Correlación Perfecta (y=x)')

# Personalizar el gráfico
ax.set_title('Comparación de log2 Fold Change: R vs. Python', fontsize=16, weight='bold')
ax.set_xlabel('log2FC (DESeq2 en R)', fontsize=12)
ax.set_ylabel('log2FC (PyDESeq2 en Python)', fontsize=12)
ax.set_aspect('equal', 'box')
ax.legend()
ax.text(0.05, 0.95, f'Correlación de Pearson (todos los genes): {corr:.4f}',
        transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

# Guardar y mostrar el gráfico
plt.savefig('comparacion_log2fc.png', dpi=300)

# Imprimir resumen final en la consola
print(f"Gráfico de dispersión guardado como 'comparacion_log2fc.png'")
print(f"Se compararon los `log2FoldChange` de {len(df_comparativo)} genes presentes en ambos archivos.")
print(f"El coeficiente de correlación de Pearson es: {corr:.4f}")