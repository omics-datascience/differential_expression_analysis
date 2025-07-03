#!/bin/bash

# Lista de datasets a procesar, con valores repetidos porque cada uno tiene distintos atributos
datasets=("brca_tcga" "brca_tcga" "brca_tcga" "acc_tcga" "acc_tcga" "acc_tcga")

# Lista de atributos a analizar para cada dataset, se corresponde por indice con "datasets"
atributos=("MENOPAUSE_STATUS" "ETHNICITY" "SEX" "CLINICAL_STATUS_WITHIN_3_MTHS_SURGERY" "TUMOR_STATUS" "LATERALITY")

# NOTAS
# brca_tcga --> MENOPAUSE_STATUS --> 5 grupos (tamaño dataset: 195MB)
# brca_tcga --> ETHNICITY --> 3 grupos (tamaño dataset: 195MB)
# brca_tcga --> SEX --> 2 grupos (tamaño dataset: 195MB)
# acc_tcga --> CLINICAL_STATUS_WITHIN_3_MTHS_SURGERY --> 5 grupos (tamaño dataset: 14MB)
# acc_tcga --> TUMOR_STATUS --> 3 grupos (tamaño dataset: 14MB)
# acc_tcga --> LATERALITY --> 2 grupos (tamaño dataset: 14MB)

script="aed_limma_cbioportal_TMP.r"

output_dir="benchmark"

logs_dir="${output_dir}/logs"

# Archivo resultados
output_table="${output_dir}/benchmark_resultados.tsv"

mkdir -p "$logs_dir"

# encabezado de resultados
echo -e "dataset\tdataset_size_MB\tatributo\ttime_real(s)\ttime_user(s)\ttime_sys(s)\tmemory_kb\ttime_DEA(s)" > "$output_table"

# Recorro cada dataset y atributo
for i in "${!datasets[@]}"; do
    dataset="${datasets[$i]}"
    atributo="${atributos[$i]}"

    # Calcular tamaño del dataset en MB
    # 'du -sm' devuelve tamaño en MB, sacamos solo el número con cut
    dataset_size=$(du -sm "datasets/${dataset}" 2>/dev/null | cut -f1)

    # se define path del archivo log (de aca se saca la información de tiempo y memoria)
    log_file="${logs_dir}/${dataset}_${atributo}.log"

    echo "Ejecutando: $dataset - $atributo"

    # Ejecuto AED con Rscript y guardo el log de tiempo y memoria
    # /usr/bin/time se usa para medir el tiempo de ejecución y uso de memoria del script R
    # -f "%e %U %S %M" especifica el formato de salida:
    # %e: tiempo real en segundos, %U: tiempo de CPU en modo usuario, %S: tiempo de CPU en 
    # modo sistema, %M: memoria máxima usada en KB.
    # El resultado se redirige a un archivo de log para su posterior análisis
    /usr/bin/time -f "%e %U %S %M" Rscript "$script" "$dataset" "$atributo" 2> "$log_file"

    # Obtener la última línea del archivo de log (debe contener la salida de /usr/bin/time)
    last_line=$(tail -n 1 "$log_file")

    # Verificar si la última línea comienza con un número
    if [[ $last_line =~ ^[0-9] ]]; then
        # Extraer cada métrica por separado usando awk: tiempo real, tiempo usuario, tiempo sistema y memoria (KB)
        time_real=$(echo "$last_line" | awk '{print $1}')
        time_user=$(echo "$last_line" | awk '{print $2}')
        time_sys=$(echo "$last_line" | awk '{print $3}')
        mem_kb=$(echo "$last_line" | awk '{print $4}')
    else
        # En caso de error al obtener las métricas, imprimir un mensaje de error y asignar "NA" a cada campo
        echo "Error al obtener metricas de tiempo para $dataset - $atributo. Ver log: $log_file"
    fi

    # Buscar dentro del log la lineaque contenga "TIEMPO_LINEA_DIF_EXP:" y extraer el valor posterior a ":"
    # Esta linea la escribe el script R para obtener el tiempo solo de la línea de análisis diferencial
    # Se usa tail -n 1 para obtener la última ocurrencia.
    # cut -d ':' -f2 para obtener el valor despues de los dos puntos.
    # Si no se encuentra la linea, se asigna "NA" para indicar que no hay dato disponible
    # grep busca la linea, tail obtiene la última ocurrencia,
    # cut extrae el valor despues de los dos puntos
    tiempo_linea=$(grep "TIEMPO_LINEA_DIF_EXP:" "$log_file" | tail -n 1 | cut -d ':' -f2)

    # Agreg resultados parseados al archivo de resultados
    echo -e "${dataset}\t${dataset_size}\t${atributo}\t${time_real}\t${time_user}\t${time_sys}\t${mem_kb}\t${tiempo_linea}" >> "$output_table"
done
echo ""

echo "Benchmark terminado."
echo "Resultados guardados en: $output_table"
echo "Logs disponibles en: $logs_dir/"
