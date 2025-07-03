library(shiny)
library(shinycssloaders)
library(dplyr)
library(shinybusy)

ui <- fluidPage(
  titlePanel("AED con datasets de CBioPortal"),
  # tableOutput("tabla_resultados") %>% withSpinner(color = "#0dc5c1"),
  add_busy_spinner(
    spin = "fading-circle",    # tipo de spinner
    position = "full-page",    # <--- centrado sobre toda la pantalla
    margins = c(0, 0),
    color = "#0dc5c1",
    timeout = 1000             # tiempo mínimo visible en ms
  ),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Seleccionar un dataset:", choices = NULL),
      selectInput("atributo", "Seleccionar un atributo clínico:", choices = NULL),
      actionButton("ejecutar", "Ejecutar análisis de expresión"),
      h4("Valores únicos del atributo seleccionado:"),
      tableOutput("valores_unicos"),
       width = 6
    ),
    
    mainPanel(
        h4("Tabla de resultados (top 5 por p-valor ajustado):"),
         tableOutput("tabla_resultados"),
        h4("Descargar Resultados (top 50 por p-valor ajustado):"),
        actionButton("generar_reporte", "Generar Reporte"),
        h4("Descargar Resultados (CSV Completo):"),
        downloadButton("descargar_csv", "Descargar CSV completo"),
        width = 6
    )
  )
)

server <- function(input, output, session) {
  datasets_dir <- "datasets"
  
    # observe({
    # carpetas <- list.dirs(datasets_dir, full.names = FALSE, recursive = FALSE)

    # # Filtrar solo las que contienen el archivo data_mrna_seq_v2_rsem.txt
    # carpetas_validas <- Filter(function(carpeta) {
    #     archivo_expr <- file.path(datasets_dir, carpeta, "data_mrna_seq_v2_rsem.txt")
    #     file.exists(archivo_expr)
    # }, carpetas)

    # updateSelectInput(session, "dataset", choices = carpetas_validas)
    # })

    observe({
      carpetas <- list.dirs(datasets_dir, full.names = FALSE, recursive = FALSE)
      
      # Función para calcular tamaño total de una carpeta (en bytes)
      tamano_carpeta <- function(carpeta) {
        ruta_carpeta <- file.path(datasets_dir, carpeta)
        archivos <- list.files(ruta_carpeta, recursive = TRUE, full.names = TRUE)
        if (length(archivos) == 0) return(0)
        sum(file.info(archivos)$size, na.rm = TRUE)
      }
      
      # Filtrar carpetas que contienen el archivo requerido
      carpetas_validas <- Filter(function(carpeta) {
        archivo_expr <- file.path(datasets_dir, carpeta, "data_mrna_seq_v2_rsem.txt")
        file.exists(archivo_expr)
      }, carpetas)
      
      # Calcular tamaño y ordenar carpetas válidas por tamaño descendente
      tamanos <- sapply(carpetas_validas, tamano_carpeta)
      carpetas_ordenadas <- carpetas_validas[order(tamanos, decreasing = TRUE)]
      
      updateSelectInput(session, "dataset", choices = carpetas_ordenadas)
  })
  
  observeEvent(input$dataset, {
    req(input$dataset)
    archivo <- file.path(datasets_dir, input$dataset, "data_clinical_patient.txt")
    
    if (file.exists(archivo)) {
      linea <- readLines(archivo, n = 5)[5]
      atributos <- unlist(strsplit(linea, "\t"))
      atributos <- atributos[atributos != "PATIENT_ID"]
      updateSelectInput(session, "atributo", choices = atributos)
    } else {
      updateSelectInput(session, "atributo", choices = NULL)
    }
  })
  
  output$valores_unicos <- renderTable({
    req(input$dataset, input$atributo)

    archivo <- file.path(datasets_dir, input$dataset, "data_clinical_patient.txt")
    if (!file.exists(archivo)) {
      return(data.frame(Mensaje = "Archivo no encontrado."))
    }

    header_line <- readLines(archivo, n = 5)[5]
    header <- strsplit(header_line, "\t")[[1]]
    
    datos <- read.delim(archivo, header = FALSE, skip = 4, stringsAsFactors = FALSE)
    colnames(datos) <- header

    if (!(input$atributo %in% colnames(datos))) {
      return(data.frame(Mensaje = paste("Atributo", input$atributo, "no encontrado.")))
    }

    # Obtener frecuencias y porcentajes
    valores <- datos[[input$atributo]]
    valores <- valores[-1] # Eliminar la primera fila que es el encabezado
    tabla <- table(trimws(valores))
    porcentajes <- round(100 * prop.table(tabla), 1)

    # Combinar en un data.frame más legible
    data.frame(
      Valor = names(tabla),
      Frecuencia = as.integer(tabla),
      Porcentaje = paste0(porcentajes, " %"),
      check.names = FALSE
    )
  })

  observeEvent(input$ejecutar, {
  req(input$dataset, input$atributo)
  
  # Ejecutar el script externo con los argumentos seleccionados
  script_path <- "aed_limma_cbioportal_TMP.r"
  dataset <- input$dataset
  atributo <- input$atributo
  output_dir <- "output"
  dir.create(output_dir, showWarnings = FALSE)

  # Llamar a Rscript con los argumentos
  system2("Rscript", args = c(script_path, dataset, atributo))
  })

  # Mostrar la tabla con los resultados
  output$tabla_resultados <- renderTable({
    req(input$ejecutar)
    dataset <- input$dataset
    atributo <- input$atributo
    archivo_tabla <- file.path("output", paste0(dataset, "_", atributo, "_top50_results.csv"))

      if (file.exists(archivo_tabla)) {
      # Leer los datos
      tabla <- read.csv(archivo_tabla, nrows = 5)
      
      # Establecer la opción de dígitos para asegurar la visualización correcta
      options(digits = 10)  # Establecer una alta precisión para los números
      
      # Asegurar que los números no sean redondeados al mostrarse
      tabla[] <- lapply(tabla, function(x) if(is.numeric(x)) format(x, scientific = TRUE) else x)
      
      # Mostrar la tabla
      return(tabla)
    } else {
      NULL
    }
  })

  output$descargar_csv <- downloadHandler(
    filename = function() {
      paste0(input$dataset, "_", input$atributo, "_results.csv")
    },
    content = function(file) {
      dataset <- input$dataset
      atributo <- input$atributo
      archivo_csv <- file.path("output", paste0(dataset, "_", atributo, "_results.csv"))
      
      file.copy(archivo_csv, file)
    }
  )

  observeEvent(input$generar_reporte, {
    req(input$dataset, input$atributo)
    
    # Nombre del archivo dinámico
    nombre_archivo <- paste0("informe_resultados_", input$dataset, "_", input$atributo, ".html")
    ruta_salida <- file.path("output", nombre_archivo)
    
    # Generar el informe
    rmarkdown::render(
      input = "report_template.Rmd",
      output_file = ruta_salida,
      params = list(dataset = input$dataset, atributo = input$atributo),
      envir = new.env(parent = globalenv())
    )

    # Abrir el HTML si estás trabajando localmente
    browseURL(ruta_salida)  # Esto abre el archivo generado en el navegador
  })
}

shinyApp(ui = ui, server = server, options = list(port = 8077))
