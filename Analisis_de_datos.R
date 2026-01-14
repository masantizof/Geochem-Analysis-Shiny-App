# ============================================================
# ============================================================

# --- 1. Cargar paquetes necesarios ---
paquetes <- c("shiny", "ggplot2", "DBI", "RSQLite", "dplyr", "data.table",
              "DT", "shinycssloaders", "corrplot", "FactoMineR", # "reshape2" eliminado
              "factoextra", "ggrepel", "MASS", "randomForest", "RColorBrewer")

instalar <- paquetes[!paquetes %in% installed.packages()]
if (length(instalar) > 0) install.packages(instalar)

lapply(paquetes, library, character.only = TRUE)

Sys.setenv(LANGUAGE = "en")

# --- 2. Conexión a la base de datos ---

db_path <- "C:/hdm_resultados.sqlite"
con <- dbConnect(RSQLite::SQLite(), db_path)
resultados_db <- tbl(con, "testMASF")


# --- 3. Interfaz de Usuario (UI) ---
ui <- fluidPage(
  titlePanel(" Análisis Espectral y Multivariado de Muestras"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Filtros Geográficos"),
      selectInput("distrito", "Distrito:", choices = NULL, multiple = TRUE),
      selectInput("departamento", "Departamento:", choices = NULL, multiple = TRUE),
      selectInput("municipio", "Municipio:", choices = NULL, multiple = TRUE),
      
      hr(),
      h4("Filtros de Muestra"),
      selectInput("proyecto", "Proyecto:", choices = NULL, multiple = TRUE),
      selectInput("tipo", "Tipo de Muestra:", choices = NULL, multiple = TRUE),
      selectInput("subtipo", "Subtipo Muestra:", choices = NULL, multiple = TRUE),
      selectInput("codigo", "Código Muestra:", choices = NULL, multiple = TRUE),
      selectInput("fechamu", "Fecha de Muestreo:", choices = NULL, multiple = TRUE),
      
      hr(),
      h4("Filtros de Laboratorio"),
      selectInput("laboratorio", "Laboratorio:", choices = NULL, multiple = TRUE),
      selectInput("tipolab", "Tipo Muestra Lab:", choices = NULL, multiple = TRUE),
      selectInput("codigolab", "Código Muestra Lab:", choices = NULL, multiple = TRUE),
      selectInput("orden", "Orden Servicio:", choices = NULL, multiple = TRUE),
      selectInput("numreporte", "Número Reporte GLQ:", choices = NULL, multiple = TRUE),
      
      hr(),
      h4("Análisis Multivariado"),
      selectInput("grupo_elementos", "Grupo de Elementos:",
                  choices = c("Todos", "Platino (PGE)", "Calcófilos", 
                              "Siderófilos", "Litófilos"),
                  selected = "Todos"),
      
      selectInput("variable_color", "Colorear por:",
                  choices = c("tipo_muestra", "proyecto", "municipio", "distrito"),
                  selected = "tipo_muestra"),
      
      checkboxInput("mostrar_etiquetas", "Mostrar etiquetas en PCA", value = FALSE),
      
      checkboxInput("imputar_na", "Imputar valores faltantes (mediana)", value = TRUE),
      
      br(),
      actionButton("run_query", " Generar Análisis", class = "btn-primary"),
      br(), br(),
      htmlOutput("info_datos"),
      br(),
      downloadButton("descargar_csv", " Descargar Resultados"),
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(" Espectro Base", 
                 withSpinner(plotOutput("espectro", height = "600px"), type = 6)),
        
        tabPanel(" Boxplot", 
                 withSpinner(plotOutput("boxplot_plot", height = "600px"), type = 6)),
        
        tabPanel(" Correlación",
                 withSpinner(plotOutput("correlacion_plot", height = "700px"), type = 6),
                 hr(),
                 helpText("Nota: Solo elementos con varianza > 0")),
        
        tabPanel(" PCA - Biplot",
                 withSpinner(plotOutput("pca_biplot", height = "700px"), type = 6),
                 hr(),
                 verbatimTextOutput("pca_summary")),
        
        tabPanel(" PCA - Individuos",
                 withSpinner(plotOutput("pca_individuos", height = "700px"), type = 6)),
        
        tabPanel(" PCA - Variables",
                 withSpinner(plotOutput("pca_variables", height = "700px"), type = 6)),
        
        tabPanel(" LDA - Proyección",
                 withSpinner(plotOutput("lda_plot", height = "600px"), type = 6),
                 hr(),
                 verbatimTextOutput("lda_summary")),
        
        tabPanel(" LDA - Confusión",
                 withSpinner(plotOutput("lda_confusion", height = "500px"), type = 6),
                 hr(),
                 verbatimTextOutput("lda_accuracy")),
        
        tabPanel(" RF - Importancia",
                 withSpinner(plotOutput("rf_importance", height = "700px"), type = 6),
                 hr(),
                 verbatimTextOutput("rf_summary")),
        
        tabPanel(" RF - Confusión",
                 withSpinner(plotOutput("rf_confusion", height = "500px"), type = 6)),
        
        tabPanel(" RF - Error",
                 withSpinner(plotOutput("rf_error", height = "500px"), type = 6)),
        
        tabPanel(" Tabla de Datos",
                 withSpinner(DTOutput("tabla"), type = 6))
      )
    )
  )
)

# --- 4. Servidor (SERVER) ---
server <- function(input, output, session) {

  preparar_datos_multivariado <- function(df_wide, variable_color, imputar_na = FALSE) {

    df_wide_df <- as.data.frame(df_wide)
    
    grupo_var <- as.factor(df_wide_df[[variable_color]])
    
    codigo_muestra <- df_wide_df$codigo_muestra

    cols_numericas <- names(df_wide_df)[sapply(df_wide_df, is.numeric)]
    data_num <- df_wide_df[, cols_numericas, drop = FALSE]
    
    if (imputar_na) {
      for (col in names(data_num)) {
        if (any(is.na(data_num[[col]]))) {
          mediana <- median(data_num[[col]], na.rm = TRUE)
          
          if (!is.na(mediana)) {
            data_num[[col]][is.na(data_num[[col]])] <- mediana
          } else {

            data_num[[col]][is.na(data_num[[col]])] <- 0
          }
        }
      }
    }
    

    filas_completas <- complete.cases(data_num)
    
    data_num_clean <- data_num[filas_completas, , drop = FALSE]
    grupo_var_clean <- grupo_var[filas_completas]
    codigo_muestra_clean <- codigo_muestra[filas_completas]
    
    cols_validas <- sapply(data_num_clean, function(x) {
      all(is.finite(x)) && var(x, na.rm = TRUE) > 1e-10
    })
    
    data_clean <- data_num_clean[, cols_validas, drop = FALSE]
    
    list(
      data = data_clean,
      grupo = grupo_var_clean,
      codigo = codigo_muestra_clean,
      n_original = nrow(df_wide_df),
      n_clean = nrow(data_clean),
      n_vars = ncol(data_clean)
    )
  }
  
  # --- 4.1. Actualización de Filtros ---
  
  observe({
    distritos <- resultados_db %>%
      distinct(distrito) %>%
      arrange(distrito) %>%
      pull()
    updateSelectInput(session, "distrito", choices = distritos)
    
    proyectos <- resultados_db %>% distinct(proyecto) %>% arrange(proyecto) %>% pull()
    updateSelectInput(session, "proyecto", choices = proyectos)
    
    labs <- resultados_db %>% distinct(laboratorio_nombre) %>% arrange(laboratorio_nombre) %>% pull()
    updateSelectInput(session, "laboratorio", choices = labs)
    
    tipos_lab <- resultados_db %>% distinct(tipo_muestra_lab) %>% arrange(tipo_muestra_lab) %>% pull()
    updateSelectInput(session, "tipolab", choices = tipos_lab)
    
    codigos_lab <- resultados_db %>% distinct(codigo_muestra_lab) %>% arrange(codigo_muestra_lab) %>% pull()
    updateSelectInput(session, "codigolab", choices = codigos_lab)
    
    ordenes <- resultados_db %>% distinct(orden_servicio) %>% arrange(orden_servicio) %>% pull()
    updateSelectInput(session, "orden", choices = ordenes)
    
    num_rep <- resultados_db %>% distinct(numero_reporte_glq) %>% arrange(numero_reporte_glq) %>% pull()
    updateSelectInput(session, "numreporte", choices = num_rep)
  })
  
  
  # --- Filtros geográficos ---
  observeEvent(input$distrito, {
    query <- resultados_db
    if (length(input$distrito)) {
      query <- query %>% filter(distrito %in% !!input$distrito)
    }
    deptos <- query %>% distinct(departamento) %>% arrange(departamento) %>% pull()
    updateSelectInput(session, "departamento", choices = deptos)
  })
  
  observeEvent(input$departamento, {
    query <- resultados_db
    if (length(input$distrito)) {
      query <- query %>% filter(distrito %in% !!input$distrito)
    }
    if (length(input$departamento)) {
      query <- query %>% filter(departamento %in% !!input$departamento)
    }
    munis <- query %>% distinct(municipio) %>% arrange(municipio) %>% pull()
    updateSelectInput(session, "municipio", choices = munis)
  })
  
  # --- Filtros de muestra ---
  observeEvent(input$proyecto, {
    query <- resultados_db
    if (length(input$proyecto)) {
      query <- query %>% filter(proyecto %in% !!input$proyecto)
    }
    
    tipos <- query %>% distinct(tipo_muestra) %>% arrange(tipo_muestra) %>% pull()
    updateSelectInput(session, "tipo", choices = tipos)
    
    subtipos <- query %>% distinct(subtipo_muestra) %>% arrange(subtipo_muestra) %>% pull()
    updateSelectInput(session, "subtipo", choices = subtipos)
  })
  
  observeEvent(input$tipo, {
    query <- resultados_db
    if (length(input$proyecto)) query <- query %>% filter(proyecto %in% !!input$proyecto)
    if (length(input$tipo)) query <- query %>% filter(tipo_muestra %in% !!input$tipo)
    
    codigos <- query %>% distinct(codigo_muestra) %>% arrange(codigo_muestra) %>% pull()
    updateSelectInput(session, "codigo", choices = codigos)
    
    fechas <- query %>% distinct(fecha_muestreo) %>% arrange(fecha_muestreo) %>% pull()
    updateSelectInput(session, "fechamu", choices = fechas)
  })
  
  
  # --- 4.2. Reactivos Principales ---
  
  # --- Consulta a la base de datos ---
  datos_filtrados <- eventReactive(input$run_query, {
    query <- resultados_db
    
    # Filtros geográficos
    if (length(input$distrito)) query <- query %>% filter(distrito %in% !!input$distrito)
    if (length(input$departamento)) query <- query %>% filter(departamento %in% !!input$departamento)
    if (length(input$municipio)) query <- query %>% filter(municipio %in% !!input$municipio)
    
    # Filtros de muestra
    if (length(input$proyecto)) query <- query %>% filter(proyecto %in% !!input$proyecto)
    if (length(input$tipo)) query <- query %>% filter(tipo_muestra %in% !!input$tipo)
    if (length(input$subtipo)) query <- query %>% filter(subtipo_muestra %in% !!input$subtipo)
    if (length(input$codigo)) query <- query %>% filter(codigo_muestra %in% !!input$codigo)
    if (length(input$fechamu)) query <- query %>% filter(fecha_muestreo %in% !!input$fechamu)
    
    # Filtros de laboratorio
    if (length(input$laboratorio)) query <- query %>% filter(laboratorio_nombre %in% !!input$laboratorio)
    if (length(input$tipolab)) query <- query %>% filter(tipo_muestra_lab %in% !!input$tipolab)
    if (length(input$codigolab)) query <- query %>% filter(codigo_muestra_lab %in% !!input$codigolab)
    if (length(input$orden)) query <- query %>% filter(orden_servicio %in% !!input$orden)
    if (length(input$numreporte)) query <- query %>% filter(numero_reporte_glq %in% !!input$numreporte)
    
    cat("Ejecutando consulta SQL...\n")
    
    df_final <- query %>% collect()
    
    cat(paste(nrow(df_final), "filas recuperadas.\n"))
    
    setDT(df_final)
    df_final[, valor_ppm := as.numeric(valor_ppm)]
    df_final <- df_final[!is.na(valor_ppm) & !is.na(elemento)]
    setorder(df_final, elemento)
    
    return(df_final)
  })
  
  # --- Preparar datos en formato ancho ---
  datos_wide <- reactive({
    df_f <- datos_filtrados()
    validate(need(nrow(df_f) > 0, "No hay datos en la base para esta selección."))
    
    df_wide <- dcast(df_f, codigo_muestra + tipo_muestra + subtipo_muestra + 
                       proyecto + municipio + departamento + distrito ~ elemento, 
                     value.var = "valor_ppm", fun.aggregate = mean, na.rm = TRUE)

    if (input$grupo_elementos != "Todos") {
      elementos_patron <- switch(input$grupo_elementos,
                                 "Platino (PGE)" = c("Ru", "Rh", "Pd", "Ir", "Pt"),
                                 "Calcófilos" = c("Cu", "Zn", "Ga", "Ge", "As", "Se", "Ag", "Cd", 
                                                  "In", "Sn", "Sb", "Te", "Hg", "Tl", "Pb", "Bi","Au"),
                                 "Siderófilos" = c("Fe", "Co", "Ni", "Ru", "Rh", "Pd", "Ir", "Pt","Au","Re"),
                                 "Litófilos" = c("Mg", "Al", "Si", "Ti", "V", "Cr", "Mn", "Zr", 
                                                 "Nb", "Ba", "La", "Ta", "W", "Th", "U","Ce")
      )
      
      cols_categoricas <- c("codigo_muestra", "tipo_muestra", "subtipo_muestra",
                            "proyecto", "municipio", "departamento", "distrito")
      cols_numericas <- names(df_wide)[sapply(df_wide, is.numeric)]
      cols_mantener <- cols_categoricas
      
      for (elem in elementos_patron) {

        patron <- paste0("^", elem, "\\d*$|", elem)
        cols_mantener <- c(cols_mantener, 
                           grep(patron, cols_numericas, value = TRUE, ignore.case = TRUE))
      }
      
      df_wide <- df_wide[, unique(cols_mantener), with = FALSE]
    }
    
    validate(need(nrow(df_wide) > 0, "No hay muestras con los elementos seleccionados."))
    
    return(df_wide)
  })
  
  # --- 4.3. Salidas ---
  
  output$info_datos <- renderUI({
    df_f <- datos_filtrados()
    df_wide <- datos_wide()
    
    n_muestras <- nrow(df_wide)
    n_elementos <- sum(sapply(df_wide, is.numeric))
    
    df_wide_df <- as.data.frame(df_wide)
    cols_num <- names(df_wide_df)[sapply(df_wide_df, is.numeric)]
    data_num <- df_wide_df[, cols_num]
    
    filas_completas <- sum(complete.cases(data_num))
    pct_completas <- if (n_muestras > 0) round(filas_completas / n_muestras * 100, 1) else 0
    
    grupo_var <- df_wide_df[[input$variable_color]]
    tabla_grupos <- table(grupo_var)
    
    HTML(paste0(
      "<div style='background-color: #f0f0f0; padding: 10px; border-radius: 5px; font-size: 12px;'>",
      "<b> Resumen de Datos:</b><br>",
      "• Muestras: ", n_muestras, "<br>",
      "• Elementos: ", n_elementos, "<br>",
      "• Filas completas: ", filas_completas, " (", pct_completas, "%)<br>",
      "<br><b>Grupos (", input$variable_color, "):</b><br>",
      paste("• ", names(tabla_grupos), ": ", tabla_grupos, collapse = "<br>"),
      "</div>"
    ))
  })
  

  output$espectro <- renderPlot({
    df_f <- datos_filtrados()
    validate(need(nrow(df_f) > 0, "No hay datos para esta selección."))
    
    df_f[, ELE_NUM := as.numeric(factor(elemento, levels = sort(unique(elemento))))]
    
    ggplot(df_f, aes(x = ELE_NUM, y = valor_ppm, 
                     group = codigo_muestra, color = codigo_muestra)) +
      geom_area(aes(fill = codigo_muestra), alpha = 0.3, position = "identity") +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      scale_x_continuous(breaks = df_f$ELE_NUM, labels = df_f$elemento) +
      labs(x = "Elemento", y = "Concentración (ppm)", 
           color = "Código muestra", title = "Espectro elemental por muestra") +
      theme_minimal(base_size = 15) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
            legend.position = "bottom")
  })
  
  # --- Boxplot ---
  output$boxplot_plot <- renderPlot({
    df_f <- datos_filtrados()
    validate(need(nrow(df_f) > 0, "No hay datos para esta selección."))
    
    medianas <- df_f %>%
      group_by(elemento) %>%
      summarise(value = median(valor_ppm, na.rm = TRUE))
    
    ggplot(df_f, aes(x = elemento, y = valor_ppm)) +
      geom_boxplot(fill = "#F7DC6F", outlier.shape = 1, outlier.size = 0.5) +
      scale_y_log10() +
      geom_line(data = medianas, aes(x = elemento, y = value, group = 1),
                color = "darkgreen", linewidth = 0.5) +
      geom_point(data = medianas, aes(x = elemento, y = value),
                 color = "grey", size = 1.5) +
      theme_light(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Elemento", y = "Valor en log(ppm)",
           title = "Distribución de valores por elemento con línea de medianas")
  })
  
  # --- Correlación ---
  output$correlacion_plot <- renderPlot({

    df_wide_dt <- datos_wide()
    validate(need(nrow(df_wide_dt) > 0, "No hay datos para esta selección."))
    
    # Convertir a data.frame
    df_wide <- as.data.frame(df_wide_dt)
    
    cols_numericas <- names(df_wide)[sapply(df_wide, is.numeric)]
    data_num <- df_wide[, cols_numericas]
    
    varianzas <- sapply(data_num, function(x) var(x, na.rm = TRUE))
    cols_con_varianza <- names(varianzas[varianzas > 0 & !is.na(varianzas)])
    data_filtered <- data_num[, cols_con_varianza]
    
    validate(need(ncol(data_filtered) >= 2, 
                  "No hay suficientes elementos con varianza."))
    
    corr_matrix <- cor(data_filtered, use = "pairwise.complete.obs")
    
    corrplot(corr_matrix, method = 'square', type = "upper", 
             order = "original", tl.cex = 0.9, tl.col = "black", 
             cl.cex = 0.9, title = "Matriz de Correlación entre Elementos",
             mar = c(0, 0, 2, 0))
  })
  
  # ============================================================
  # PCA
  # ============================================================
  pca_result <- reactive({
    df_wide <- datos_wide()
    validate(need(nrow(df_wide) >= 3, "Se necesitan al menos 3 muestras para PCA."))
    
    df_wide_df <- as.data.frame(df_wide)
    
    cols_numericas <- names(df_wide_df)[sapply(df_wide_df, is.numeric)]
    data_num <- df_wide_df[, cols_numericas]
    
    varianzas <- sapply(data_num, function(x) var(x, na.rm = TRUE))
    cols_con_varianza <- names(varianzas[varianzas > 0 & !is.na(varianzas)])
    data_filtered <- data_num[, cols_con_varianza]
    
    validate(need(ncol(data_filtered) >= 2, "Insuficientes variables para PCA."))
    
    rownames(data_filtered) <- df_wide_df$codigo_muestra
    
    pca <- PCA(log(data_filtered + 1), scale.unit = TRUE, 
               ncp = min(5, ncol(data_filtered)-1), graph = FALSE)
    
    list(pca = pca, df_wide = df_wide_df, data_filtered = data_filtered)
  })
  
  # --- PCA Biplot ---
  output$pca_biplot <- renderPlot({
    res <- pca_result()
    df_wide <- res$df_wide
    
    var_color <- as.factor(df_wide[[input$variable_color]])
    
    p <- fviz_pca_biplot(res$pca,
                         geom.ind = "point",
                         pointshape = 16, 
                         pointsize = 4,
                         addEllipses = TRUE,
                         mean.point = FALSE,
                         circle = TRUE,
                         label = "var",
                         repel = TRUE,
                         title = paste("PCA Biplot -", input$grupo_elementos),
                         col.ind = var_color,
                         palette = "Set2")
    
    if (input$mostrar_etiquetas) {
      coord_ind <- as.data.frame(res$pca$ind$coord)
      coord_ind$etiquetas <- df_wide$codigo_muestra
      
      p <- p + ggrepel::geom_text_repel(data = coord_ind,
                                        aes(x = Dim.1, y = Dim.2, label = etiquetas),
                                        size = 2.5, max.overlaps = 20)
    }
    
    p + theme(legend.text = element_text(size = 12),
              legend.title = element_blank())
  })
  
  # --- PCA Individuos ---
  output$pca_individuos <- renderPlot({
    res <- pca_result()
    df_wide <- res$df_wide
    
    var_color <- as.factor(df_wide[[input$variable_color]])
    
    fviz_pca_ind(res$pca,
                 geom.ind = "point",
                 pointshape = 16,
                 pointsize = 4,
                 addEllipses = TRUE,
                 mean.point = FALSE,
                 title = paste("PCA - Individuos -", input$grupo_elementos),
                 col.ind = var_color,
                 palette = "Set2",
                 legend.title = input$variable_color) +
      theme(legend.text = element_text(size = 12))
  })
  
  # --- PCA Variables ---
  output$pca_variables <- renderPlot({
    res <- pca_result()
    
    fviz_pca_var(res$pca,
                 col.var = "contrib",
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE,
                 title = paste("PCA - Variables -", input$grupo_elementos)) +
      theme_minimal(base_size = 14)
  })
  
  output$pca_summary <- renderPrint({
    res <- pca_result()
    cat("=== RESUMEN PCA ===\n\n")
    cat("Varianza explicada por dimensión:\n")
    print(round(res$pca$eig, 2))
  })
  
  # ============================================================
  # LDA
  # ============================================================
  lda_result <- reactive({
    df_wide <- datos_wide()
    validate(need(nrow(df_wide) >= 5, "Se necesitan al menos 5 muestras para LDA."))
    
    datos_prep <- preparar_datos_multivariado(df_wide, input$variable_color, input$imputar_na)
    
    validate(need(datos_prep$n_vars >= 2, "Insuficientes variables válidas para LDA."))
    validate(need(datos_prep$n_clean >= 5, "Insuficientes observaciones completas para LDA."))
    validate(need(length(unique(datos_prep$grupo)) >= 2, "Se necesitan al menos 2 grupos."))
    
    tabla_grupos <- table(datos_prep$grupo)
    grupos_pequenos <- names(tabla_grupos[tabla_grupos < 2])
    validate(need(length(grupos_pequenos) == 0, 
                  paste("Grupos con muy pocas observaciones:", 
                        paste(grupos_pequenos, collapse = ", "))))
    
    data_lda <- datos_prep$data
    data_lda$grupo <- datos_prep$grupo
    
    tryCatch({
      lda_model <- lda(grupo ~ ., data = data_lda)
      pred <- predict(lda_model)
      
      # Validación cruzada
      lda_cv <- lda(grupo ~ ., data = data_lda, CV = TRUE)
      
      list(
        model = lda_model, 
        pred = pred, 
        data = data_lda,
        cv = lda_cv, 
        codigo = datos_prep$codigo,
        n_removed = datos_prep$n_original - datos_prep$n_clean
      )
    }, error = function(e) {
      validate(need(FALSE, paste("Error en LDA (revise si hay colinealidad):", e$message)))
    })
  })
  
  output$lda_plot <- renderPlot({
    res <- lda_result()
    
    df_plot <- as.data.frame(res$pred$x)
    df_plot$grupo <- res$data$grupo
    df_plot$etiquetas <- res$codigo  
    
    n_dims <- ncol(res$pred$x)
    
    if(n_dims >= 2) {
      var_exp1 <- round(res$model$svd[1]^2 / sum(res$model$svd^2) * 100, 1)
      var_exp2 <- round(res$model$svd[2]^2 / sum(res$model$svd^2) * 100, 1)
      
      ggplot(df_plot, aes(x = LD1, y = LD2, color = grupo)) +
        geom_point(size = 3, alpha = 0.7) +
        stat_ellipse(level = 0.95, linewidth = 1) +
        theme_minimal(base_size = 14) +
        labs(
          title = paste("Análisis Discriminante Lineal (LDA) -", input$grupo_elementos),
          subtitle = if(res$n_removed > 0) 
            paste("(", res$n_removed, "muestras removidas por datos faltantes)") else NULL,
          x = paste0("LD1 (", var_exp1, "%)"),
          y = paste0("LD2 (", var_exp2, "%)")
        ) +
        theme(legend.position = "bottom", legend.title = element_blank())
    } else {
      ggplot(df_plot, aes(x = LD1, y = 0, color = grupo)) +
        geom_point(size = 3, position = position_jitter(height = 0.1), alpha = 0.7) +
        theme_minimal(base_size = 14) +
        labs(
          title = "Análisis Discriminante Lineal (LDA) - 1D",
          subtitle = if(res$n_removed > 0) 
            paste("(", res$n_removed, "muestras removidas por datos faltantes)") else NULL,
          x = "LD1", y = ""
        ) +
        theme(legend.position = "bottom", legend.title = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
  })
  
  output$lda_summary <- renderPrint({
    res <- lda_result()
    
    cat("=== RESUMEN LDA ===\n\n")
    print(res$model)
    
    cat("\n\nProporción de varianza explicada:\n")
    var_exp <- res$model$svd^2 / sum(res$model$svd^2)
    print(round(var_exp, 3))
  })
  
  output$lda_confusion <- renderPlot({
    res <- lda_result()
    
    tabla <- table(Predicho = res$pred$class, Real = res$data$grupo)
    conf_df <- as.data.frame(tabla)
    
    ggplot(conf_df, aes(x = Real, y = Predicho, fill = Freq)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Freq), color = "black", size = 6) +
      scale_fill_gradient(low = "white", high = "steelblue") +
      theme_minimal(base_size = 14) +
      labs(title = "Matriz de Confusión - LDA")
  })
  
  output$lda_accuracy <- renderPrint({
    res <- lda_result()
    
    cat("=== EXACTITUD LDA ===\n\n")
    
    tabla <- table(Predicho = res$pred$class, Real = res$data$grupo)
    accuracy <- sum(diag(tabla)) / sum(tabla)
    cat(paste("Exactitud global:", round(accuracy * 100, 2), "%\n"))
    
    aer <- 1 - accuracy
    cat(paste("Tasa de error aparente:", round(aer * 100, 2), "%\n\n"))
    
    # Validación cruzada
    tabla_cv <- table(Predicho = res$cv$class, Real = res$data$grupo)
    accuracy_cv <- sum(diag(tabla_cv)) / sum(tabla_cv)
    cat(paste("Exactitud con validación cruzada:", round(accuracy_cv * 100, 2), "%\n"))
  })
  
  # ============================================================
  # RANDOM FOREST
  # ============================================================
  rf_result <- reactive({
    df_wide <- datos_wide()
    validate(need(nrow(df_wide) >= 5, "Se necesitan al menos 5 muestras para RF."))
    
    datos_prep <- preparar_datos_multivariado(df_wide, input$variable_color, input$imputar_na)
    
    validate(need(datos_prep$n_vars >= 2, "Insuficientes variables válidas para RF."))
    validate(need(datos_prep$n_clean >= 5, "Insuficientes observaciones completas para RF."))
    validate(need(length(unique(datos_prep$grupo)) >= 2, "Se necesitan al menos 2 grupos."))
    
    tabla_grupos <- table(datos_prep$grupo)
    grupos_pequenos <- names(tabla_grupos[tabla_grupos < 2])
    validate(need(length(grupos_pequenos) == 0, 
                  paste("Grupos con muy pocas observaciones:", 
                        paste(grupos_pequenos, collapse = ", "))))
    

    data_rf <- datos_prep$data
    data_rf$grupo <- datos_prep$grupo
    
    tryCatch({
      set.seed(123)
      rf_model <- randomForest(
        grupo ~ ., 
        data = data_rf, 
        importance = TRUE, 
        ntree = 500,
        na.action = na.fail 
      )
      
      list(
        model = rf_model, 
        data = data_rf,
        codigo = datos_prep$codigo,
        n_removed = datos_prep$n_original - datos_prep$n_clean
      )
    }, error = function(e) {
      validate(need(FALSE, paste("Error en RF:", e$message)))
    })
  })
  
  output$rf_importance <- renderPlot({
    res <- rf_result()
    
    validate(need(!is.null(res$model$importance), 
                  "No se pudo calcular la importancia de las variables."))
    
    imp_df <- as.data.frame(importance(res$model))
    
    validate(need("MeanDecreaseAccuracy" %in% names(imp_df),
                  "No se encontró MeanDecreaseAccuracy en importance."))
    
    # Remover NAs
    imp_df <- imp_df[is.finite(imp_df$MeanDecreaseAccuracy), ]
    
    validate(need(nrow(imp_df) > 0, "No hay datos de importancia disponibles."))
    
    imp_df$Variable <- rownames(imp_df)
    imp_df <- imp_df[order(-imp_df$MeanDecreaseAccuracy), ]
    
    n_top <- min(20, nrow(imp_df))
    imp_df <- head(imp_df, n_top)
    
    ggplot(imp_df, aes(x = reorder(Variable, MeanDecreaseAccuracy), 
                       y = MeanDecreaseAccuracy)) +
      geom_col(fill = "steelblue", alpha = 0.8) +
      coord_flip() +
      labs(
        title = paste("Importancia de Variables - Random Forest"),
        subtitle = paste("Top", n_top, "de", nrow(importance(res$model)), "variables"),
        x = NULL, 
        y = "Mean Decrease Accuracy"
      ) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.major.y = element_blank())
  })
  
  output$rf_summary <- renderPrint({
    res <- rf_result()
    
    cat("=== RESUMEN RANDOM FOREST ===\n\n")
    cat(paste("Tipo de modelo:", res$model$type, "\n"))
    cat(paste("Número de árboles:", res$model$ntree, "\n"))
    cat(paste("Variables por split (mtry):", res$model$mtry, "\n"))
    cat(paste("Muestras utilizadas:", nrow(res$data), "\n"))
    if(res$n_removed > 0) {
      cat(paste("Muestras removidas por NAs:", res$n_removed, "\n"))
    }
    cat("\n")
    
    cat("=== MATRIZ DE CONFUSIÓN OOB ===\n")
    conf <- res$model$confusion[, 1:(ncol(res$model$confusion)-1)]
    print(conf)
    
    cat("\n=== EXACTITUD POR CLASE (OOB) ===\n")
    for(i in 1:nrow(conf)) {
      clase <- rownames(conf)[i]
      total <- sum(conf[i, ])
      correctos <- conf[i, i]
      acc <- correctos / total * 100
      cat(sprintf("%s: %.2f%% (%d/%d)\n", clase, acc, correctos, total))
    }
    
    cat("\n=== EXACTITUD GLOBAL (OOB) ===\n")
    accuracy <- sum(diag(conf)) / sum(conf)
    cat(sprintf("Exactitud: %.2f%%\n", accuracy * 100))
    cat(sprintf("Error OOB: %.2f%%\n", (1 - accuracy) * 100))
  })
  
  output$rf_confusion <- renderPlot({
    res <- rf_result()
    
    conf <- as.table(res$model$confusion[, 1:nlevels(res$data$grupo)])
    conf_df <- as.data.frame(conf)
    names(conf_df) <- c("Real", "Predicho", "Frecuencia")
    
    ggplot(conf_df, aes(x = Real, y = Predicho, fill = Frecuencia)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Frecuencia), color = "black", size = 5) +
      scale_fill_gradient(low = "white", high = "steelblue") +
      theme_minimal(base_size = 14) +
      labs(title = "Matriz de Confusión - Random Forest (OOB)")
  })
  
  output$rf_error <- renderPlot({
    res <- rf_result()
    
    validate(need(!is.null(res$model$err.rate), 
                  "No se encontró información de error en el modelo."))
    
    error_matrix <- res$model$err.rate
    
    validate(need(nrow(error_matrix) > 0 && ncol(error_matrix) > 0,
                  "Matriz de error vacía."))

    error_df <- data.frame(
      arboles = 1:nrow(error_matrix),
      OOB = error_matrix[, "OOB"]
    )
    
    clases <- setdiff(colnames(error_matrix), "OOB")
    for (clase in clases) {
      error_df[[clase]] <- error_matrix[, clase]
    }
    

    error_long <- data.table::melt(setDT(error_df), id.vars = "arboles", 
                                   variable.name = "Tipo", value.name = "Error")
    
    error_long <- error_long[is.finite(Error)]
    
    validate(need(nrow(error_long) > 0, 
                  "No hay datos de error válidos para graficar."))
    
    n_clases <- length(clases)
    if(n_clases <= 8) {
      colores <- c("OOB" = "black", 
                   setNames(brewer.pal(max(3, n_clases), "Set2")[1:n_clases], clases))
    } else {
      colores <- c("OOB" = "black", 
                   setNames(rainbow(n_clases), clases))
    }
    
    ggplot(error_long, aes(x = arboles, y = Error, color = Tipo)) +
      geom_line(linewidth = 1, alpha = 0.8) +
      scale_color_manual(values = colores) +
      theme_minimal(base_size = 14) +
      labs(
        title = "Evolución del Error en Random Forest",
        subtitle = paste("Basado en", res$model$ntree, "árboles"),
        x = "Número de Árboles",
        y = "Tasa de Error",
        color = "Tipo de Error"
      ) +
      theme(legend.position = "bottom") +
      scale_y_continuous(labels = scales::percent_format())
  })
  
  # --- 4.4. Tabla y Descarga ---
  
  output$tabla <- renderDT({
    df_f <- datos_filtrados()
    datatable(df_f, 
              options = list(
                pageLength = 15, 
                scrollX = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
              ),
              rownames = FALSE,
              filter = 'top')
  })
  
  output$descargar_csv <- downloadHandler(
    filename = function() {
      paste0("resultados_", input$grupo_elementos, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- datos_filtrados()
      fwrite(df, file, sep = ",")
    }
  )
  
  onStop(function() {
    dbDisconnect(con)
    cat("Conexión a la base de datos cerrada.\n")
  })
}

# --- 5. Lanzar aplicación ---
shinyApp(ui = ui, server = server)