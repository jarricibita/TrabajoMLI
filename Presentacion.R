library(shiny)
library(readr)
library(caTools)
# Cosas que he añadido: multivariante con HER2, sin HER2, tabla de tipos de cáncer + frec
# He cambiado el sourcing de las funciones
# He añadido un return de Remision reales de cada crossvalidation
# FALTA: añadir tablas? otros resultados? en la seccion de cross validations, no sé qué poner
# REMOVE COLUMNS LEFT FOR LOCALIZACION PRIMARIA

# Funciona estando en una carpeta igual a la del repositorio
source('./contar_NAs.R')
source('./escalar.R')
source('./filtrar_datos.R')
# source('./bivariante.R')
source("./src/fun_resumen_dataset.R")
source("./src/fun_dummy.R")
source("./src/fun_normalizacion.R")
source("./src/fun_univar_pvalue.R")
source("./src/fun_rerun_multivariant.R")
source("./src/fun_accuracy.R")
source("./src/fun_training_test_first.R")
source("./src/fun_crossval.R")

# Tranformaciones y funciones
datos <- read_csv("pacientes_cancer3.csv")
filtrados <- filtrar_datos(datos)
datos_finales <- nuevo_dataset_dummy(filtrados)

# P-valores bivariante
tabla_pvalores_univar <- table_univar_sig(datos_finales |> select(-Remision), datos_finales$Remision)
tabla_significativos <- tabla_pvalores_univar[tabla_pvalores_univar['P.Valores']<0.05, ]
colnames_signif <- tabla_pvalores_univar['Nombres'][tabla_pvalores_univar['P.Valores']<0.05]

# Datos con col de HER2 y modelo multivariante ----
## Quitar col y NA
datos_HER2 <- datos_finales |> select(-CRPlevel) |> remove_na()
## Modelo logistico
m.logistico <- glm(datos_HER2$Remision~., data = datos_HER2[colnames_signif], family = "binomial")
P.valores_multi <- summary(m.logistico)$coefficients[, 4]
P.valores_multi <- find_sig_variables(datos_HER2 |> select(-Remision), datos_HER2$Remision, umbral = 0.75, P.valores_multi)
sig_variables_HER2 <- names(P.valores_multi[-1])
## Tabla con Nombres y P-valor significativos (con HER2)
tabla_significativos_multi <- as.data.frame(P.valores_multi[-1]) |> reframe(P.valor = P.valores_multi[-1]) |> mutate(Nombres = names(P.valores_multi[-1])) |> select(Nombres, P.valor)

# Datos SIN col de HER2, con más pacientes y modelo multivariante ----
## Quitar col y NA
datos_pac <- datos_finales |> select(-c(CRPlevel, Mut_HER2)) |> remove_na() 
## Columnas significativas sin HER2
colnames_signif_pac <- colnames_signif[!colnames_signif%in%"Mut_HER2"]
## Modelo logístico
m.logistico_pac <- glm(datos_pac$Remision~., data = datos_pac[colnames_signif_pac], family = "binomial")
P.valores_multi_pac <- summary(m.logistico_pac)$coefficients[, 4]
P.valores_multi_pac <- find_sig_variables(datos_pac |> select(-Remision), datos_pac$Remision, umbral = 0.75, P.valores_multi_pac)
sig_variables_pac <- names(P.valores_multi_pac[-1])
## Tabla con Nombres y P-valor significativos (sin HER2)
tabla_significativos_multi_pac <- as.data.frame(P.valores_multi_pac[-1]) |> reframe(P.valor = P.valores_multi_pac[-1]) |> mutate(Nombres = names(P.valores_multi_pac[-1])) |> select(Nombres, P.valor)

# Cross validation con HER2 ----
datos_HER2_modelo <- datos_HER2[sig_variables_HER2]
datos_HER2_modelo$Remision <- datos_HER2$Remision
modelo_inicial_HER2 <- train_and_test_first(datos_HER2_modelo, datos_HER2_modelo['Remision'], size_training = 0.8, sig_variable_names = sig_variables_HER2, threshold = 0.5)
train_HER2 <- modelo_inicial_HER2[[1]]
test_HER2 <- modelo_inicial_HER2[[2]]
accuracy_ini_HER2 <- modelo_inicial_HER2[[3]]
predictions_HER2 <- modelo_inicial_HER2[[4]]
TrueValues_HER2 <- modelo_inicial_HER2[[5]]

# Cross validation: Lista de precisiones, variables y predicciones
return_crossvalHER2 <- auto_cross_val(train_HER2, train_HER2['Remision'], num_batch = 4, num_threshold = 0.5)
precision_crossvalHER2 <- return_crossvalHER2[[1]]
variables_crossvalHER2 <- return_crossvalHER2[[2]]
predicciones_crossvalHER2 <- return_crossvalHER2[[3]]
realidad_crossvalHER2 <- return_crossvalHER2[[4]] # Remision REAL
tabla_HER2 <- as.data.frame(sort(table(unlist(return_crossvalHER2[[2]])), decreasing=TRUE)) |> 
  rename('Variable' = Var1,
         'Splits' = Freq)



# Cross validation pac, sin HER2 ----
datos_pac_modelo <- datos_pac[sig_variables_pac]
datos_pac_modelo$Remision <- datos_pac$Remision
modelo_inicial_pac <- train_and_test_first(datos_pac_modelo, datos_pac_modelo['Remision'], size_training = 0.8, sig_variable_names = sig_variables_pac, threshold = 0.5)
train_pac <- modelo_inicial_pac[[1]]
test_pac <- modelo_inicial_pac[[2]]
accuracy_ini_pac <- modelo_inicial_pac[[3]] # Precision inicial
predictions_pac <- modelo_inicial_pac[[4]]
TrueValues_pac <- modelo_inicial_pac[[5]]

# Cross validation: Lista de precisiones, variables y predicciones
return_crossvalpac <- auto_cross_val(train_pac, train_pac['Remision'], num_batch = 4, num_threshold = 0.5)
precision_crossvalpac <- return_crossvalpac[[1]]
variables_crossvalpac <- return_crossvalpac[[2]]
predicciones_crossvalpac <- return_crossvalpac[[3]]
realidad_crossvalpac <- return_crossvalpac[[4]] # Remision REAL
tabla_pac <- as.data.frame(sort(table(unlist(return_crossvalpac[[2]])), decreasing=TRUE)) |> 
  rename('Variable' = Var1,
         'Splits' = Freq)


# Modelos por tipo de cáncer
datos_prostata <- datos_HER2 |> filter(Próstata_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))
datos_colon <- datos_HER2 |> filter(Colon_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))
datos_pulmon <- datos_HER2 |> filter(Pulmón_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))
datos_linfoma <- datos_HER2 |> filter(Linfoma_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))
datos_leucemia <- datos_HER2 |> filter(Leucemia_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))
datos_mama <- datos_HER2 |> filter(Mama_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))
datos_melanoma <- datos_HER2 |> filter(Melanoma_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))

lista_datos_localizacion <- list(datos_prostata, datos_colon, datos_pulmon, datos_linfoma, datos_mama, datos_melanoma)
names(lista_datos_localizacion) <- c("Próstata", "Colon", "Pulmón", "Linfoma", "Mama", "Melanoma")
lista_precisiones_tipo <- NULL


#####################################################################################

# SHINY APP ----
ui <- fluidPage(
  titlePanel("Trabajo final MLI"),
  h4('Cristina Li y Javier Arricibita'),
  
  sidebarLayout(
    sidebarPanel(
      h3("Controles Interactivos"),
      selectInput("Modelo", "Elige qué quieres ver:", 
                  choices = c("Distribución de variables", "Bivariante", "Multivariante", "CV + predicciones HER2", "CV + predicciones sin HER2", "Tipos de cáncer")),
      
      conditionalPanel(
        condition = "input.Modelo == 'Distribución de variables'",
        selectInput("Variable", "Selecciona una variable:", choices = NULL)
      ),
      
      conditionalPanel(
        condition = "input.Modelo == 'CV + predicciones HER2'",
        selectInput("NumCrossVal", "Selecciona una predicción:", choices = NULL)
      ),
      
      conditionalPanel(
        condition = "input.Modelo == 'CV + predicciones sin HER2'",
        selectInput("NumCrossVal2", "Selecciona una predicción:", choices = NULL)
      )
    ),
    
    mainPanel(
      fluidRow(
        column(8, offset = 2,  
               div(style = "text-align: center;", uiOutput('Tabla')),
               div(style = "text-align: center;", uiOutput("TextoPrecision")),
               div(style = "text-align: center;", plotOutput("grafico"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  observe({
    updateSelectInput(session, "Variable", choices = names(filtrados))
    updateSelectInput(session, "NumCrossVal", choices = c('General', names(predicciones_crossvalHER2)))
    updateSelectInput(session, "NumCrossVal2", choices = c('General', names(predicciones_crossvalpac)))
  })
  
  output$Tabla <- renderUI({
    req(input$Modelo)
    
    if (input$Modelo == 'Distribución de variables') {
      tableOutput("tablaDistribucion")
      
    } else if (input$Modelo == 'Bivariante') {
      tableOutput("tablaBivariante")
      
    } else if (input$Modelo == 'Multivariante') {
      tagList(
        h4("Multivariante HER2"),
        tableOutput("tablaHER2"),
        paste("Precisión:", accuracy_ini_HER2, "%"),
        h4("Multivariante sin HER2"),
        tableOutput("tablaSinHER2"),
        paste("Precisión: ", accuracy_ini_pac, "%")
      )
      
    } else if (input$Modelo == 'CV + predicciones HER2'){
      tableOutput('tablaCVHER2')
      
    } else if (input$Modelo == 'CV + predicciones sin HER2'){
      tableOutput('tablaCVsinHER2')
    
    } else if (input$Modelo == 'Tipos de cáncer') {
      tagList(
        h4("HER2 mantenida"),
        tableOutput("tablaTiposCancerHER2"),
        h4("HER2 sin mantener"),
        tableOutput("tablaTiposCancer")
      )
    }
  })
  
  # Función para formatear tablas
  formatear_tabla <- function(tabla) {
    as.data.frame(lapply(tabla, function(x) {
      if (is.numeric(x)) {
        format(x, digits = 4)  # Formato numérico
      } else {
        x
      }
    }))
  }
  
  # Renderizado de las tablas con formato
  output$tablaDistribucion <- renderTable({
    formatear_tabla(contar_NAs(datos))
  })
  
  output$tablaBivariante <- renderTable({
    formatear_tabla(tabla_significativos)
  })
  
  output$tablaHER2 <- renderTable({
    formatear_tabla(tabla_significativos_multi)
  })
  
  output$tablaSinHER2 <- renderTable({
    formatear_tabla(tabla_significativos_multi_pac)
  })
  
  output$tablaCVHER2 <- renderTable({
    req(input$NumCrossVal)  
    
    # Según lo que seleccione el usuario en el desplegable, muestra el vector correspondiente
    if (input$NumCrossVal == 'General'){
      formatear_tabla(tabla_HER2)
    } else if(input$NumCrossVal == "Prediction1") {
      variables_crossvalHER2$p_val1
    } else if (input$NumCrossVal == "Prediction2") {
      variables_crossvalHER2$p_val2  
    } else if (input$NumCrossVal == "Prediction3") {
      variables_crossvalHER2$p_val3  
    } else if (input$NumCrossVal == "Prediction4") {
      variables_crossvalHER2$p_val4  
    }
  })
  
  
  output$tablaCVsinHER2 <- renderTable({
    req(input$NumCrossVal2)  
    
    # Según lo que seleccione el usuario en el desplegable, muestra el vector correspondiente
    if (input$NumCrossVal2 == 'General'){
      formatear_tabla(tabla_pac)
    } else if(input$NumCrossVal2 == "Prediction1") {
      variables_crossvalpac$p_val1
    } else if (input$NumCrossVal2 == "Prediction2") {
      variables_crossvalpac$p_val2  
    } else if (input$NumCrossVal2 == "Prediction3") {
      variables_crossvalpac$p_val3  
    } else if (input$NumCrossVal2 == "Prediction4") {
      variables_crossvalpac$p_val4  
    }
  })
  
  output$tablaTiposCancerHER2 <- renderTable({
    formatear_tabla(
      filtrados |> select(-CRPlevel) |> 
        group_by(Localizacion_Primaria) |> drop_na() |> summarise(Num = n()) |> mutate(Percentage = Num / sum(Num) * 100)
    )
  })
  
  output$tablaTiposCancer <- renderTable({
    formatear_tabla(
      filtrados |> select(-CRPlevel) |>
        group_by(Localizacion_Primaria) |> summarise(Num = n()) |> mutate(Percentage = Num / sum(Num) * 100)
    )
  })
  
  # Precisión en las pestañas de CV
  output$TextoPrecision <- renderUI({
    req(input$Modelo)
    if (input$Modelo == 'CV + predicciones HER2') {
      textOutput('PrecisionesHER2')
    } else if (input$Modelo == 'CV + predicciones sin HER2') {
      textOutput('Precisionespac')
    } else {
      NULL
    }
  })
  
  output$PrecisionesHER2 <- renderText({
    req(input$NumCrossVal)
    if (input$NumCrossVal == 'General'){
      paste("Precisión:", accuracy_ini_HER2, "%")
    } else if(input$NumCrossVal == "Prediction1") {
      paste("Precisión:", precision_crossvalHER2$accuracy1, "%")
    } else if (input$NumCrossVal == "Prediction2") {
      paste("Precisión:", precision_crossvalHER2$accuracy2, "%")
    } else if (input$NumCrossVal == "Prediction3") {
      paste("Precisión:", precision_crossvalHER2$accuracy3, "%")
    } else if (input$NumCrossVal == "Prediction4") {
      paste("Precisión:", precision_crossvalHER2$accuracy4, "%")
    }
  })
  
  output$Precisionespac <- renderText({
    req(input$NumCrossVal2)
    if (input$NumCrossVal2 == 'General'){
      paste("Precisión:", accuracy_ini_pac, "%")
    } else if(input$NumCrossVal2 == "Prediction1") {
      paste("Precisión:", precision_crossvalpac$accuracy1, "%")
    } else if (input$NumCrossVal2 == "Prediction2") {
      paste("Precisión:", precision_crossvalpac$accuracy2, "%")
    } else if (input$NumCrossVal2 == "Prediction3") {
      paste("Precisión:", precision_crossvalpac$accuracy3, "%")
    } else if (input$NumCrossVal2 == "Prediction4") {
      paste("Precisión:", precision_crossvalpac$accuracy4, "%")
    }
  })
  
  output$grafico <- renderPlot({
    req(input$Modelo, input$Variable, input$NumCrossVal, input$NumCrossVal2)
    if (input$Modelo == 'Distribución de variables') {
      var <- filtrados[[input$Variable]]
      if (is.numeric(var)) {
        hist(var, main = paste("Histograma de:", input$Variable), col = "lightblue", border = "black")
      } else {
        freqs <- table(var, exclude = 'NA')
        par(mar = c(10, 10, 4, 7))
        barplot(freqs, main = paste('Diagrama de barras de', input$Variable), names.arg = names(freqs), horiz = TRUE, las=2)
      }
    } else if (input$Modelo == 'CV + predicciones HER2') {
      num_crossval <- predicciones_crossvalHER2[[input$NumCrossVal]]
      #hist(num_crossval, main = paste("Frecuencias probabilidad de", input$NumCrossVal), xlim=c(0, 1), xlab="Probabilidad predictiva")
      colores <- c("red", "blue") 
      colores_HER2 <- colores[as.integer(TrueValues_HER2)+1]
      plot(predictions_HER2, col=colores_HER2, pch=16, main='Probabilidades')
      legend('topright', 
             legend = c('Remisión 0', 'Remisión 1'),
             col = c('red', 'blue'),
             pch=16,
             box.lwd = 0)
    } else if (input$Modelo == 'CV + predicciones sin HER2') { 
      num_crossval <- predicciones_crossvalpac[[input$NumCrossVal2]]
      colores <- c("red", "blue") 
      colores_sinHER2 <- colores[as.integer(TrueValues_pac)+1]
      plot(predictions_pac, col=colores_sinHER2, pch=16, main='Probabilidades')
      legend('topright', 
             legend = c('Remisión 0', 'Remisión 1'),
             col = c('red', 'blue'),
             pch=16,
             box.lwd=0)
    } else if (input$Modelo == 'Tipos de cáncer') {
    }
  })
}

shinyApp(ui = ui, server = server)
