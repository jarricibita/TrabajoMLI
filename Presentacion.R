library(shiny)
library(readr)
library(caTools)
# Cosas que he añadido: multivariante con HER2, sin HER2, tabla de tipos de cáncer + frec
# He cambiado el sourcing de las funciones

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


# Cross validation pac, sin HER2 ----


#####################################################################################

# SHINY APP
ui <- fluidPage(
  titlePanel("Trabajo final MLI"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Controles Interactivos"),
      selectInput("Modelo", "Elige qué quieres ver:", 
                  choices = c("Distribución de variables", "Bivariante", "Multivariante HER2", "Multivariante sin HER2", "CV + predicciones HER2", "CV + predicciones sin HER2", "Tipos de cáncer")),
      
      conditionalPanel(
        condition = "input.Modelo == 'Distribución de variables'",
        selectInput("Variable", "Selecciona una variable:", choices = NULL)
      )
    ),
    
    mainPanel(
      fluidRow(
        column(8, offset = 2,  
               div(style = "text-align: center;", tableOutput('Tabla')),
               div(style = "text-align: center;", plotOutput("grafico"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  observe({
    updateSelectInput(session, "Variable", choices = names(filtrados))
  })
  
  output$Tabla <- renderTable({
    req(input$Modelo)
    if (input$Modelo == 'Distribución de variables') {
      return(contar_NAs(datos))
    } else if (input$Modelo == 'Bivariante') {
      return(tabla_significativos)
    } else if (input$Modelo == 'Multivariante HER2') {
      return(tabla_significativos_multi[-1, ])
    } else if (input$Modelo == 'Multivariante sin HER2') {
      return(tabla_significativos_multi_pac[-1, ])
    } else if (input$Modelo == 'CV + predicciones') {
      
    } else if (input$Modelo == 'Tipos de cáncer') {
      frec_localizacion <- filtrados |> select(-CRPlevel) |>
        group_by(Localizacion_Primaria) |> summarise(Num = n()) |> mutate(Percentage = Num/(sum(Num))*100)
    }
  })
  
  output$grafico <- renderPlot({
    req(input$Modelo, input$Variable)
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
      
    } else if (input$Modelo == 'CV + predicciones sin HER2') { 
      
    } else if (input$Modelo == 'Tipos de cáncer') {
    }
  })
}

shinyApp(ui = ui, server = server)
