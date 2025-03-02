library(shiny)
library(readr)

datos <- read_csv("pacientes_cancer3.csv")
filtrados <- filtrar_datos(datos)
datos_HER2 <- filtrados |> select(-CRPlevel) |> remove_na()

# Dummies y escalado
datos_HER2dummy <- nuevo_dataset_dummy(datos_HER2)
datos_HER2escalado <- nuevo_dataset_normalizado(datos_HER2dummy)

# Bivariante HER2
datos_HER2 <- datos_HER2escalado
tabla_pvalores_univar <- table_univar_sig(datos_HER2 |> select(-Remision), datos_HER2$Remision)
tabla_pvalores_univar[tabla_pvalores_univar['P.Valores']<0.05, ]
colnames_signif <- tabla_pvalores_univar['Nombres'][tabla_pvalores_univar['P.Valores']<0.05]

ui <- fluidPage(
  titlePanel("Trabajo final MLI"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Controles Interactivos"),
      selectInput("Modelo", "Elige qué quieres ver:", 
                  choices = c('Distribución de variables', "Bivariante", "Multivariante", "CV + predicciones", 'Tipos de cáncer')),
      
      conditionalPanel(
        condition = "input.Modelo == 'Distribución de variables'",
        selectInput("Variable", "Selecciona una variable:", choices = NULL))
# AQUI NO SÉ CÓMO PONER LO DE LAS PESTAÑAS DE SHINY      
      conditionalPanel(
        condition = "input.Modelo == 'Bivariante'",
        selectInput("Visualización", "Selecciona opción", choices = NULL))
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
    }
  })
}

# Ejecutar la aplicación
shinyApp(ui = ui, server = server)

