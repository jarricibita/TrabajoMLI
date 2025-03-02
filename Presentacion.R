library(shiny)
library(readr)

source('./contar_NAs.R')
source('./escalar.R')
source('./filtrar_datos.R')
source('./bivariante.R')
source('./nuevo_dataset_dummy.R')

datos <- read_csv("pacientes_cancer3.csv")
filtrados <- filtrar_datos(datos)
datos_finales <- nuevo_dataset_dummy(filtrados)

ui <- fluidPage(
  titlePanel("Trabajo final MLI"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Controles Interactivos"),
      selectInput("Modelo", "Elige qué quieres ver:", 
                  choices = c('Distribución de variables', "Bivariante", "Multivariante", "CV + predicciones", 'Tipos de cáncer')),
      
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
      return(bivariante(datos_finales))
    } else if (input$Modelo == 'Multivariante') {
      
    } else if (input$Modelo == 'CV + predicciones') {
      
    } else if (input$Modelo == 'Tipos de cáncer') {
      
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
    } else if (input$Modelo == 'CV + predicciones') {
      
    } else if (input$Modelo == 'Tipos de cáncer') {
      
    }
  })
}

shinyApp(ui = ui, server = server)
