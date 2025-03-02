filtrar_datos <- function (datos) {
  library(tidyverse)
  datos <- datos |> 
    select(-c(1,2, 'Fiabilidad'))
  datos$FechaNacimiento <- as.Date(datos$FechaNacimiento)
  datos$edad <- as.numeric(as.Date('2025-02-13')-datos$FechaNacimiento)/365
  datos$Tiempo_Diagnostico <- as.numeric(as.Date('2025-02-28')-datos$FechaDiagnostico)/365
  datos$FechaNacimiento <- NULL
  datos$FechaDiagnostico <- NULL
  datos <- datos[datos$edad>0 & datos$edad<100, ]
  datos <- datos[datos$Tiempo_Diagnostico>0 & datos$Tiempo_Diagnostico<20, ]
  datos$BMI <- datos$Peso_kg/(datos$Altura_cm/100)^2
  nombres <- names(datos)
  for (nombre in nombres){
    columna <- datos[[nombre]]
    NAs <- is.na(columna)
    sum_NAs <- sum(NAs)
    if (sum_NAs < 10 & sum_NAs>0){
      datos <- datos[!NAs, ]
    }
  }
  return (datos)
}
