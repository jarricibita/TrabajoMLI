contar_NAs <- function(datos){
  library(readr)
  library(tidyverse)
  nombres <- names(datos)
  NAs <- numeric(length(nombres))
  i <- 1
  for (nombre in nombres){
    columna <- datos[[nombre]]
    NAs[i] <- sum(is.na(columna))
    i <- i+1
  }
  Tabla_NAs <- data.frame(nombres, NAs)
  Tabla_NAs |> 
    filter(NAs>0)
} 

