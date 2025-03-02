#'@title Resumen de un dataset...
#'
#'@description Da el tipo así como otras características de cada columna, quitar NA y contar NA
#'
#'@param datos Dataset
#'
#'@return Resumen + característica según tipo
#'
#'@example
#'
#'resumen_dataset(datos = 1)
#'

is.convertible.to.date <- function(x){
  !is.na(as.Date(as.character(x), tz = 'UTC', format = '%Y-%m-%d'))}

columna_numerica <- function(x){
  print("Resumen numérico de la columna:")
  print(paste0("El mínimo es ", min(x, na.rm = T), ", el máximo es ", max(x, na.rm = T), "."))
  print(paste0("La media es: ", mean(x, na.rm = T)))
  print(paste0("Cantidad de NA: ", sum(is.na(x))))
}

columna_caracter <- function(x){
  print("Elementos de la columna:")
  print(unique(x))
}

columna_fechas <- function(x){
  print("Columna de fechas con rango:")
  print(range(x))
  print(paste0("Cantidad de NA: ", sum(is.na(x))))
  
}

graficas_categoricas <- function(x, num, dataset){
  barplot(prop.table(table(x)), ylab = "Proportion", main = colnames(dataset[num]))
}

graficas_numericas <- function(x, num, dataset){
  hist(x, ylab = "Frequency", xlab = colnames(dataset[num]), main = paste0("Distribution of ", colnames(dataset[num]), " ", num))
}

calc_edad <- function(dataset){
  fecha_actual <- Sys.Date()
  dataset$Edad <- round((Sys.Date() - as.Date(dataset$FechaNacimiento))/365)
}

resumen_dataset <- function(datos){
  lista_numericos <- NULL
  lista_fechas <- NULL
  lista_categoricos <- NULL
  print("La base de datos tiene la siguientes columnas:")
  for(i in 1:length(datos)){
    if((is.convertible.to.date(datos[[i]][[1]])) == TRUE){
      datos[[i]] <- as.Date(datos[[i]])
    }
    tipo = class(datos[[i]])
    print(paste("Columna", i, "-", colnames(datos[i]), ":", tipo, sep = " "))
    if(tipo=="integer"){
      columna_numerica(datos[[i]])
      lista_numericos <- c(lista_numericos, i)
    }
    if(tipo=="numeric"){
      columna_numerica(datos[[i]])
      lista_numericos <- c(lista_numericos, i)
    }
    else if(tipo=="character"){
      columna_caracter(datos[[i]])
      lista_categoricos <- c(lista_categoricos, i)
    }
    else if(tipo=="Date"){
      columna_fechas(datos[[i]])
      lista_fechas <- c(lista_fechas, i)
    }
    
    cat("#######################################################################################\n")
  }
  for(j in lista_categoricos){
    valores_categoricos <- datos[[j]]
    graficas_categoricas(valores_categoricos, j, datos)
  }
  for(k in lista_numericos){
    valores_numericos <- datos[[k]]
    graficas_numericas(valores_numericos, k, datos)
  }
}

remove_na <- function(dataset){
  library(tidyverse)
  dataset <- dataset |> drop_na()
  print(paste0("Nuevas dimensiones del dataset: ", dim(dataset)[1], " ", dim(dataset)[2]))
  return(dataset)
}

count_na <- function(dataset){
  df_na <- tibble(Name = colnames(dataset), NA_count = sapply(dataset, function(x) sum(is.na(x)))) |> 
    subset(NA_count>0)
  return(df_na)
}
