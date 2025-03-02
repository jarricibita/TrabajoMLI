#'@title Otros calculos de dataset, edad y BMI...
#'
#'@description Calcula la edad dada una columna de nacimiento y BMI dadas peso y altura, y las añade al dataset
#'
#'@param datos Dataset
#'
#'@return Resumen + característica según tipo
#'


calc_edad <- function(dataset, col_fecha){
  fecha_actual <- Sys.Date()
  dataset[paste0(colnames(col_fecha))] <- as.Date(col_fecha)
  dataset[paste0("Edad")] <- as.numeric(round((Sys.Date() - as.Date(col_fecha))/365))
  dataset <- dataset[!(is.na(dataset[paste0("Edad")])), ]
  if(sum((dataset[paste0("Edad")]<0) | (dataset[paste0("Edad")]>100))>0){
    dataset <- dataset[!(dataset[paste0("Edad")]<0 | dataset[paste0("Edad")]>100), ]
  }
  print(paste0("Nuevas dimensiones del dataset con edad: ", dim(dataset)[1], " ", dim(dataset)[2]))
  return(dataset)
}


calc_bmi <- function(dataset, col_peso, col_altura){
  dataset[paste0("BMI")] <- round(col_peso/(col_altura/100), 2)
  print(paste0("Nuevas dimensiones del dataset con BMI: ", dim(dataset)[1], " ", dim(dataset)[2]))
  return(dataset)
}

calc_otrasfechas <- function(dataset, col_fecha){
  fecha_actual <- Sys.Date()
  dataset[paste0(colnames(col_fecha))] <- as.Date(col_fecha)
  dataset[paste0("FechaAlgo")] <- as.numeric(round((Sys.Date() - as.Date(col_fecha))/365))
  dataset <- dataset[!(is.na(dataset[paste0("FechaAlgo")])), ]
  if(sum((dataset[paste0("FechaAlgo")]<0) | (dataset[paste0("FechaAlgo")]>100))>0){
    dataset <- dataset[!(dataset[paste0("FechaAlgo")]<0 | dataset[paste0("FechaAlgo")]>100), ]
  }
  print(paste0("Nuevas dimensiones del dataset con edad: ", dim(dataset)[1], " ", dim(dataset)[2]))
  return(dataset)
}
