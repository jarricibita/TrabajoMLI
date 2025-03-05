#'@title Generando variables dummy para columnas categóricas...
#'
#'@description Es para el dataset de hipertensión que empieza con columna de fecha que se ignora...
#' Para cada columna que es categórica, cuenta las categorías y según el tipo de categoría, crea variables dummy de un tipo u otro
#'
#'@param dataset
#'
#'@return Dataset categórico a dummy
#'
#'@example
#'
#'nuevo_dataset_dummy(dataset = 1)
#'
#'
nuevo_dataset_dummy <- function(dataset){
  nuevo_dataset <- dataset
  for(i in 1:length(dataset)){
    if(class(dataset[[i]])=="character"){
      factores <- unique(dataset[[i]])
      num_factores <- unique(dataset[[i]]) |> table() |> sum()
      
      if(num_factores==2){
        factores <- unique(dataset[[i]])
        if("Sí" %in% factores){
          col_temp <- dataset[[i]]
          col_temp <- as.numeric(col_temp == "Sí")
          nuevo_dataset[paste0(colnames(dataset[i]))] <- col_temp
        }else{
          col_temp <- dataset[[i]]
          col_temp <- as.numeric(col_temp == paste0(factores[1]))
          nuevo_dataset[paste0(colnames(dataset[i]))] <- col_temp
        }
      }
      
      if(num_factores>2){
        for(j in 1:length(factores)){
          col_temp <- dataset[[i]]
          col_temp <- as.numeric(col_temp == paste0(factores[j]))
          nuevo_dataset[paste0(factores[j], "_",colnames(dataset[i]))] <- col_temp
        }
        nuevo_dataset <- nuevo_dataset[, !(colnames(nuevo_dataset) %in% colnames(dataset)[i])]
      }
    }
  }
  print(paste0("Nuevas dimensiones del dataset con dummies: ", dim(nuevo_dataset)[1], " ", dim(nuevo_dataset)[2]))
  print("####################################################################")
  print("Nombres de las nuevas columnas: ")
  print(colnames(nuevo_dataset))
  return(nuevo_dataset)
}