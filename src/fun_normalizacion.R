#'@title Normalización de dataset...
#'
#'@description Pasa por cada columna numérica realiza un Shapiro test
#' y aquellos significativos son transformados a log y escalados
#'
#'@param datos dataset
#'
#'@return Dataset normalizado
#'
#'@example
#'
#'nuevo_dataset_normalizado(dataset = 1)
#'
#'

# Escalar variable, z
normalizacion_escalado <- function(x) {
  temp <- log(x + 1)
  mean_temp <- mean(temp, na.rm = TRUE)
  sd_temp <- sd(temp, na.rm = TRUE)
  z <- (temp - mean_temp) / sd_temp
  return(z)
}

nuevo_dataset_normalizado <- function(dataset){
  nuevo_dataset <- dataset
  for(i in 1:length(dataset)){
    if(class(dataset[[i]])=="integer"){
      p.valor <- shapiro.test(dataset[[i]])$p.value
      if(p.valor<0.05){
        nuevo_dataset[paste0(colnames(dataset[i]), "_scale")] <- normalizacion_escalado(dataset[[i]])
        nuevo_dataset <- nuevo_dataset[, !(colnames(nuevo_dataset) %in% colnames(dataset)[i])]}
    }
    if(class(dataset[[i]])=="numeric"){
      p.valor <- shapiro.test(dataset[[i]])$p.value
      if(p.valor<0.05){
        nuevo_dataset[paste0(colnames(dataset[i]), "_scale")] <- normalizacion_escalado(dataset[[i]])
        nuevo_dataset <- nuevo_dataset[, !(colnames(nuevo_dataset) %in% colnames(dataset)[i])]}
    }
  }
  return(nuevo_dataset)
}