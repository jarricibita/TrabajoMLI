#'@title Modelización repetida con las columnas más significativas hasta llegar a un modelo cuyas variables son todas significativas
#'
#'@description Cuando el modelo multivariante no cuenta con una lista completa de variables significativas, se repite el modelo con menos variables, escogiendo la cantidad de variables acorde a un umbral (0-1) hasta que todas las variables son significativas.
#'
#'@param dataset
#'@param umbral
#'@return p-valores de variables más significativas
#'
#'@example
#'
#'table_univar_sig(umbral = 0.5)
#'
#'
find_sig_variables <- function(dataset_indep, var_dep, umbral, p.values){
  while(sum(p.values[-1]>0.05)>0){
    num_colnames <- round(length(p.values[-1])*umbral)
    colnames_sorted <- p.values[-1] |> sort() |> names()
    colnames_keep <- colnames_sorted[1:num_colnames]
    
    if(is.null(colnames_keep)) {
      print("No se han encontrado variables significativas en el modelo multivariante.")
      break}
    
    m.logistico <- glm(var_dep~., data = dataset_indep[colnames_keep], family = "binomial")
    p.values <- summary(m.logistico)$coefficients[, 4]
    
  }
  print("Lista de variables más significativas:")
  print(p.values[-1])
  return(p.values)
}