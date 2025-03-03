#'@title Calculo de precisión (accuracy)
#'
#'@description Dado predicción y datos de testeo
#'
#'@param model_prediction
#'@param testSet_and_dep
#'@return accuracy
#'
calc_accuracy <- function(model_prediction, testSet_and_dep, threshold){
  clasificacion <- ifelse(model_prediction > threshold,"Sí","No")
  clasificacion <- ordered(clasificacion, levels = c("Sí", "No"))
  testSet_and_dep <- ordered(testSet_and_dep, levels = c("Sí", "No"))
  
  mat_res <- table(Predicted = clasificacion, Actual = testSet_and_dep)
  mat_res
  accuracy <- (mat_res[1, 1]+mat_res[2, 2])/(sum(mat_res))
  accuracy <- round(100*accuracy, 2)
  print(paste0("Precisión del modelo: ", accuracy, '%'))
  return(accuracy)
}
