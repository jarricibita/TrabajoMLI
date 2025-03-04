#'@title Busqueda de p-valores significativos en modelo univariante
#'
#'@description Crea modelo logístico entre variable dependiente y todas las variables dependientes, muestra un resumen del modelo y saca una tabla con cada nombre de variable y p-valor
#'
#'@param dataset
#'
#'@return Dataset categórico a dummy
#'
#'@example
#'
#'table_univar_sig(dataset_independent = c(0, 1), dependent_variable = c(1, 2))
#'
#'
table_univar_sig <- function(dataset_independent, dependent_variable){
  Nombres <- NULL
  P.Valores <- NULL
  y <- dependent_variable
  for (i in 1:length(dataset_independent)){
    x <- dataset_independent[[i]]
    if(class(x) == "character"){
      x <- as.factor(x)
    }
    m.logistico <- glm(y ~ x, family = "binomial")
    # print(paste0(colnames(dataset_independent[i])))
    # print(summary(m.logistico))
    # print("####################################################################")
    if(!is.null(summary(m.logistico)$coefficients[2, ])){
    # if(nrow(summary(m.logistico)$coefficients)>1){  
      p.valor <- summary(m.logistico)$coefficients[2, 4]
    }else{
      p.valor <- NULL
    }
    Nombres <- c(Nombres, colnames(dataset_independent[i]))
    P.Valores <- c(P.Valores, p.valor)
  }
  tabla_pvalor <- data.frame(Nombres, P.Valores)
  # print(tabla_pvalor)
  return(tabla_pvalor)
}

