#'@title Ejecuta modelización multivariante y cálculo de precisión con una lista de varios datasets
#'
#'@description En el caso del trabajo, era con divisiones de datasets por tipo de localización del cáncer. Los datasets se encuentran divididos antes de meterlos en esta función en una lista, cada elemento de la lista con un nombre identificativo.
#'
#'@param list_of_datasets: list of datasets with names
#'@param name_dep_variable: str
#'@param accuracy_threshold: threshold for accuracy calculation for each 
#'@param chosen_seed: seed for training testing split
#'
#'@return Lista de listas: precisiones, predicciones (probabilidad), verdaderos; se puede cambiar
#'
#'@example
#'
#'estudio_datasets_multiples(accuracy_threshold = 0.5)
#'

estudio_datasets_multiples <- function(list_of_datasets, name_dep_variable, accuracy_threshold, chosen_seed = 123){
  lista_precisiones_tipo <- NULL
  lista_predicciones_tipo <- NULL
  lista_verdaderos_tipo <- NULL
  for(i in 1:length(list_of_datasets)){
    nombre_dataset <- names(list_of_datasets)[i]
    print(nombre_dataset)
    dataset <- list_of_datasets[[i]]
    nombre_dep <- name_dep_variable
    
    tabla_pvalores_univar <- table_univar_sig(dataset[, !names(dataset)%in%nombre_dep], dataset[[nombre_dep]])
    sig_variables_bi <- tabla_pvalores_univar['Nombres'][tabla_pvalores_univar['P.Valores']<0.05]
    
    model_multi <- glm(dataset[[nombre_dep]]~., data = dataset[sig_variables_bi], family = "binomial")
    P.valores_multi <- summary(model_multi)$coefficients[, 4]
    print(P.valores_multi)
    P.valores_multi <- find_sig_variables(dataset[, !names(dataset)%in%nombre_dep], dataset[[nombre_dep]], umbral = 0.75, P.valores_multi)
    sig_variables_multi <- names(P.valores_multi[-1])
    
    tabla_significativos_multi <- as.data.frame(P.valores_multi[-1]) |> reframe(P.valor = P.valores_multi[-1]) |> mutate(Nombres = names(P.valores_multi[-1])) |> select(Nombres, P.valor)
    
    dataset_modelo <- dataset[sig_variables_multi]
    dataset_modelo[nombre_dep] <- dataset[nombre_dep]
    modelo_inicial <- train_and_test_first(dataset_modelo, dataset_modelo[nombre_dep], size_training = 0.8, sig_variable_names = sig_variables_multi, threshold = accuracy_threshold, seed = chosen_seed)
    train_set <- modelo_inicial[[1]]
    test_set <- modelo_inicial[[2]]
    accuracy_ini <- modelo_inicial[[3]]
    predictions <- modelo_inicial[[4]]
    TrueValues <- modelo_inicial[[5]]
    print("################################################################")
    lista_precisiones_tipo[paste0("Precisión_Dataset_", nombre_dataset)] <- accuracy_ini
    lista_predicciones_tipo[paste0("Predicciones_Dataset_", nombre_dataset)] <- list(predictions)
    lista_verdaderos_tipo[paste0("Verdaderos_Dataset_", nombre_dataset)] <- list(TrueValues)
  }
  return(list(lista_precisiones_tipo, lista_predicciones_tipo, lista_verdaderos_tipo))
}