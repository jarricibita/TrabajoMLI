#'@title Cross validation of training portion of dataset
#'
#'@description Tras realizar el modelo multivariante con la parte del dataset completo para training, se requiere la columna dependiente y el dataset restante de variables independientes, la cantidad de batches (k-folds) y threshold para calculo de precision
#'
#'@param training_dataset_crossval_indep
#'@param training_col_crossval_dep
#'@param num_batch
#'@param num_threshold
#'@return list_accuracy, list_pval
#'
#'@example
#'
#'auto_cross_val(num_batch = 4, num_threshold = 0.5)
#'
#'
auto_cross_val <- function(training_dataset_crossval_indep, training_col_crossval_dep, num_batch, num_threshold){
  nombre_dep <- colnames(col_crossval_dep)
  datos_cross_val <- training_dataset_crossval_indep
  datos_cross_val[paste0(nombre_dep)] <- training_col_crossval_dep
  num_batch <- num_batch
  
  num_rows <- nrow(datos_cross_val)
  idx_sampling <- sample(1:num_rows) 
  
  # k-fold, 4 batches
  batches <- list()
  num_idx <- round(((1/num_batch)*length(idx_sampling)))
  for(i in 1:num_batch){
    if(i==num_batch){
      batches[[paste0("idx_batch", i)]] <- idx_sampling
    }else{
      batches[[paste0("idx_batch", i)]] <- idx_sampling[1:num_idx]
      idx_sampling <- idx_sampling[(num_idx+1):length(idx_sampling)]
    }
  }
  
  df_batches_test <- list()
  for(i in 1:num_batch){
    df_batches_test[[paste0("batch", i)]] <- datos_cross_val[batches[[i]], ]
  }
  
  df_batches_train <- list()
  for(i in 1:num_batch){
    df_batches_train[[paste0("batch", i)]] <- datos_cross_val[-batches[[i]], ]
  }
  
  list_pval <- list()
  list_accuracy <- list()
  for(i in 1:num_batch){
    trainingSet_cv <- df_batches_train[[i]]
    testSet_cv <- df_batches_test[[i]]
    # testSet_cv$Remision <- ifelse(testSet_cv$Remision == 1,"Sí","No")
    
    # preparando df
    trainingSet_sin_dependiente <- trainingSet_cv[!paste0(nombre_dep)]
    
    # Bivariante
    print("Lista de variables bivariantes más significativas: ")
    p_value_cv <- table_univar_sig(trainingSet_sin_dependiente, trainingSet_cv[paste0(nombre_dep)])
    colnames_signif_cv <- p_value_cv['Nombres'][p_value_cv['P.Valores']<0.05]
    
    # Primer multivariante
    print("Llevand a cabo modelo multivariante...")
    cancer_model_temp <- glm(trainingSet_cv[paste0(nombre_dep)]~., data = trainingSet_sin_dependiente[colnames_signif_cv], family = "binomial")
    P.valores <- summary(cancer_model_temp)$coefficients[, 4]
    P.valores
    
    # Multivariante rerun
    P.valores <- find_sig_variables(trainingSet_sin_dependiente, trainingSet_cv[paste0(nombre_dep)], umbral = 0.5, P.valores)
    sig_variables <- names(P.valores[-1])
    
    cancer_model_cv <- glm(trainingSet_cv[paste0(nombre_dep)]~., data = trainingSet_sin_dependiente[sig_variables], family = "binomial")
    cancer_model_predict_cv <- predict(cancer_model_cv, testSet_cv, type = "response")
    accuracy_cv <- calc_accuracy(cancer_model_predict_cv, ifelse(trainingSet_cv[paste0(nombre_dep)]==1, "Sí", "No"), threshold = num_threshold)
    list_accuracy[paste0("accuracy", i)] <- accuracy_cv
    list_pval[paste0("p_val", i)] <- list(sig_variables)
    print("########################################################")
  }
  return(list_accuracy, list_pval)
}
