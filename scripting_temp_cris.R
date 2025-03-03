library(caTools)
# Cosas que he añadido: multivariante con HER2, sin HER2, tabla de tipos de cáncer + frec
# He cambiado el sourcing de las funciones

# Funciona estando en una carpeta igual a la del repositorio
source('./contar_NAs.R')
source('./escalar.R')
source('./filtrar_datos.R')
# source('./bivariante.R')
source("./src/fun_resumen_dataset.R")
source("./src/fun_dummy.R")
source("./src/fun_normalizacion.R")
source("./src/fun_univar_pvalue.R")
source("./src/fun_rerun_multivariant.R")
source("./src/fun_accuracy.R")


# Tranformaciones y funciones
datos <- read_csv("pacientes_cancer3.csv")
filtrados <- filtrar_datos(datos)
datos_finales <- nuevo_dataset_dummy(filtrados)

# P-valores bivariante
tabla_pvalores_univar <- table_univar_sig(datos_finales |> select(-Remision), datos_finales$Remision)
tabla_significativos <- tabla_pvalores_univar[tabla_pvalores_univar['P.Valores']<0.05, ]
colnames_signif <- tabla_pvalores_univar['Nombres'][tabla_pvalores_univar['P.Valores']<0.05]

# Datos con col de HER2 y modelo multivariante ----
## Quitar col y NA
datos_HER2 <- datos_finales |> select(-CRPlevel) |> remove_na()
## Modelo logistico
m.logistico <- glm(datos_HER2$Remision~., data = datos_HER2[colnames_signif], family = "binomial")
P.valores_multi <- summary(m.logistico)$coefficients[, 4]
P.valores_multi <- find_sig_variables(datos_HER2 |> select(-Remision), datos_HER2$Remision, umbral = 0.75, P.valores_multi)
sig_variables_HER2 <- names(P.valores_multi[-1])
## Tabla con Nombres y P-valor significativos (con HER2)
tabla_significativos_multi <- as.data.frame(P.valores_multi[-1]) |> reframe(P.valor = P.valores_multi[-1]) |> mutate(Nombres = names(P.valores_multi[-1])) |> select(Nombres, P.valor)

# Datos SIN col de HER2, con más pacientes y modelo multivariante ----
## Quitar col y NA
datos_pac <- datos_finales |> select(-c(CRPlevel, Mut_HER2)) |> remove_na() 
## Columnas significativas sin HER2
colnames_signif_pac <- colnames_signif[!colnames_signif%in%"Mut_HER2"]
## Modelo logístico
m.logistico_pac <- glm(datos_pac$Remision~., data = datos_pac[colnames_signif_pac], family = "binomial")
P.valores_multi_pac <- summary(m.logistico_pac)$coefficients[, 4]
P.valores_multi_pac <- find_sig_variables(datos_pac |> select(-Remision), datos_pac$Remision, umbral = 0.75, P.valores_multi_pac)
sig_variables_pac <- names(P.valores_multi_pac[-1])
## Tabla con Nombres y P-valor significativos (sin HER2)
tabla_significativos_multi_pac <- as.data.frame(P.valores_multi_pac[-1]) |> reframe(P.valor = P.valores_multi_pac[-1]) |> mutate(Nombres = names(P.valores_multi_pac[-1])) |> select(Nombres, P.valor)


# Model training early ----

prueba <- train_and_test_first(datos_pac, datos_pac['Remision'], size_training = 0.8, sig_variable_names = sig_variables_pac, threshold = 0.5)



# prueba cross val ----
auto_cross_val <- function(training_dataset_crossval, training_col_crossval_dep, num_batch, num_threshold){
  nombre_dep <- colnames(training_col_crossval_dep)
  datos_cross_val <- training_dataset_crossval
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
    cancer_model_temp <- glm(trainingSet_cv[[paste0(nombre_dep)]]~., data = trainingSet_sin_dependiente[colnames_signif_cv], family = "binomial")
    P.valores <- summary(cancer_model_temp)$coefficients[, 4]
    P.valores
    
    # Multivariante rerun
    P.valores <- find_sig_variables(trainingSet_sin_dependiente, trainingSet_cv[paste0(nombre_dep)], umbral = 0.5, P.valores)
    sig_variables <- names(P.valores[-1])
    
    cancer_model_cv <- glm(trainingSet_cv[[paste0(nombre_dep)]]~., data = trainingSet_sin_dependiente[sig_variables], family = "binomial")
    cancer_model_predict_cv <- predict(cancer_model_cv, testSet_cv, type = "response")
    accuracy_cv <- calc_accuracy(cancer_model_predict_cv, ifelse(trainingSet_cv[paste0(nombre_dep)]==1, "Sí", "No"), threshold = num_threshold)
    list_accuracy[paste0("accuracy", i)] <- accuracy_cv
    list_pval[paste0("p_val", i)] <- list(sig_variables)
    print("########################################################")
  }
  return(list(list_accuracy, list_pval))
}

test <- prueba[[2]]
train <- prueba[[1]]
auto_cross_val()