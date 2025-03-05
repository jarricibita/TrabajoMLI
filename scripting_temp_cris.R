library(readr)
library(caTools)
library(plotROC)
# Cosas que he añadido: multivariante con HER2, sin HER2, tabla de tipos de cáncer + frec
# He cambiado el sourcing de las funciones
# He añadido un return de Remision reales de cada crossvalidation
# FALTA: añadir tablas? otros resultados? en la seccion de cross validations, no sé qué poner
# REMOVE COLUMNS LEFT FOR LOCALIZACION PRIMARIA

# Funciona estando en una carpeta igual a la del repositorio
source('./contar_NAs.R')
source('./escalar.R')
source('./filtrar_datos.R')
# source('./bivariante.R')
source("./src/fun_resumen_dataset.R")
source("./src/fun_normalizacion.R")
source("./src/fun_dummy.R")
source("./src/fun_normalizacion.R")
source("./src/fun_univar_pvalue.R")
source("./src/fun_rerun_multivariant.R")
source("./src/fun_accuracy.R")
source("./src/fun_training_test_first.R")
source("./src/fun_crossval.R")
source("./src/fun_estudio_datasets_multiples.R")

# Tranformaciones y funciones
datos <- read_csv("pacientes_cancer3.csv")
filtrados <- filtrar_datos(datos)
datos_normalizados <- nuevo_dataset_normalizado(filtrados |> select(-c(CRPlevel, Remision)))
datos_normalizados$Remision <- filtrados$Remision
datos_finales <- nuevo_dataset_dummy(datos_normalizados)
datos_finales$CRPlevel <- filtrados$CRPlevel

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

# Cross validation con HER2 ----
datos_HER2_modelo <- datos_HER2[sig_variables_HER2]
datos_HER2_modelo$Remision <- datos_HER2$Remision
modelo_inicial_HER2 <- train_and_test_first(datos_HER2_modelo, datos_HER2_modelo['Remision'], size_training = 0.8, sig_variable_names = sig_variables_HER2, threshold = 0.5)
train_HER2 <- modelo_inicial_HER2[[1]]
test_HER2 <- modelo_inicial_HER2[[2]]
accuracy_ini_HER2 <- modelo_inicial_HER2[[3]]
predictions_HER2 <- modelo_inicial_HER2[[4]]
TrueValues_HER2 <- modelo_inicial_HER2[[5]]

# Cross validation: Lista de precisiones, variables y predicciones
return_crossvalHER2 <- auto_cross_val(train_HER2, train_HER2['Remision'], num_batch = 4, num_threshold = 0.5)
precision_crossvalHER2 <- return_crossvalHER2[[1]]
variables_crossvalHER2 <- return_crossvalHER2[[2]]
predicciones_crossvalHER2 <- return_crossvalHER2[[3]]
realidad_crossvalHER2 <- return_crossvalHER2[[4]] # Remision REAL
tabla_HER2 <- as.data.frame(sort(table(unlist(return_crossvalHER2[[2]])), decreasing=TRUE)) |> 
  rename('Variable' = Var1,
         'Splits' = Freq)



# Cross validation pac, sin HER2 ----
datos_pac_modelo <- datos_pac[sig_variables_pac]
datos_pac_modelo$Remision <- datos_pac$Remision
modelo_inicial_pac <- train_and_test_first(datos_pac_modelo, datos_pac_modelo['Remision'], size_training = 0.8, sig_variable_names = sig_variables_pac, threshold = 0.5)
train_pac <- modelo_inicial_pac[[1]]
test_pac <- modelo_inicial_pac[[2]]
accuracy_ini_pac <- modelo_inicial_pac[[3]] # Precision inicial
predictions_pac <- modelo_inicial_pac[[4]]
TrueValues_pac <- modelo_inicial_pac[[5]]

# Cross validation: Lista de precisiones, variables y predicciones
return_crossvalpac <- auto_cross_val(train_pac, train_pac['Remision'], num_batch = 4, num_threshold = 0.5)
precision_crossvalpac <- return_crossvalpac[[1]]
variables_crossvalpac <- return_crossvalpac[[2]]
predicciones_crossvalpac <- return_crossvalpac[[3]]
realidad_crossvalpac <- return_crossvalpac[[4]] # Remision REAL
tabla_pac <- as.data.frame(sort(table(unlist(return_crossvalpac[[2]])), decreasing=TRUE)) |> 
  rename('Variable' = Var1,
         'Splits' = Freq)


# Modelos por tipo de cáncer
datos_prostata <- datos_HER2 |> filter(Próstata_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))
datos_colon <- datos_HER2 |> filter(Colon_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))
datos_pulmon <- datos_HER2 |> filter(Pulmón_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))
datos_linfoma <- datos_HER2 |> filter(Linfoma_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))
datos_leucemia <- datos_HER2 |> filter(Leucemia_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))
datos_mama <- datos_HER2 |> filter(Mama_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))
datos_melanoma <- datos_HER2 |> filter(Melanoma_Localizacion_Primaria == "1") |> select(-contains("Localizacion_Primaria"))

lista_datos_localizacion <- list(datos_prostata, datos_colon, datos_pulmon, datos_linfoma, datos_mama, datos_melanoma)
names(lista_datos_localizacion) <- c("Próstata", "Colon", "Pulmón", "Linfoma", "Mama", "Melanoma")
nombre_dep <- 'Remision'

return_tipos_1 <- estudio_datasets_multiples(lista_datos_localizacion, nombre_dep, accuracy_threshold = 0.5, chosen_seed = 123)
return_tipos_2 <- estudio_datasets_multiples(lista_datos_localizacion, nombre_dep, accuracy_threshold = 0.5, chosen_seed = 1818)
return_tipos_3 <- estudio_datasets_multiples(lista_datos_localizacion, nombre_dep, accuracy_threshold = 0.5, chosen_seed = 8787)

return_tipos_1[1]

medias_precision_tipo <- mapply(function(x, y, z) mean(c(x, y, z)), return_tipos_1[[2]], return_tipos_2[[2]], return_tipos_3[[2]])


ifelse(testSet[[paste0(name_dep)]]==1, "Sí", "No")


lista_thresholds <- c(0.3, 0.4, 0.5, 0.6, 0.7)
tabla_precisiones_tipo <- data.frame(Tipo_Cancer = names(lista_datos_localizacion))
for (i in 1:length(lista_thresholds)){
  lista_medias_temp <- NULL
  for (j in 1:length(lista_datos_localizacion)){
    prec1 <- calc_accuracy(return_tipos_1[[3]][[j]], ifelse(return_tipos_1[[4]][[j]]==1, "Sí", "No"), threshold = lista_thresholds[i])
    prec2 <- calc_accuracy(return_tipos_2[[3]][[j]], ifelse(return_tipos_2[[4]][[j]]==1, "Sí", "No"), threshold = lista_thresholds[i])
    prec3 <- calc_accuracy(return_tipos_3[[3]][[j]], ifelse(return_tipos_3[[4]][[j]]==1, "Sí", "No"), threshold = lista_thresholds[i])
    media_temp <- round(mean(c(prec1, prec2, prec3)), 2)
    lista_medias_temp <- c(lista_medias_temp, media_temp)
  }
  tabla_precisiones_tipo[paste0("Threshold_", lista_thresholds[i])] <- lista_medias_temp
}
tabla_precisiones_tipo



roc(return_tipos_2[[4]][[4]], return_tipos_3[[3]][[4]], plot = TRUE, legacy.axes = TRUE,
    percent = TRUE, xlab = "Porcentaje Falsos positivos",
    ylab = "Porcentaje verdaderos postivios", col = "#377eb8", lwd = 2,
    print.auc = TRUE, print.auc.x =45, partial.auc = c(100, 90), # en terminos de especificidad
    auc.polygon = TRUE, auc.polygon.col = "#377eb850")

library(ggplot2)
# install.packages("plotROC")



ggplot(data = data.frame(True=return_tipos_2[[4]][[4]], Pred= return_tipos_3[[3]][[4]]), aes(d = True, m = Pred)) + geom_roc()
ggplot(data = data.frame(True = TrueValues_pac, Pred = predictions_pac), aes(d = True, m = Pred)) +
  geom_roc() + labs(x = "1-Especificidad", y = "Sensibilidad") + ggtitle("Modelo sin HER2") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")




