# Script para setup

# Source de funciones y otras preparaciones para ejecutar Presentacion.R
source("./filtrar_datos.R")
source("./contar_NAs.R")
source("./src/fun_resumen_dataset.R")
source("./src/fun_dummy.R")
source("./src/fun_normalizacion.R")
source("./src/fun_univar_pvalue.R")


# Script ----
datos <- read_csv("pacientes_cancer3.csv")
filtrados <- filtrar_datos(datos)

datos_HER2 <- filtrados |> select(-CRPlevel) |> remove_na()

# Dummies y escalado
datos_HER2dummy <- nuevo_dataset_dummy(datos_HER2)
datos_HER2escalado <- nuevo_dataset_normalizado(datos_HER2dummy)

# Bivariante HER2
datos_HER2 <- datos_HER2escalado
tabla_pvalores_univar <- table_univar_sig(datos_HER2 |> select(-Remision), datos_HER2$Remision)
tabla_pvalores_univar[tabla_pvalores_univar['P.Valores']<0.05, ]
colnames_signif <- tabla_pvalores_univar['Nombres'][tabla_pvalores_univar['P.Valores']<0.05]
