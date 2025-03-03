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

# Training 
library(caTools)
set.seed(123)
size_training <- 0.8

# splitting data
split <- sample.split(datos_HER2$Remision, SplitRatio = size_training)

# creating training dataset
trainingSet <- subset(datos_HER2, split == TRUE)

# creating test data set
testSet <- subset(datos_HER2, split == FALSE)

# model
cancer_model <- glm(trainingSet$Remision~.,
                    data = trainingSet[sig_variables_HER2],
                    family = "binomial")

# prediction
cancer_model_predict <- predict(cancer_model, testSet, type = "response")

# accuracy
accuracy <- calc_accuracy(cancer_model_predict, ifelse(testSet$Remision==1, "Sí", "No"), threshold = 0.6)
