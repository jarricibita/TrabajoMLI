#'@title Primer training y testing con split
#'
#'@description col_dataset_dep en formato dataset['NombreCol'], size_training de split (0-1), threshold de precision (0-1)
#'
#'@param dataset
#'@param col_dataset_dep
#'@param size_training
#'@param sig_variable_names
#'@param threshold
#'@param seed
#'@return trainingSet, testSet, accuracy
#'
#'@example
#'train_and_test_first(size_training = 0.8, threshold = 0.5, seed = 123)
#'
#'
#'
train_and_test_first <- function(dataset, col_dataset_dep, size_training, sig_variable_names, threshold, seed = 123){
  set.seed(seed)
  name_dep <- colnames(col_dataset_dep)
  size_training <- size_training
  
  # splitting data
  split <- sample.split(dataset[[paste0(name_dep)]], SplitRatio = size_training)
  
  # creating training dataset
  trainingSet <- subset(dataset, split == TRUE)
  
  # creating test data set
  testSet <- subset(dataset, split == FALSE)
  
  # model
  glm_model <- glm(trainingSet[[paste0(name_dep)]]~.,
                   data = trainingSet[sig_variable_names],
                   family = "binomial")
  table <- data.frame('Variables'=sig_variable_names,
                      'p-valores' = summary(glm_model)$coefficients[-1,4])
  
  # prediction
  glm_model_predict <- predict(glm_model, testSet, type = "response")
  TrueValues <- testSet$Remision
  # accuracy
  accuracy <- calc_accuracy(glm_model_predict, ifelse(testSet[[paste0(name_dep)]]==1, "Sí", "No"), threshold = threshold)
  return(list(trainingSet, testSet, accuracy, glm_model_predict, TrueValues, table))
}
