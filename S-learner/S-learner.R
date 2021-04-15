install.packages("xgboost", repos=c("http://dmlc.ml/drat/", getOption("repos")), type="source")
vec.pac= c("SuperLearner", "gbm", "glmnet","ranger","nnet","caret")

lapply(vec.pac, require, character.only = TRUE) 




#Learner Library:
learners <- c( "SL.glmnet","SL.xgboost", "SL.ranger","SL.nnet")

#CV Control for the SuperLearner
control <- SuperLearner.CV.control(V=5)



S_learner <- function(data,covariates,learners){
  
  data$ID <- c(1:nrow(data))
  
  score_S <- matrix(0,nrow(data),1)
  
  set.seed(1234)
  folds <- createFolds(data$d,k=5)
  
  
  for(f in 1:(length(folds))){
    
    if(f == 1){
      data1 <- data[c(folds[[5]],folds[[2]],folds[[3]],folds[[4]]),]
      df_main <- data[folds[[1]],]
    } 
    if(f == 2){
      data1 <- data[c(folds[[1]],folds[[5]],folds[[3]],folds[[4]]),]
      df_main <- data[folds[[2]],]
    } 
    
    if(f == 3){
      data1 <- data[c(folds[[1]],folds[[2]],folds[[5]],folds[[4]]),]
      df_main <- data[folds[[3]],]
    } 
    
    if(f == 4){
      data1 <- data[c(folds[[1]],folds[[2]],folds[[3]],folds[[5]]),]
      df_main <- data[folds[[4]],]
    } 
    
    if(f == 5){
      data1 <- data[c(folds[[1]],folds[[2]],folds[[3]],folds[[4]]),]
      df_main <- data[folds[[5]],]
    } 
    
    df_aux <- data1
    
    X_train <- (df_aux[,c(covariates,"d")])
    
    
    # Train a regression model using the covariates and the treatment variable
    m_mod <- SuperLearner(Y = df_aux$y, X = X_train, SL.library = learners,
                           verbose = FALSE, method = "method.NNLS",cvControl = control)
    
    # Set treatment variable to 0
    X_test_0 <- (df_main[,c(covariates,"d")])
    X_test_0$d <- 0
    
    # Set treatment variable to 1
    X_test_1 <- (df_main[,c(covariates,"d")])
    X_test_1$d <- 1
    
    # Estimate the CATE as the difference between the model with different treatment status
    score_S[,1][df_main$ID] = predict(m_mod,X_test_1)$pred - predict(m_mod,X_test_0)$pred
    
  }
  
  return(score_S)
}
