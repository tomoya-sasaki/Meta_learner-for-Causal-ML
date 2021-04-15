install.packages("xgboost", repos=c("http://dmlc.ml/drat/", getOption("repos")), type="source")
vec.pac= c("SuperLearner", "gbm", "glmnet","ranger","caret")

lapply(vec.pac, require, character.only = TRUE) 




#Learner Library:
learners <- c( "SL.glmnet","SL.xgboost", "SL.ranger","SL.nnet")

#CV Control for the SuperLearner
control <- SuperLearner.CV.control(V=5)


IPW_learner <- function(data,covariates,learners){

	data$ID <- c(1:nrow(data))
  
 ##### # 5-fold sample splitting
# Sample splitting
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
  
  
  ## IPW-learner 
  # Train a classification model to get the propensity scores
  p_mod <- SuperLearner(Y = df_aux$d, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                        verbose = FALSE, method = "method.NNLS", family = binomial(),cvControl = control)
  
  p_hat <- p_mod$SL.predict
  p_hat = ifelse(p_hat<0.025, 0.025, ifelse(p_hat>.975,.975, p_hat)) # Overlap bounding
  
  # Transformed outcome
  y_to <- df_main$y * df_main$d/p_hat - df_main$y * (1-df_main$d)/(1-p_hat)

  ## Collect all pseudo outcomes
  pseudo_all[,1][df_main$ID] <- y_to

}
  
  
nhalf = floor(0.5*nrow(data))
n = nrow(data)


to_mod_oob <- SuperLearner(Y = pseudo_all[1:nhalf,1], X = data[1:nhalf,covariates], newX = data[(nhalf+1):n,covariates], SL.library = learners,
                           verbose = FALSE, method = "method.NNLS",cvControl = control)

score_to_1 <- to_mod_oob$SL.predict

to_mod_oob <- SuperLearner(Y = pseudo_all[(nhalf+1):n,1], X = data[(nhalf+1):n,covariates], newX = data[1:nhalf,covariates], SL.library = learners,
                           verbose = FALSE, method = "method.NNLS",cvControl = control)

score_to_0 <- to_mod_oob$SL.predict


score_to_oob= rbind(score_to_0,score_to_1)
 
  
  return(score_to_oob)
  }
