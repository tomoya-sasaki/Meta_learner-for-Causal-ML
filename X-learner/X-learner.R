install.packages("xgboost", repos=c("http://dmlc.ml/drat/", getOption("repos")), type="source")
vec.pac= c("SuperLearner", "gbm", "glmnet","ranger","caret")

lapply(vec.pac, require, character.only = TRUE) 





#Learner Library:
learners <- c( "SL.glmnet","SL.xgboost", "SL.ranger","SL.nnet")

#CV Control for the SuperLearner
control <- SuperLearner.CV.control(V=5)


X_learner <- function(data,covariates,learners){

  data$ID <- c(1:nrow(data))

pseudo_all <- matrix(NA,nrow(data),2)
ID_pseudo <- 1:nrow(data)
pseudo_all <- cbind(pseudo_all,ID_pseudo)




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
  
  
  
  # Train a classification model to get the propensity scores
  p_mod <- SuperLearner(Y = df_aux$d, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                        verbose = FALSE, method = "method.NNLS", family = binomial(),cvControl = control)
  
  p_hat <- p_mod$SL.predict
  p_hat = ifelse(p_hat<0.025, 0.025, ifelse(p_hat>.975,.975, p_hat)) # Overlap bounding
  
  
  ## Collect propensity-score

  pseudo_all[,2][df_main$ID] <- p_hat
  
  # Split the training data into treatment and control observations
  aux_1 <- df_aux[which(df_aux$d==1),]
  aux_0 <- df_aux[which(df_aux$d==0),]
  
  # Train a regression model for the treatment observations
  m1_mod <- SuperLearner(Y = aux_1$y, X = aux_1[,covariates], newX = df_main[,covariates], SL.library = learners,
                         verbose = FALSE, method = "method.NNLS",cvControl = control)
  
  m1_hat <- m1_mod$SL.predict
  
  # Train a regression model for the control observations
  m0_mod <- SuperLearner(Y = aux_0$y, X = aux_0[,covariates], newX = df_main[,covariates], SL.library = learners,
                         verbose = FALSE, method = "method.NNLS",cvControl = control)
  
  m0_hat <- m0_mod$SL.predict
  
  
  
  
  
  ######################
  
  
  
  
  
  ###  X-learner
  
  tau1 <- df_main[which(df_main$d==1),"y"] - m0_hat[which(df_main$d==1),]
  
  tau0 <- m1_hat[which(df_main$d==0),]  -  df_main[which(df_main$d==0),"y"]
  
  ## Collect all pseudo outcomes
  pseudo_all[,1][ (df_main$ID[df_main$d==1])] <- tau1
  pseudo_all[,1][ (df_main$ID[df_main$d==0])] <- tau0
  
 

  
}




# X-learner final estimate
a1  <- tryCatch({
  tau1_mod <- SuperLearner(Y = pseudo_all[,1][data$d==1], X = data[which(data$d==1),covariates], newX = data[,covariates], SL.library = learners,
                           verbose = FALSE, method = "method.NNLS",cvControl = control)
  
  score_tau1 <- tau1_mod$SL.predict
  a1 <- score_tau1
  
  
},error=function(e){
  
  mean_score <- mean(pseudo_all[,1])
  score_tau1 <- rep.int(mean_score, times = nrow(data))
  a1 <- score_tau1
  return(a1)
})

score_tau1 <- a1

a0  <- tryCatch({
  tau0_mod <- SuperLearner(Y =pseudo_all[,1][data$d==0], X = data[which(data$d==0),covariates], newX = data[,covariates], SL.library = learners,
                           verbose = FALSE, method = "method.NNLS",cvControl = control)
  
  score_tau0 <- tau0_mod$SL.predict
  a0 <- score_tau0
  
  
},error=function(e){
  
  mean_score <- mean(pseudo_all[,1])
  score_tau0 <- rep.int(mean_score, times = nrow(data))
  a0 <- score_tau0
  return(a0)
})

score_tau0 <- a0



score_X <- pseudo_all[,2]*score_tau0 + (1-pseudo_all[,2])*score_tau1


  
  return(score_X)
  }
