install.packages("xgboost", repos=c("http://dmlc.ml/drat/", getOption("repos")), type="source")
vec.pac= c("SuperLearner", "gbm", "glmnet","ranger","caret")

lapply(vec.pac, require, character.only = TRUE) 





#Learner Library:
learners <- c( "SL.glmnet","SL.xgboost", "SL.ranger","SL.nnet")

#CV Control for the SuperLearner
control <- SuperLearner.CV.control(V=5)



DR_learner <- function(data,covariates,learners){

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
  
  
  ## DR-learner 
  # Train a classification model to get the propensity scores
  p_mod <- SuperLearner(Y = df_aux$d, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                        verbose = FALSE, method = "method.NNLS", family = binomial(),cvControl = control)
  
  p_hat <- p_mod$SL.predict
  p_hat = ifelse(p_hat<0.025, 0.025, ifelse(p_hat>.975,.975, p_hat)) # Overlap bounding
  
  
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
  
  
  
  # Apply the doubly-robust estimator 
  y_mo <- (m1_hat - m0_hat) + ((df_main$d*(df_main$y -m1_hat))/p_hat) - ((1-df_main$d)*(df_main$y - m0_hat)/(1-p_hat))
  
  
  
  ## Collect all pseudo outcomes
  pseudo_all[,1][df_main$ID] <- y_mo
  pseudo_all[,2][df_main$ID] <- p_hat
  
  
  
  ######################
  
}


# DR-learner final estimate  
res_combined_dr <- matrix(NA,nrow(data),5)
# loop
for(l in 1:10){
  
  set <- seq(from=1, to=nrow(data)+1,by=(nrow(data)/10))
  set
  
  if(l<=5){
    dr_mod_cf <- SuperLearner(Y = pseudo_all[(set[l]:set[l+1]-1),1],X=data[(set[l]:set[l+1]-1),covariates], 
                              newX = data[(set[6]:set[11]-1),covariates], SL.library = learners,
                               verbose = FALSE, method = "method.NNLS",cvControl = control)
    
    score_dr_1_cf <- dr_mod_cf$SL.predict
    res_combined_dr[(set[6]:set[11]-1),l] <- score_dr_1_cf
  }
  if(l>5){
    dr_mod_cf <- SuperLearner(Y =pseudo_all[(set[l]:set[l+1]-1),1],X=data[(set[l]:set[l+1]-1),covariates], 
                              newX =data[(set[1]:set[6]-1),covariates], SL.library = learners,
                              verbose = FALSE, method = "method.NNLS",cvControl = control)
    
    score_dr_0_cf <- dr_mod_cf$SL.predict
    res_combined_dr[(set[1]:set[6]-1),(l-5)] <- score_dr_0_cf
  }
  
}
  score_dr_oob <- rowMeans(res_combined_dr)

return(score_dr_oob)
}

