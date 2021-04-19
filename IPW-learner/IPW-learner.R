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
  
  
res_combined_ipw <- matrix(NA,nrow(data),5)
# loop
for(l in 1:10){
  
  set <- seq(from=1, to=nrow(data)+1,by=(nrow(data)/10))
  set
  
  if(l<=5){
    ipw_mod_cf <- SuperLearner(Y = pseudo_all[(set[l]:set[l+1]-1),1],X=data[(set[l]:set[l+1]-1),covariates], 
                              newX = data[(set[6]:set[11]-1),covariates], SL.library = learners,
                               verbose = FALSE, method = "method.NNLS",cvControl = control)
    
    score_ipw_1_cf <- ipw_mod_cf$SL.predict
    res_combined_ipw[(set[6]:set[11]-1),l] <- score_dr_1_cf
  }
  if(l>5){
    ipw_mod_cf <- SuperLearner(Y =pseudo_all[(set[l]:set[l+1]-1),1],X=data[(set[l]:set[l+1]-1),covariates], 
                              newX =data[(set[1]:set[6]-1),covariates], SL.library = learners,
                              verbose = FALSE, method = "method.NNLS",cvControl = control)
    
    score_ipw_0_cf <- ipw_mod_cf$SL.predict
    res_combined_ipw[(set[1]:set[6]-1),(l-5)] <- score_dr_0_cf
  }
  
}
 
score_ipw_oob <- rowMeans(res_combined_ipw)
  
  return(score_ipw_oob)
  }
