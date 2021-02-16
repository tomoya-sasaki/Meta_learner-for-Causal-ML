install.packages("xgboost", repos=c("http://dmlc.ml/drat/", getOption("repos")), type="source")
vec.pac= c("SuperLearner", "gbm", "glmnet","ranger")

lapply(vec.pac, require, character.only = TRUE) 


# df_aux is the training set
# df_main is the test set where the CATE is estimated on
# covariates are the pre-treatment covariates 


#Learner Library for the SuperLearner:
learners <- c( "SL.glmnet","SL.xgboost", "SL.ranger","SL.lm","SL.mean")

#CV Control for the SuperLearner
control <- SuperLearner.CV.control(V=5) # The cross-validation parameter 

R <- 10 # The number of repetitions for cross-fitting and median averaging.


R_learner <- function(df_aux,df_main,covariates,learners){
  
  est_cate_R <- matrix(0,nrow(df_main),R)
  
  
  
    for(r in 1:R){
    
    ind <- createDataPartition(df_aux[,"d"], p = .5, list = FALSE)
    
    df_aux1<- as.data.frame(df_aux[ ind,])
    df_aux2 <- as.data.frame(df_aux[-ind,])

m_mod <- SuperLearner(Y = df_aux1$y, X = df_aux1[,covariates], newX = df_aux2[,covariates], SL.library = learners,
                        verbose = FALSE, method = "method.NNLS",cvControl = control)
  
  m_hat <- m_mod$SL.predict
  
  # Apply the R-learner (residual-on-residual approach)
  y_tilde = df_aux2$y - m_hat
  w_tilde = df_aux2$d - p_hat
  pseudo_outcome = y_tilde/w_tilde
  
  weights = w_tilde^2
  
  
  a  <- tryCatch({
    R_mod <- SuperLearner(Y = pseudo_outcome, X = df_aux2[,covariates], newX = df_main[,covariates], SL.library = learners,
                          verbose = FALSE, method = "method.NNLS",obsWeights = weights[,1],cvControl = control)
    score_R <- R_mod$SL.predict
    a <- score_R
    
    
  },error=function(e){
    
    mean_score <- weighted.mean(pseudo_outcome, w = weights)
    score_R <- rep.int(mean_score, times = nrow(df_main))
    a <- score_R
    return(a)
  })
  
  est_cate_R[,r] <- score_dr 
  }
  
  return(apply(est_cate_R,1,median))
}
