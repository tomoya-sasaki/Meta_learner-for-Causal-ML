[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **T-learner** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml


Name of Quantlet: T-learner

Published in: 'CATE meets MLE - A tutorial'

Description: Two model method to estimate the conditional average treatment effect (CATE) via a variety of machine learning (ML) methods. 

The data structure has to be the following: 

- y = outcome variable
- d = treatment variable
- covariates = the covariates to map on d or y
- learners = the machine learning methods to use in the [SuperLearner package](https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html).
- df_aux and df_main = the train and test-set. df_aux is used for training the nuisance functions while df_main is used to estimate the CATE. 


Keywords: 'CATE, ML, two model, causal-inference, treatment'

Author: 'Daniel Jacob'

See also: ''

Submitted:  '10.02.2021'

```

### R Code
```r

install.packages("xgboost", repos=c("http://dmlc.ml/drat/", getOption("repos")), type="source")
vec.pac= c("SuperLearner", "gbm", "glmnet","ranger")

lapply(vec.pac, require, character.only = TRUE) 




#Learner Library:
learners <- c( "SL.glmnet","SL.xgboost", "SL.ranger","SL.lm","SL.mean")

#CV Control for the SuperLearner
control <- SuperLearner.CV.control(V=5)



T_learner <- function(df_aux,df_main,covariates,learners){
  
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
  
  # Estimate the CATE as the difference between the two models
  score_t <- m1_hat - m0_hat
  
  return(score_t)
}

```

automatically created on 2021-02-21