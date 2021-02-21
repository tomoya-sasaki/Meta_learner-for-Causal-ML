[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **TO-learner** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml


Name of Quantlet: TO-learner

Published in: 'CATE meets MLE - A tutorial'

Description: Transformed-outcome method to estimate the conditional average treatment effect (CATE) via a variety of machine learning (ML) methods.

The data structure has to be the following: 

- y = outcome variable
- d = treatment variable
- covariates = the covariates to map on d or y
- learners = the machine learning methods to use in the [SuperLearner package](https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html).
- df_aux and df_main = the train and test-set. df_aux is used for training the nuisance functions while df_main is used to estimate the CATE. 


Keywords: 'CATE, ML, transformed outcome, causal-inference, treatment'

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


TOM <- function(df_aux, df_main,covariates,learners){
  
  # Propensity score
  p <- rep((nrow(df_aux[df_aux[,"d"]==1,]) + nrow(df_main[df_main[,"d"]==1,]))/(nrow(df_main) + nrow(df_main)), nrow(df_aux)) 
  
  # Transformed outcome
  Y_TO <- df_aux$y * df_aux$d/p - df_aux$y * (1-df_aux$d)/(1-p)
  
  m_TO <- SuperLearner(Y =Y_TO, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                           verbose = FALSE, method = "method.NNLS",cvControl=control)
  
 
  
  tau_TO <- m_TO$SL.predict
 
  
  return(tau_TO)
  }

```

automatically created on 2021-02-21