[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **Causal_Boosting** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of Quantlet: Causal_Boosting

Published in: 'CATE meets MLE - A tutorial'

Description: Causal Boosting (regression tree) method to estimate the conditional average treatment effect (CATE) via a variety of machine learning (ML) methods. 

The data structure has to be the following: 

- y = outcome variable
- d = treatment variable
- covariates = the covariates to map on d or y


Keywords: 'CATE, ML, causal Boosting, causal-inference, treatment'

Author: 'Daniel Jacob'

See also: 'https://github.com/saberpowers/causalLearning'

Submitted:  '10.02.2021'

```

### R Code
```r

remotes::install_github("saberpowers/causalLearning")
library(causalLearning)

Causal_Boost <- function(df_aux,df_main){
fit_cb = causalBoosting(y=df_aux$y,tx=df_aux$d,x=df_aux[,covariates], num.trees = 500)

score_CBoost = predict(fit_cb, newx =df_main[,covariates] , num.trees = 500)

return(score_CBoost)
}

```

automatically created on 2021-02-21