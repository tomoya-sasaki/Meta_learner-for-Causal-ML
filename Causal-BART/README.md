[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **Causal_BART** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of Quantlet: Causal_BART

Published in: 'CATE meets MLE - A tutorial'

Description: Causal Bayesian additive regression tree method to estimate the conditional average treatment effect (CATE) via a variety of machine learning (ML) methods. 

The data structure has to be the following: 

- y = outcome variable
- d = treatment variable
- covariates = the covariates to map on d or y


Keywords: 'CATE, ML, causal BART, causal-inference, treatment'

Author: 'Daniel Jacob'

See also: 'https://rdrr.io/github/vdorie/bartCause/man/bartc.html'

Submitted:  '10.02.2021'

```

### R Code
```r

install.packages("remotes")
remotes::install_github("vdorie/bartCause")

# Fits a collection of treatment and response models using the Bayesian Additive
# Regression Trees (BART) algorithm, producing estimates of treatment effects.


Causal_BART <- function(df_aux,df_main,covariates){

fit_bart <- bartc(df_aux$y,df_aux$d,df_aux[,covariates],keepTrees = TRUE)
score_bart <- predict(fit_bart,newdata = df_main,type="icate")
score_bart_m <- apply(score_bart,2,mean)
ite.sd <- apply(score_bart, 2, sd)
ite.lb <- score_bart_m - 2 * ite.sd
ite.ub <- score_bart_m + 2 * ite.sd

cate_CBART <- as.data.frame(cbind(score_bart_m,ite.lb,ite.ub))
colnames(cate_CBART) <- c("pred","X5.","X95.")

return(cate_CBART)

}


```

automatically created on 2021-02-21