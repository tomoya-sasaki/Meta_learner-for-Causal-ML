[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **Bootstrap_CI** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of Quantlet: Bootstrap_CI

Published in: 'CATE meets MLE - A tutorial'

Description: Estimates confidence intervals based on the bootstrapping and a gaussian assumption



Keywords: 'CATE, ML, bootstrap,confidence-intervals, causal-inference, treatment'

Author: 'Daniel Jacob'

See also: 'https://www.pnas.org/content/116/10/4156'

Submitted:  '10.02.2021'

```

### R Code
```r




# Here we use bootstrapped data (pred_B) to estimate the standard deviation and a non bootstrapped but cross-fitted data (score) as the point estimate.
                 
bootCI <- function(pred_B,score){
  
  
  # the the 5% and 95% CI from the bootstrapped procedure
  CI_b <- data.frame(
    X5. =  apply(pred_B, 1, function(x)
      quantile(x, c(.025))),
    X95. = apply(pred_B, 1, function(x)
      quantile(x, c(.975))),
    sd = apply(pred_B, 1, function(x) sd(x))
  )
  
  return(data.frame(
    pred = score,
    X5. =  score - 1.96 * CI_b$sd,
    X95. = score + 1.96 * CI_b$sd
  ))
}

                 
 
# An alternative is to use a bootstrapped version for the point of interest.               
# This function uses the mean of the bootstrapped estimates as the point of interest ( a smoothed bootstrap version).
bootCI <- function(pred_B){
  

  # the the 5% and 95% CI from the bootstrapped procedure
  CI_b <- data.frame(
    X5. =  apply(pred_B, 1, function(x)
      quantile(x, c(.025))),
    X95. = apply(pred_B, 1, function(x)
      quantile(x, c(.975))),
    sd = apply(pred_B, 1, function(x) sd(x))
  )
  
  return(data.frame(
    pred = apply(pred_B,1,mean),
    X5. =  apply(pred_B,1,mean) - 1.96 * CI_b$sd,
    X95. = apply(pred_B,1,mean) + 1.96 * CI_b$sd
  ))
}


```

automatically created on 2021-04-21
