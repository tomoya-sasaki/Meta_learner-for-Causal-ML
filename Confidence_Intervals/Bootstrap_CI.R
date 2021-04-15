
# This function used the mean of the bootstrapped estimates as the point of interest (pred=).
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

# An alternative is to not use a bootstrapped version for the point of interest.
# This means we estimate the CATE based on sample-splitting and use these estimates as the prediction.
                 
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
