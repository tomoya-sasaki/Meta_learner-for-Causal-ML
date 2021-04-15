### This generates the Simulation Data as a dataframe
# by Daniel Jacob (daniel.jacob@hu-berlin.de) 

# Arguments to specify are: 

# N = Number of observations (real number)
# k = Number of covariates (real number)
# random_d = treatment assignment: (Either T for random assignment or other values for confounding on X)
# theta = treatment effect: (Either real number for only one theta, or low dimensional {-0.5,0,+0.5} or "con" for continuous values)
# var = Size of the variance (Noise-level)

#Required Packages
if(!require("clusterGeneration")) install.packages("clusterGeneration"); library("clusterGeneration")
if(!require("mvtnorm")) install.packages("mvtnorm"); library("mvtnorm")


datagen <- function(N=N,y,k,random_d,theta,var,iter) {
  
  # fixed components
  b = 1 / (1:k)
  
  # = Generate covariance matrix of z = #
  set.seed(123*iter)
  sigma <- genPositiveDefMat(k, "unifcorrmat")$Sigma
  set.seed(123*iter)
  sigma <- cov2cor(sigma)
  
  set.seed(123*iter)
  z_fix <- rmvnorm(N, sigma = sigma) # = Generate z = #
  
  
  ### Options for theta
  
  # create subgroup effect
  z_fix[,5] <- rep(c(-0.2,0,0.2,0.6),length.out=nrow(z_fix))
  
  set.seed(123*iter)
  theta_con_lin <- as.vector(0.6*z_fix[,1] + 0.6*z_fix[,2] + 0.6*z_fix[,3] + 0.6*z_fix[,4] + z_fix[,5] + rnorm(N,0,0.5))
  #set.seed(123*iter)
  theta_con_non <- as.vector(sin(z_fix[,1:3] %*% b[1:3]) + 1.5*cos(z_fix[,4])) + z_fix[,5]
  #set.seed(123*iter)
  theta_lowdim <- as.vector(ifelse(z_fix[,1]>0,0.5,-0.5)) + ifelse(z_fix[,5]>0,0.5,-0.5)
  
  theta_interaction <- as.vector(z_fix[,1]*z_fix[,2] + z_fix[,2] + z_fix[,3]*z_fix[,4] + z_fix[,5])
  
  
  
  if(theta=="con_lin")
  {theta <- theta_con_lin}
  else if(theta=="con_non")
  {theta <- theta_con_non}
  else if(theta=="low_dim")
  {theta <- theta_lowdim}
  else if(theta=="interaction")
  {theta <- theta_interaction}
  else 
  {theta <- 0}
  
  
  z <- z_fix
  
  ### Options for D (m_(X))
  if (random_d == T) {
    d <- rep(c(0, 1), length.out = N)
  } else 
    if(random_d =="imbalanced"){
      set.seed(123+iter)
      d <-  as.numeric(rbinom(N,prob=0.2,size=1))
      
    }
  else
    
    
    if(random_d == "linear"){
      d_prop <- pnorm( z[,1] + z[,2] + z[,3] - z[,4]) # D is dependent on Za
      d <- as.numeric(rbinom(N, prob = d_prop, size = 1))
    }
  else 
    if(random_d == "interaction"){
      d_prop <- pnorm( z[,1]*z[,2] + z[,3]*z[,4]) # D is dependent on Za
      d <- as.numeric(rbinom(N, prob = d_prop, size = 1))
    }
  else{
    d_prop <- pnorm( sin(z[,1]) + sin(z[,2]) + cos(z[,3]*z[,4])) # D is dependent on Za
    d <- as.numeric(rbinom(N, prob = d_prop, size = 1))
    
  }
  
  
  
  
  
  
  g <- as.vector(z[,1]*z[,2] + z[,3]*z[,4] + z[,5])
  
  
  if(y=="binary") {
    y1 <- theta * d + g 
    y1.1 <- rbinom(N,prob=pnorm(scale(y1)),size=1)
    #y1.1 <- (y1 - min(y1)) * (1) / (max(y1) - min(y1)) + 0
    #y <-  rbinom(N,prob=y1.1,size=1)
    y <- y1.1
  } else {y <- theta * d + g + rnorm(N,0,var)}
  
  data <- as.data.frame(y)
  data <- cbind(data, theta, d, z)
  colnames(data) <- c("y", "theta", "d", c(paste0("V", 1:k)))
  
  return(data)
}


### Example
df_aux <- datagen(y="con",N = 1000, k = 20, random_d = "linear", theta = "con_non", var = 1,iter=2)
df_main <- datagen(y="con",N = 1000, k = 20, random_d = "linear", theta = "con_non", var = 1,iter=2)

all.equal(df_aux,df_main)
