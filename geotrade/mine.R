###############################################################
# preamble
###############################################################

vec.pac <- c("foreign", "quantreg", "gbm", "glmnet",
            "MASS", "rpart", "nnet", "matrixStats",
            "xtable", "readstata13","grf","remotes",
            "caret",  "multcomp","cowplot","SuperLearner",
            "ranger","reshape2","gridExtra","bartCause",
            "xgboost","bartMachine","nnet")
install.packages(vec.pac)
remotes::install_github("vdorie/bartCause")

lapply(vec.pac, require, character.only = TRUE)

# load data
data <- readRDS("df_merged.rds")

data$ID <- c(1:nrow(data))

covariates <- colnames(data)[! colnames(data) %in% c("y", "d", "ID")]


# erase some missing values
# NO NEED b/c already removed
# d2 <- data[complete.cases(data),]



B <- 4 # Number of Bootstrap repetitions

#Learner Library:

#SL.ranger_td = create.Learner("SL.ranger", params = list(num.trees = 1000, min.node.size = 10))
learners <- c("SL.ranger")

#CV Control for the SuperLearner
control <- SuperLearner.CV.control(V=5)





# Create a matrix to store the CATE results from each method
results_cate_DR <- matrix(0,nrow(data),B)
results_cate_R <- matrix(0,nrow(data),B)
results_cate_T <- matrix(0,nrow(data),B)
results_cate_X <- matrix(0,nrow(data),B)
cate_CBART <- as.data.frame(matrix(0,nrow(data),3))
cate_CForest <- as.data.frame(matrix(0,nrow(data),3))


createbootstrappedData <- function(df_boot) {

  smpl_0 <- sample((1:nrow(df_boot))[df_boot$d == 0],
                   replace = TRUE,
                   size = sum(1 - df_boot$d))
  smpl_1 <- sample((1:nrow(df_boot))[df_boot$d == 1],
                   replace = TRUE,
                   size = sum(df_boot$d))
  smpl <- sample(c(smpl_0, smpl_1))

  return(df_boot[smpl,])
}

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



  for(b in 1:B){

    set.seed(1011+b)

    ### Apply the meta-learners with Bootstrapping

    df_aux <- createbootstrappedData(data1)

    ## DR-learner
    # Train a classification model to get the propensity scores
    p_mod <- SuperLearner(Y = df_aux$d, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                          verbose = FALSE, method = "method.NNLS", family = binomial(),cvControl = control)

    p_hat <- p_mod$SL.predict
    p_hat = ifelse(p_hat<0.025, 0.025, ifelse(p_hat>.975,.975, p_hat)) # Overlap bounding

    # Prop-Score for df_main
    #p_hat_main <- predict(p_mod,df_main[,covariates])$pred

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





    a  <- tryCatch({
      dr_mod <- SuperLearner(Y = y_mo, X = df_main[,covariates], newX = df_main[,covariates], SL.library = learners,
                             verbose = FALSE, method = "method.NNLS",cvControl = control)

      score_dr <- dr_mod$SL.predict
      a <- score_dr


    },error=function(e){

      mean_score <- mean(y_mo)
      score_dr <- rep.int(mean_score, times = nrow(df_main))
      a <- score_dr
      return(a)
    })

    score_dr <- a

    ######################

    results_cate_DR[df_main$ID,b] <- score_dr
    results_cate_T[df_main$ID,b] <- (predict(m1_mod,df_main[,covariates])$pred - predict(m0_mod,df_main[,covariates])$pred)

    ### R-learner

    # Train a regression model
    m_mod <- SuperLearner(Y = df_aux$y, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                          verbose = FALSE, method = "method.NNLS",cvControl = control)

    m_hat <- m_mod$SL.predict

    # Apply the R-learner (residual-on-residual approach)
    y_tilde = df_main$y - m_hat
    w_tilde = df_main$d - p_hat
    pseudo_outcome = y_tilde/w_tilde

    weights = w_tilde^2


    a  <- tryCatch({
      R_mod <- SuperLearner(Y = pseudo_outcome, X = df_main[,covariates], newX = df_main[,covariates], SL.library = learners,
                            verbose = FALSE, method = "method.NNLS",obsWeights = weights[,1],cvControl = control)
      score_R <- R_mod$SL.predict
      a <- score_R


    },error=function(e){

      mean_score <- weighted.mean(pseudo_outcome, w = weights)
      score_R <- rep.int(mean_score, times = nrow(df_main))
      a <- score_R
      return(a)
    })

    score_R <- a



    ###########

    results_cate_R[df_main$ID,b] <- score_R

    ###  X-learner

    tau1 <- df_main[which(df_main$d==1),"y"] - m0_hat[which(df_main$d==1),]

    tau0 <- m1_hat[which(df_main$d==0),]  -  df_main[which(df_main$d==0),"y"]


    a1  <- tryCatch({
      tau1_mod <- SuperLearner(Y = tau1, X = df_main[which(df_main$d==1),covariates], newX = df_main[,covariates], SL.library = learners,
                               verbose = FALSE, method = "method.NNLS",cvControl = control)

      score_tau1 <- tau1_mod$SL.predict
      a1 <- score_tau1


    },error=function(e){

      mean_score <- mean(tau1)
      score_tau1 <- rep.int(mean_score, times = nrow(df_main))
      a1 <- score_tau1
      return(a1)
    })

    score_tau1 <- a1

    a0  <- tryCatch({
      tau0_mod <- SuperLearner(Y =tau0, X = df_main[which(df_main$d==0),covariates], newX = df_main[,covariates], SL.library = learners,
                               verbose = FALSE, method = "method.NNLS",cvControl = control)

      score_tau0 <- tau0_mod$SL.predict
      a0 <- score_tau0


    },error=function(e){

      mean_score <- mean(tau0)
      score_tau0 <- rep.int(mean_score, times = nrow(df_main))
      a0 <- score_tau0
      return(a0)
    })

    score_tau0 <- a0



    score_X <- p_hat*score_tau0 + (1-p_hat)*score_tau1

    results_cate_X[df_main$ID,b] <- score_X




    cat("This is Iteration: ", b, "out of", B,"\n")
  }





  ### Causal BART
  fit_bart <- bartc(data1$y,data1$d,data1[,covariates],keepTrees = TRUE,ntree=2000)
  score_bart <- predict(fit_bart,newdata = df_main,type="icate")
  score_bart_m <- apply(score_bart,2,mean)
  ite.sd <- apply(score_bart, 2, sd)
  ite.lb <- score_bart_m - 1.96 * ite.sd
  ite.ub <- score_bart_m + 1.96 * ite.sd


  cate_CBART[df_main$ID,] <- as.data.frame(cbind(score_bart_m,ite.lb,ite.ub))
  colnames(cate_CBART) <- c("pred","X5.","X95.")

  ### Causal Forest

  forest.D <- regression_forest(data1[,covariates], data1$d, tune.parameters = "all")
  W.hat <- predict(forest.D)$predictions

  forest.Y <- regression_forest(data1[,covariates], data1$y, tune.parameters = "all")
  Y.hat <- predict(forest.Y)$predictions

  tau.forest <- causal_forest(data1[,covariates], data1$y, data1$d,tune.parameters = "all",num.trees = 8000,
                              min.node.size=10,W.hat = W.hat, Y.hat = Y.hat)
  tau.hat <- predict(tau.forest, df_main[,covariates],estimate.variance = TRUE)
  summary(tau.hat$predictions)
  sigma.hat <- sqrt(tau.hat$variance.estimates)

  pred = tau.hat$predictions
  X5. =  tau.hat$predictions - 1.96 * sigma.hat
  X95. = tau.hat$predictions + 1.96 * sigma.hat

  cate_CForest[df_main$ID,] <- as.data.frame(cbind(pred,X5.,X95.))
  colnames(cate_CForest) <- c("pred","X5.","X95.")

  print(f)

}


# Estimate CATE for each observation
learners <- c("SL.ranger","SL.xgboost","SL.nnet")


score_T <- matrix(0,nrow(data),1)

pseudo_all <- matrix(NA,nrow(data),5)
ID_pseudo <- 1:nrow(data)
pseudo_all <- cbind(pseudo_all,ID_pseudo)


res_learner_p <- matrix(0,length(learners),5*2)
res_learner_m1 <- matrix(0,length(learners),5*2)
res_learner_m0 <- matrix(0,length(learners),5*2)
res_learner_dr <- matrix(0,length(learners),5*2)
res_learner_m <- matrix(0,length(learners),5*2)
res_learner_R <- matrix(0,length(learners),5*2)




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
  pseudo_all[,5][df_main$ID] <- p_hat



  score_T[,1][df_main$ID] = predict(m1_mod,df_main[,covariates])$pred - predict(m0_mod,df_main[,covariates])$pred



  ######################



  ### R-learner

  # Train a regression model
  m_mod <- SuperLearner(Y = df_aux$y, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                        verbose = FALSE, method = "method.NNLS",cvControl = control)

  m_hat <- m_mod$SL.predict

  # Apply the R-learner (residual-on-residual approach)
  y_tilde = df_main$y - m_hat
  w_tilde = df_main$d - p_hat
  pseudo_outcome = y_tilde/w_tilde

  weights = w_tilde^2
  ## Collect all pseudo outcomes
  pseudo_all[,2][df_main$ID] <- pseudo_outcome
  pseudo_all[,3][df_main$ID] <- weights




  ###########



  ###  X-learner

  tau1 <- df_main[which(df_main$d==1),"y"] - m0_hat[which(df_main$d==1),]

  tau0 <- m1_hat[which(df_main$d==0),]  -  df_main[which(df_main$d==0),"y"]

  ## Collect all pseudo outcomes
  pseudo_all[,4][ (df_main$ID[df_main$d==1])] <- tau1
  pseudo_all[,4][ (df_main$ID[df_main$d==0])] <- tau0

  print(f)


  a <- c(seq(1,2*length(folds),by=2))
  a <- a[f]
  res_learner_p[,a] <- p_mod$cvRisk
  res_learner_p[,a+1] <- p_mod$coef

  res_learner_m1[,a] <- m1_mod$cvRisk
  res_learner_m1[,a+1] <- m1_mod$coef

  res_learner_m0[,a] <- m0_mod$cvRisk
  res_learner_m0[,a+1] <- m0_mod$coef


  res_learner_m[,a] <- m_mod$cvRisk
  res_learner_m[,a+1] <- m_mod$coef
}
