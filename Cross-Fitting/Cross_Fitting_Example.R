if(!require("SuperLearner")) install.packages("SuperLearner"); library("SuperLearner")
install.packages("xgboost", repos=c("http://dmlc.ml/drat/", getOption("repos")), type="source")
vec.pac= c("mvtnorm","clusterGeneration","SuperLearner", "gbm", "glmnet",
           "MASS", "rpart", "doParallel", "sandwich", "randomForest",
           "nnet", "matrixStats", "xtable", "car", "lfe","dplyr",
           "caret", "multcomp", "dplyr","ranger","cowplot", "ggplot2" ,"reshape","e1071","foreach","doParallel")

lapply(vec.pac, require, character.only = TRUE) 
cbp <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



#Learner Library
learners <- c( "SL.glmnet", "SL.ranger","SL.lm","SL.mean") # SVM excluded since it does not offer "obs.weights" variable 
learners <- c( "SL.glmnet")

#CV Control
control <- SuperLearner.CV.control(V=5)







N <- 2000
k <- 10

# create matrix of DGP settings
settings <- as.data.frame(matrix(c(N,N,N,N,N,N,N,N,N,N,N,N,
                                   k,k,k,k,k,k,k,k,k,k,k,k), ncol=2))
settings$V3 <- c(T,"imbalanced","linear","interaction","non-lin","linear",T,"imbalanced","linear","interaction","non-lin","linear") # treatment assigment T = random; lin = linear; interaction = non linear in coefficients, imbalanced = P(D=1) = 0.8; non_lin = sin and cos functions
settings$V4 <- c("con_lin","con_lin","con_non","binary","con_non","zero","con_non","binary","con_non","zero","con_lin","con_lin") # Settings for theta


settings <- settings[12,]





all_results_DR <- list()

M <- 50
d <- 1
  

    
    
    
    
    
    
    
  
    
    # CATE 
    cate_50_50 <- matrix(NA,500,M)
  
    
    cate_50_50_cross <- matrix(NA,500,M)

    
    mse_all <- matrix(0,2,M,dimnames = list(c( "CATE_50_50", "CATE_50_50_cross"),
                                            NULL))
    
    for(b in 1:M) {
 
    data <- datagen(y="con", N=settings[d,1],k=settings[d,2],random_d=settings[d,3],theta=settings[d,4],var=1)
    test_index <- createDataPartition(data$d,p=0.25, list=FALSE)
    
    test_data <- data[test_index,]
    data <- data[-test_index,]
    
    
    ID <- c(1:nrow(test_data))
    test_data$ID <- ID
    data$ID <- c(1:nrow(data))
    
    

    
    
    
    
    

    
    res_folds_cate_50_50 <- matrix(0,nrow(test_data),2)

    
    
    
    
    k <- ncol(data)-4
    covariates <- c(paste0("V", 1:k)) 
    
    
    
    
    
    ##### SuperLearner_DR #####
  
      
      
      ##### 50:50 split ######
      trainIndex <- createDataPartition(data$d, p = .5, list = FALSE) 
      
      for(t in 1:2){
        
        if(t == 1){
          df_aux <- data[trainIndex,] 
          df_main  <- data[-trainIndex,]
          
        } 
        if(t == 2){
          df_main <- data[trainIndex,] 
          df_aux  <- data[-trainIndex,]
        }
        
        
        ## Normal two split
        
        
        p_mod <- SuperLearner(Y = df_aux$d, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                              verbose = FALSE, method = "method.NNLS", family = binomial())
        
        p_hat <- p_mod$SL.predict
        p_hat = ifelse(p_hat<0.025, 0.025, ifelse(p_hat>.975,.975, p_hat)) # Overlap bounding
        
        
        m_mod <- SuperLearner(Y = df_aux$y, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                              verbose = FALSE, method = "method.NNLS")
        
        m_hat <- m_mod$SL.predict
        
        y_tilde = df_main$y - m_hat
        w_tilde = df_main$d - p_hat
        pseudo_outcome = y_tilde/w_tilde
        
        weights = w_tilde^2
        
        
        a  <- tryCatch({
          R_mod <- SuperLearner(Y = pseudo_outcome, X = df_main[,covariates], newX = test_data[,covariates], SL.library = learners,
                                verbose = FALSE, method = "method.NNLS",obsWeights = weights[,1])
          score_R <- R_mod$SL.predict
          a <- score_R
          
          
        },error=function(e){
          
          mean_score <- weighted.mean(pseudo_outcome, w = weights)
          score_R <- rep.int(mean_score, times = nrow(test_data))
          a <- score_R
          return(a)
        })
        
        score_R <- a
        
        
        
        
        res_folds_cate_50_50[,t] <- score_R
        
      }
      
      cate_50_50[,b] <- score_R
      cate_50_50_cross[,b]<- rowMeans(res_folds_cate_50_50)
      
     
      ##########################################################
      
      
      mse_all[1,b] <- mean((as.numeric(cate_50_50[,b]) - as.numeric(test_data$theta))^2)
      mse_all[2,b] <- mean((as.numeric(cate_50_50_cross[,b]) - as.numeric(test_data$theta))^2)
      
   
      
      print(paste0("................................... ","The current iteration is: ", b, " out of " ,M))
      
    
    
    
    
    
    
    
    
    
    }
    
 
    

  
  
  
  mm <- melt(mse_all)
  colnames(mm) <- c("Estimator","repetitions","value")
  mm$Estimator <- rep(c("single","cross-fit"),length.out=nrow(mm))
  mm$Estimator <- factor(mm$Estimator,levels=c("single","cross-fit"))
  
  
  
  ggplot(mm, aes(x=repetitions, y=value,group=Estimator))+
    geom_point(aes(color=Estimator, shape=Estimator),size=3)+
    #geom_line(data=mm[mm$Model=="DR_3folds_cross (median)"| mm$Model=="DR_2folds (median)"| mm$Model=="T_cross (median)" | mm$Model=="DR_3folds (median)", ],
    #         aes(color=factor(Model))) +
    theme_cowplot() +
    #facet_wrap( ~ DGP, scales="free", nrow=4) + # Facet wrap with common scales 
    labs(y = "MSE ", x = "Replication") +
    scale_color_manual(values = cbp) +
    theme(legend.position="bottom", legend.justification = 'center') +
    scale_shape_manual(values = c(16,17,15,3,7,8,42,43))
  
  
  
  
  
  
  

