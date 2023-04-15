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

