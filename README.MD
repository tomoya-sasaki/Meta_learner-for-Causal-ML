# CATE meets ML - Collection of machine learning (ML) methods to estimate the conditional average treatment effect (CATE)

The repository contains code for meta-learners and ML-methods designed for heterogeneout treatment effect estimation 

* Causal BART
* Causal Boosting
* Causal Forest (based on GRF)
* DR-learner
* R-learner
* S-learner
* T-learner
* TO-leanrer
* X-learner

# Estimation procedure

We aim to estimate the CATE on the whole sample. We proceed as follows:

1. Split data into *5* folds.
2. Use all but fold *k* to train the nuisance functions and fold *k* predict the nuisance paramters.
3. Create the pseudo-outcome for fold *k*.
4. Loop through all the folds such that every subsample is used as sample *k*. 
5. Collect all the nuisance parameters and pseudo-outcomes.
6. Create a new dataset that contains the pseudo-outcomes and covariates (for the R-learner also the weights).
7. Split this dataset into two folds (A and M).
8. Use A to regress the covariates on the pseudo-outcome and M to predict the *CATE_M*.
9. Reverse the role of the samples to obtain *CATE_A*.
10. Combine *CATE_A* and *CATE_M* (according to their ID in the dataset) and obtain the final *CATE*.

# Reference

The empirical examples and simulation study contain replecation code from the paper
[Jacob (2020) - CATE Meets ML](http://ssrn.com/abstract=3816558)
