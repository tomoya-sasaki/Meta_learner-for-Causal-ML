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
1. Split data into 5 subsamples
2. Use all but sample $k$ to train the nuisance functions and sample $k$ predict the nuisance paramters.
3. Create the pseudo-outcome for sample $k$.


The empirical examples and simulation study contain replecation code from the paper
[Jacob (2020) - CATE Meets ML](http://ssrn.com/abstract=3816558)
