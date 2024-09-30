Train machine learning models for the classification of aromatic amine
or organophosphorus structural alerts
================

This is an R Markdown document with relevant chunks from the code used
to train machine learning models for the classification of MS2 spectra
with an organophosphorus or aromatic amine structural alert.

``` r
# Load required packages
library(doParallel) # for parallel computing
library(caret) # for machine learning
library(xgboost) # for machine learning

# define seed
seed <- 234
```

### Load data and prepare parallel computing

``` r
# Load training data
data.train <- read.csv(file = "data\\data_train_aromaticamine_fragments.csv") # example

# Check the number of cores available and reserve one for operations
cores <- parallel::detectCores() - 1

# Set up parallel computing
cl <- makePSOCKcluster(cores)
registerDoParallel(cl)
```

### Model training - support vector machine

``` r
set.seed(seed)

supportvectormachine <- caret::train(alert ~ ., # set variable to predict, use all features for prediction
  data = data.train, # specify training set
  method = "svmRadial", # specify algorithm
  preProcess = NULL, # don't do preprocessing
  metric = "ROC", # optimization metric, in this case AUC-ROC
  maximize = TRUE, # the metric should be maximized
  trControl = trainControl(
    method = "cv", # resampling method
    verboseIter = TRUE, # print training log
    allowParallel = TRUE, # allow parallel processing
    savePredictions = "all", # saves predictions for the optimal tuning parameters
    classProbs = TRUE, # class probabilities are computed in each resample
    summaryFunction = twoClassSummary
  ),
  tuneLength = 5 # number for mtry hyperparameter
)
```

### Model training - extreme gradient boosting

``` r
set.seed(seed)

xgb <- caret::train(alert ~ ., # set variable to predict, use all features for prediction
  data = data.train, # specify training set
  method = "xgbTree", # specify algorithm
  preProcess = NULL, # don't do preprocessing
  metric = "ROC", # optimization metric, in this case AUC-ROC
  maximize = TRUE, # the metric should be maximized
  trControl = trainControl(
    method = "cv", # resampling method
    verboseIter = TRUE, # print training log
    allowParallel = TRUE, # allow parallel processing
    savePredictions = "all", # saves predictions for the optimal tuning parameters
    classProbs = TRUE # class probabilities are computed in each resample
  ),
  tuneLength = 5 # number for mtry hyperparameter
)
```

### Model training - random forest

``` r
set.seed(seed)

randomforest <- caret::train(alert ~ ., # set variable to predict, use all features for prediction
  data = data.train, # specify training set
  method = "rf", # specify algorithm
  preProcess = NULL, # don't do preprocessing
  metric = "ROC", # optimization metric, in this case AUC-ROC
  maximize = TRUE, # the metric should be maximized
  trControl = trainControl(
    method = "cv", # resampling method
    verboseIter = TRUE, # print training log
    allowParallel = TRUE, # allow parallel processing
    savePredictions = "all", # saves predictions for the optimal tuning parameters
    classProbs = TRUE, # class probabilities are computed in each resample
    summaryFunction = twoClassSummary
  ),
  tuneLength = 5 # number for mtry hyperparameter
)
```

### Model training - neural network

``` r
set.seed(seed)

neuralnetwork <- caret::train(alert ~ ., # set variable to predict, use all features for prediction
  data = data.train, # specify training set
  method = "nnet", # specify algorithm
  preProcess = NULL, # don't do preprocessing
  metric = "ROC", # optimization metric, in this case AUC-ROC
  maximize = TRUE, # the metric should be maximized
  trControl = trainControl(
    method = "cv", # resampling method
    verboseIter = TRUE, # print training log
    allowParallel = TRUE, # allow parallel processing
    savePredictions = "all", # saves predictions for the optimal tuning parameters
    classProbs = TRUE, # class probabilities are computed in each resample
    summaryFunction = twoClassSummary
  ),
  MaxNWts = 6000
)
```

### End parallel computing

``` r
stopCluster(cl)
```
