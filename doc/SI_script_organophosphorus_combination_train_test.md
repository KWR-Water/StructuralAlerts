Prepare training and test set for organophosphorus classification -
combination
================

This is an R Markdown document with relevant chunks from the code used
to generate a training and test set for the classification of MS2
spectra with an organophosphorus structural alert. Here, the combination
of fragments and neutral losses is made.

``` r
# Load required packages
library(dplyr) # for data handling
library(caret) # for machine learning

# define seed
seed <- 234
```

### Load and prepare data

``` r
# select data type (combination of neutral losses and fragments)
data.type <- "fragments_neutrallosses"

# load both datasets
data.nl <- read.csv(file = "data\\MassBankEU_202206_neutrallosses_pos_hrms.csv")
data.f <- read.csv(file = "data\\MassBankEU_202206_fragments_pos_hrms.csv")

# add rank to data.f to make sure that the order remains the same
data.f$rank <- 1:nrow(data.f)

# combine both datasets and remove duplicate columns (default)
data <- merge(data.nl, data.f)

# order the rows in the way they're originally ordered (like in data.nl, and data.frag)
data <- data[order(data$rank),]
data$rank <- NULL

# select spectra of molecules with the alert of interest, some might contain other alerts as well.
data.subset.alert <- data %>% filter(grepl("TA11484", data$alert) | grepl("TA11485", data$alert)) 
row.names(data.subset.alert) <- NULL
data.subset.alert$alert <- "alert"
data.subset.noalert <- data %>% filter(!grepl("TA11484", data$alert) & !grepl("TA11485", data$alert))
row.names(data.subset.noalert) <- NULL

# get a random sample of spectra without the organophosphorus alert
set.seed(seed)
x <- data.subset.noalert[sample(nrow(data.subset.noalert), size = 7600, replace = FALSE), ]
x$alert <- "none"

# combine both datasets (spectra with alert, and random subset of spectra without the alert)
data.subset <- rbind(data.subset.alert, x)
data.subset <- data.subset %>% relocate(alert, .after = last_col())

# Remove original data set so we have more memory space free
rm(data, data.subset.alert, data.subset.noalert, x)

# load NIST data with organophosphorus spectra
data.nl.nist <- readRDS(file = "data//NIST_organophosphorus_neutrallosses_pos_largerbinsize_hrms.RDS")
data.f.nist <- readRDS(file = "data//NIST_organophosphorus_fragments_pos_largerbinsize_hrms.RDS")

# combine both datasets and remove duplicate columns (default)
data <- merge(data.nl.nist, data.f.nist)

# combine MassBankEU and NIST data
data.subset <- bind_rows(data.subset, data)
data.subset <- data.subset %>% relocate(alert, .after = last_col())

# Remove original data set so we have more memory space free
rm(data, data.f.nist, data.nl.nist)

# replace NA values in the remaining part of the data frame
data.subset[is.na(data.subset)] <- 0
```

### Removal of empty columns, near zero variance predictors and duplicates

``` r
set.seed(seed) 

# Remove empty features (columns), near zero variance predictors and duplicates ----------
x <- colSums(data.subset[, 2:(ncol(data.subset) - 1)], na.rm = TRUE) # count sums per column, do not consider column 'id' and 'alert'
zeros <- x[which(x == 0)] # only select non-occurring fragments / neutral losses
data.subset <- data.subset[, names(zeros) := NULL] # remove these from the data set
rm(x, zeros) # clean environment

# identify near zero variance predictors
zerovariance <- nearZeroVar(
  x = data.subset,
  freqCut = 95/5, # the cutoff for the ratio of the most common value to the second most common value
  saveMetrics = TRUE, # return data frame with specs per predictor
)

# remove near zero variance predictors from data set
zerovariance$freqRatioRounded <- round(zerovariance$freqRatio, digits = 0)
zerovariance.name <- rownames(zerovariance[which(zerovariance$nzv == TRUE & zerovariance$freqRatioRounded == (nrow(data.subset)-1)), ])
data.subset <- data.subset %>%
  select(-all_of(zerovariance.name))

# remove duplicate predictors (i.e. exact the same columns) from data set
tokeep <- which(duplicated(as.list(data.subset)) == FALSE)
data.subset <- data.subset[, ..tokeep]

# clean environment
rm(zerovariance, tokeep, zerovariance.name, alert.count)
```

### Split data into test and training set, with unique chemicals

``` r
# Get metadata; smiles, inchikeys and accession of every spectrum
metadata.mb <- read.csv(file = "data\\MassBankEU_202206_metadata_subset.csv")
metadata.mb.subset <- metadata.mb %>%
  select("accession", "inchikey")
metadata.mb.subset$inchikey.pt1 <- substr(metadata.mb.subset$inchikey, start = 1, stop = 14)

# read metadata from NIST
metadata.nist <- readRDS(file = "data\\2024-04-16_NIST23_metadata_organophosphorus.RDS")
metadata.nist.subset <- as.data.frame(matrix(nrow = nrow(metadata.nist), ncol = 2))
metadata.nist.subset$V1 <- unlist(metadata.nist$filename)
metadata.nist.subset$V2 <- unlist(metadata.nist$inchikey)
colnames(metadata.nist.subset) <- c("accession", "inchikey")
metadata.nist.subset$inchikey.pt1 <- substr(metadata.nist.subset$inchikey, start = 1, stop = 14)

# combine both metadata variables
metadata.combined <- bind_rows(metadata.mb.subset, metadata.nist.subset)

# get unique compounds in the subset, retrieve accession, smiles and inchikey and couple with alert
compounds.unique <- metadata.combined[which(metadata.combined$accession %in% data.subset$id),]
row.names(compounds.unique) <- NULL
alerts <- c()
for (i in 1:nrow(compounds.unique)){
  accession <- compounds.unique$accession[i]
  alert <- data.subset$alert[which(data.subset$id == accession)]
  alerts <- c(alerts, alert)
  rm(accession, alert)
}
compounds.unique$alert <- alerts
row.names(compounds.unique) <- NULL # now we have a data frame with all chemicals in the subset and their smiles, inchikey and accession

# get unique inchikeys and create data.frame with these and their alerts
info.df <- compounds.unique %>%
  select("inchikey.pt1", "alert")
info.df <- distinct(info.df)
check.length <- unique(info.df$inchikey.pt1) # the length of check should be the same as the number of rows of info.df! If not, contradictory alert codes are present.

# get sampling
set.seed(seed) # make sure that the analysis is repeatable (& reproducible)
train.index <- createDataPartition(info.df$alert, p = 0.7, list = FALSE) # random sampling but with attempt to balance the class distributions within the splits.
inchikeys.train <- info.df[train.index, ] # create training set
inchikeys.test <- info.df[-train.index, ] # create test set

# Select spectra that correspond to the inchikeys
spectra.train <- compounds.unique[which(compounds.unique$inchikey.pt1 %in% inchikeys.train$inchikey.pt1), "accession"]
spectra.test <- compounds.unique[which(compounds.unique$inchikey.pt1 %in% inchikeys.test$inchikey.pt1), "accession"]

# make data frame and transfer alert column into factor
data.subset <- as.data.frame(data.subset)
data.subset$alert <- as.factor(data.subset$alert)

# create row names using the spectrum id
row.names(data.subset) <- data.subset$id
data.subset$id <- NULL # remove column

# create test and training set
data.train <- data.subset[spectra.train, ] # create training set
data.test <- data.subset[spectra.test, ] # create test set
```

### Check unique chemicals in training and test set

``` r
# make groups of spectra with the same InChIKeys
inchikey.groups <- split(x = metadata.combined, f = metadata.combined$inchikey)

# check overlapping inchikeys
training <- metadata.combined[which(metadata.combined$accession %in% row.names(data.train)),]
testing <- metadata.combined[which(metadata.combined$accession %in% row.names(data.test)),]
intersect(training$inchikey, testing$inchikey) # if this value is 0, pass

# check overlapping inchikeys (only first part of the inchikey)
training <- metadata.combined[which(metadata.combined$accession %in% row.names(data.train)),]
testing <- metadata.combined[which(metadata.combined$accession %in% row.names(data.test)),]
intersect(training$inchikey.pt1, testing$inchikey.pt1) # if this value is 0, pass
```

### Save data

``` r
# save training and test set
write.csv(x = data.train, file = paste0("data_train_organophosphorus_", data.type, ".csv"))
write.csv(x = data.test, file = paste0("data_test_organophosphorus_", data.type, ".csv"))
```
