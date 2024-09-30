Prepare training and test set for aromatic amine classification -
combination
================

This is an R Markdown document with relevant chunks from the code used
to generate a training and test set for the classification of MS2
spectra with an aromatic amine structural alert. Here, the combination
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

# select spectra of molecules with the alert of interest, some might contain other alerts as well. And include spectra with 'NA'.
data.subset.alert <- data %>% filter(grepl("TA430", data$alert) | grepl("TA386", data$alert) | grepl("TA394", data$alert) |
  grepl("TA322", data$alert) | grepl("TA330", data$alert) | grepl("TA341", data$alert) |
  grepl("TA398", data$alert) | grepl("TA355", data$alert) | grepl("TA385", data$alert) |
  grepl("TA441", data$alert)) #  |is.na(data$alert))
row.names(data.subset.alert) <- NULL
data.subset.alert$alert <- "alert"
data.subset.noalert <- data %>% filter(!grepl("TA430", data$alert) & !grepl("TA386", data$alert) & !grepl("TA394", data$alert) &
  !grepl("TA322", data$alert) & !grepl("TA330", data$alert) & !grepl("TA341", data$alert) &
  !grepl("TA398", data$alert) & !grepl("TA355", data$alert) & !grepl("TA385", data$alert) &
  !grepl("TA441", data$alert))
row.names(data.subset.noalert) <- NULL

# get a random sample of spectra without the aromatic amine alert
set.seed(seed)
x <- data.subset.noalert[sample(nrow(data.subset.noalert), size = 7600, replace = FALSE), ]
x$alert <- "none"

# combine both datasets (spectra with alert, and random subset of spectra without the alert)
data.subset <- rbind(data.subset.alert, x)
data.subset <- data.subset %>% relocate(alert, .after = last_col())

# Remove original data set so we have more memory space free
rm(data, data.subset.alert, data.subset.noalert, x, data.f, data.nl)
```

### Tweak mislabelled spectra

The last code chunk (‘check unique chemicals in training and test set’)
of this document checks whether there is no overlap between training and
test sets in terms of SMILES, InChIKey and the first part of the
InChIKey. During processing, it appeared that for some spectra with
identical first fourteen characters of the InChIKey, ended up in both
data sets due to different labels of the alert. These were visually
inspected and corrected in this code chunck.

``` r
# Adjust instance with InChIkey NFYYATWFXNPTRM-ICFOCEFXSA-N -> should be coded with alert instead of no alert
spectra <- c(
  "MSBNK-RIKEN-PR300239", "MSBNK-RIKEN-PR300249", "MSBNK-RIKEN-PR300254", "MSBNK-RIKEN-PR300259",
  "MSBNK-RIKEN-PR300269", "MSBNK-RIKEN-PR300274", "MSBNK-RIKEN-PR300289", "MSBNK-RIKEN-PR300294"
)
data.subset[which(data.subset$id %in% spectra), "alert"] <- "alert"

# Adjust instance with InChIkey JMIAZDVHNCCPDM-DQDWJNSRSA-N, JMIAZDVHNCCPDM-DAFCLMLCSA-N, JMIAZDVHNCCPDM-QLMFUGSGSA-N,
# JMIAZDVHNCCPDM-QLMFUGSGSA-N, JMIAZDVHNCCPDM-ZUNJVLJPSA-N, JMIAZDVHNCCPDM-PMJXBNNDSA-N -> should be coded with alert instead of no alert
spectra <- c(
  "MSBNK-RIKEN-PR300002", "MSBNK-RIKEN-PR300017", "MSBNK-RIKEN-PR300032", "MSBNK-RIKEN-PR300037", "MSBNK-RIKEN-PR300051", "MSBNK-RIKEN-PR300063",
  "MSBNK-RIKEN-PR300073", "MSBNK-RIKEN-PR300083", "MSBNK-RIKEN-PR300093", "MSBNK-RIKEN-PR300098", "MSBNK-RIKEN-PR300103", "MSBNK-RIKEN-PR300113",
  "MSBNK-RIKEN-PR300241", "MSBNK-RIKEN-PR300256", "MSBNK-RIKEN-PR300281", "MSBNK-RIKEN-PR300296", "MSBNK-RIKEN-PR300329", "MSBNK-RIKEN-PR300339",
  "MSBNK-RIKEN-PR300359", "MSBNK-RIKEN-PR300374", "MSBNK-RIKEN-PR300379", "MSBNK-RIKEN-PR300394", "MSBNK-RIKEN-PR300399", "MSBNK-RIKEN-PR300424",
  "MSBNK-RIKEN-PR300434", "MSBNK-RIKEN-PR300438", "MSBNK-RIKEN-PR300443", "MSBNK-RIKEN-PR300450", "MSBNK-RIKEN-PR300454", "MSBNK-RIKEN-PR300466",
  "MSBNK-RIKEN-PR300474", "MSBNK-RIKEN-PR300482", "MSBNK-RIKEN-PR300486", "MSBNK-RIKEN-PR310618", "MSBNK-RIKEN-PR310619", "MSBNK-RIKEN-PR310620"
)
data.subset[which(data.subset$id %in% spectra), "alert"] <- "alert"

# Adjust instance with InChIkey DAXYUDFNWXHGBE-VKCGGMIFSA-N and DAXYUDFNWXHGBE-KAXDATADSA-N -> should be coded with alert instead of no alert
spectra <- c(
  "MSBNK-RIKEN-PR300186", "MSBNK-RIKEN-PR300262", "MSBNK-RIKEN-PR300297"
)
data.subset[which(data.subset$id %in% spectra), "alert"] <- "alert"
```

### Removal of empty columns, near zero variance predictors and duplicates

``` r
set.seed(seed) 

# Remove empty features (columns), near zero variance predictors and duplicates ----------
x <- colSums(data.subset[, 2:(ncol(data.subset) - 1)], na.rm = TRUE) # count sums per column, do not consider column 'id' and 'alert'
zeros <- x[which(x == 0)] # only select non-occurring fragments / neutral losses
data.subset <- data.subset[, names(zeros) := NULL] # remove these from the data set
rm(x, zeros, spectra) # clean environment

# identify near zero variance predictors
zerovariance <- nearZeroVar(
  x = data.subset,
  freqCut = 95 / 5, # the cutoff for the ratio of the most common value to the second most common value
  saveMetrics = TRUE, # return data frame with specs per predictor
)

# remove near zero variance predictors from data set
zerovariance$freqRatioRounded <- round(zerovariance$freqRatio, digits = 0)
zerovariance.name <- rownames(zerovariance[which(zerovariance$nzv == TRUE & zerovariance$freqRatioRounded == (nrow(data.subset) - 1)), ])
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
metadata <- read.csv(file = "data\\MassBankEU_202206_metadata_subset.csv")
metadata.id.smiles.inchikey <- metadata %>%
  select("accession", "smiles", "inchikey")
metadata.id.smiles.inchikey$inchikey.pt1 <- substr(metadata.id.smiles.inchikey$inchikey, start = 1, stop = 14)

# get unique compounds in the subset, retrieve accession, smiles and inchikey and couple with alert
compounds.unique <- metadata.id.smiles.inchikey[which(metadata.id.smiles.inchikey$accession %in% data.subset$id), ]
row.names(compounds.unique) <- NULL
alerts <- c()
for (i in 1:nrow(compounds.unique)) {
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
# make groups of spectra with the same smiles
smiles.groups <- split(x = metadata, f = metadata$smiles)
inchikey.groups <- split(x = metadata, f = metadata$inchikey)

# check overlapping smiles
training <- metadata.id.smiles.inchikey[which(metadata.id.smiles.inchikey$accession %in% row.names(data.train)), ]
testing <- metadata.id.smiles.inchikey[which(metadata.id.smiles.inchikey$accession %in% row.names(data.test)), ]
intersect(training$smiles, testing$smiles) # if this value is 0, pass

# check overlapping inchikeys
training <- metadata.id.smiles.inchikey[which(metadata.id.smiles.inchikey$accession %in% row.names(data.train)), ]
testing <- metadata.id.smiles.inchikey[which(metadata.id.smiles.inchikey$accession %in% row.names(data.test)), ]
intersect(training$inchikey, testing$inchikey) # if this value is 0, pass

# check overlapping inchikeys (only first part of the inchikey)
training <- metadata.id.smiles.inchikey[which(metadata.id.smiles.inchikey$accession %in% row.names(data.train)), ]
testing <- metadata.id.smiles.inchikey[which(metadata.id.smiles.inchikey$accession %in% row.names(data.test)), ]
intersect(training$inchikey.pt1, testing$inchikey.pt1) # if this value is 0, pass
```

### Save data

``` r
# save training and test set
write.csv(x = data.train, file = paste0("data_train_aromaticamine_", data.type, ".csv"))
write.csv(x = data.test, file = paste0("data_test_aromaticamine_", data.type, ".csv"))
```
