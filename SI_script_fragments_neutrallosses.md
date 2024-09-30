Data pre-processing: obtain fragments and neutral losses
================

This is an R Markdown document with relevant chunks from the code used
to format [MassBankEU](https://massbank.eu/MassBank/) data into data
sets suitable for machine learning, containing fragments and neutral
losses.

``` r
# Load required packages
library(dplyr) # for data handling
```

### Define functions for obtaining neutral losses and fragments

``` r
# Function to calculate neutral losses
getNL <- function(input) {
  
  # retrieve MS2 fragments, get corresponding parent mass
  ms2 <- input$ms2
  precursor.mass <- input$precursor.mz
  fragments <- ms2$mz[which(ms2$rel.intensity > 50)] # only include fragments with a relative intensity > 50
  
  # calculate neutral loss, remove neutral losses or 0 (meaning that only the precursor has been measured)
  neutral.loss <- precursor.mass - fragments
  neutral.loss <- neutral.loss[which(neutral.loss > 0)]
  
  # round to 1 digits, make sure unique neutral losses are present 
  output <- unique(round(neutral.loss, digits = 1))
  
  return(output)
}

# Define function to get fragments
getFragments <- function(input) {
  
  # retrieve MS2 fragments, get corresponding parent mass
  ms2 <- input$ms2
  precursor.mass <- input$precursor.mz
  fragments <- ms2$mz[which(ms2$rel.intensity > 50)] # only include fragments with a relative intensity > 50
  
  # remove fragments larger than precursor
  fragments <- fragments[which(fragments < precursor.mass)]
  
  # round to 1 digits, make sure unique fragments are present 
  output <- unique(round(fragments, digits = 1))
  
  return(output)
}
```

### Load data

MassBank data was filtered and formatted in a long list with accession
code, precursor m/z, exact mass and MS2 data.

``` r
mbank.ms2 <- readRDS(file = "MassBank_2022-06_subset_MS2only_adductH_prec.RDS")
mbank.ms2[[1]] # inspect first element of list
```

    ## $accession
    ## [1] "MSBNK-AAFC-AC000001"
    ## 
    ## $precursor.mz
    ## [1] 179.0697
    ## 
    ## $exactmass
    ## [1] 178.063
    ## 
    ## $ms2
    ##         mz intensity rel.intensity
    ## 1 133.0648 21905.332           225
    ## 2 151.0754  9239.897            94
    ## 3 155.9743 10980.890           112
    ## 4 161.0597 96508.438           999
    ## 5 179.0703 72563.875           750

### Obtain neutral losses

``` r
# Apply get neutral loss function to list (progressbar included)
temp.nl <- pbapply::pbsapply(mbank.ms2, FUN = getNL, simplify = FALSE, USE.NAMES = TRUE)

# Remove spectra without any neutral losses
temp.nl <- temp.nl[lengths(temp.nl) > 0]

# Create empty matrix with neutral losses as column names and mbank accession as rownames
cnl <- matrix(nrow = length(temp.nl), ncol = length(nl.unique))
colnames(cnl) <- as.character(nl.unique)
rownames(cnl) <- names(temp.nl)

# clean environment to make memory free
rm(nl.plot, mbank.ms2, nl.unique)

# Fill matrix
for (n in 1:length(temp.nl)) {
  
  # retrieve data frame and spectrum name
  clipb <- as.character(temp.nl[[n]])
  name <- names(temp.nl[n])
  
  # assign '1' to cells that correspond to the neutral loss
  cnl[name, clipb] <- 1
  
  # clean up environment and print progress so that user can keep track of the loop
  rm(clipb, name)
  print(n)
}

# Replace all NA values with zero
cnl[is.na(cnl)] <- 0
```

### Obtain fragments

``` r
# Apply get fragments function to list (progressbar included)
temp.frag <- pbapply::pbsapply(mbank.ms2, FUN = getFragments, simplify = FALSE, USE.NAMES = TRUE)

# Remove spectra without any fragments
temp.frag <- temp.frag[lengths(temp.frag) > 0]

# Create empty matrix with neutral losses as column names and mbank accession as rownames
fragments <- matrix(nrow = length(temp.frag), ncol = length(frag.unique))
colnames(fragments) <- as.character(frag.unique)
rownames(fragments) <- names(temp.frag)

# clean environment to make memory free
rm(frag.plot, mbank.ms2, frag.unique)

# Fill matrix
for (n in 1:length(temp.frag)) {
  
  # retrieve data frame and spectrum name
  clipb <- as.character(temp.frag[[n]])
  name <- names(temp.frag[n])
  
  # assign '1' to cells that correspond to the fragments
  fragments[name, clipb] <- 1
  
  # clean up environment and print progress so that user can keep track of the loop
  rm(clipb, name)
  print(n)
}

# Replace all NA values with zero
fragments[is.na(fragments)] <- 0
```

### Combine with ToxAlerts search results

To each matrix, an extra column was added with the ToxAlerts alert code,
using an overview data frame with MassBank metadata.

``` r
library(data.table) # here we're adding another type of data to the matrix, so data.table needs to be used
mbank.metadata <- readRDS(file = "MassBank_metadata_2022-06_subset_alerts.RDS")

# Combine ToxAlerts results with neutral losses matrix
cnl <- as.data.table(cnl, keep.rownames = "id", na.rm = FALSE)
cnl[, alert := "none"] # add extra column were alerts can be specified

# Add alert to data.table with neutral losses
for (i in 1:nrow(cnl)) {
  accession <- cnl[i, id]
  alert.id <- mbank.metadata$alert[which(mbank.metadata$accession == accession)]
  cnl[i, alert := alert.id]
  print(i)
  rm(accession, alert.id)
}

# Combine ToxAlerts with fragments matrix
fragments <- as.data.table(fragments, keep.rownames = "id", na.rm = FALSE)
fragments[, alert := "none"] # add extra column were alerts can be specified

# Add alert to data.table with fragments
for (i in 1:nrow(fragments)) {
  accession <- fragments[i, id]
  alert.id <- mbank.metadata$alert[which(mbank.metadata$accession == accession)]
  fragments[i, alert := alert.id]
  print(i)
  rm(accession, alert.id)
}
```

### Split both data sets in positive and negative ionization

``` r
# Select accession codes for positive and negative ionization mode
pos <- mbank.metadata$accession[which(mbank.metadata$adduct.prec == "[M+H]+" | mbank.metadata$adduct.focus == "[M+H]+")]
neg <- mbank.metadata$accession[which(mbank.metadata$adduct.prec == "[M-H]-" | mbank.metadata$adduct.focus == "[M-H]-")]
```

Split neutral losses data set.

``` r
## NEUTRAL LOSS
colnames(cnl) <- paste0("nl", colnames(cnl)) # add 'nl' since R cannot handle column names starting with numbers in machine learning
colnames(cnl)[colnames(cnl) == "nlalert"] <- "alert" # change column 'alert' to itself.
colnames(cnl)[colnames(cnl) == "nlid"] <- "id" # change column 'id' to itself.

# Split neutral loss dataset into positive and negative mode
cnl.pos <- cnl[which(cnl$id %in% pos),]
row.names(cnl.pos) <- NULL
cnl.neg <- cnl[which(cnl$id %in% neg),]
row.names(cnl.neg) <- NULL

# Remove columns for neutral losses that do not occur in this subset
# positive
n <- colSums(cnl.pos[, 2:(which(colnames(cnl.pos) == "alert")-1)], na.rm = TRUE) #count sums per column, do not consider column 'id' and 'alert'
cnl.zero <- n[which(n == 0)] # only select non-occurring nls
cnl.pos <- cnl.pos[, names(cnl.zero):= NULL] # remove these from the data set
rm(n, cnl.zero) # clean environment

# negative
n <- colSums(cnl.neg[, 2:(which(colnames(cnl.neg) == "alert")-1)], na.rm = TRUE) #count sums per column, do not consider column 'id' and 'alert'
cnl.zero <- n[which(n == 0)] # only select non-occurring nls
cnl.neg <- cnl.neg[, names(cnl.zero):= NULL] # remove these from the data set
rm(n, cnl.zero) # clean environment
```

Split fragments data set.

``` r
## FRAGMENTS
colnames(fragments) <- paste0("f", colnames(fragments)) # add 'f' since R cannot handle column names starting with numbers in machine learning
colnames(fragments)[colnames(fragments) == "falert"] <- "alert" # change column 'alert' to itself.
colnames(fragments)[colnames(fragments) == "fid"] <- "id" # change column 'id' to itself.

# Split fragments dataset into positive and negative mode
fragments.pos <- fragments[which(fragments$id %in% pos),]
row.names(fragments.pos) <- NULL
fragments.neg <- fragments[which(fragments$id %in% neg),]
row.names(fragments.neg) <- NULL

# Remove columns for fragments that do not occur in this subset
# positive
n <- colSums(fragments.pos[, 2:(which(colnames(fragments.pos) == "alert")-1)], na.rm = TRUE) #count sums per column, do not consider column 'id' and 'alert'
fragments.zero <- n[which(n == 0)] # only select non-occurring fragments
fragments.pos <- fragments.pos[, names(fragments.zero):= NULL] # remove these from the data set
rm(n, fragments.zero) # clean environment

# negative
n <- colSums(fragments.neg[, 2:(which(colnames(fragments.neg) == "alert")-1)], na.rm = TRUE) #count sums per column, do not consider column 'id' and 'alert'
fragments.zero <- n[which(n == 0)] # only select non-occurring fragments
fragments.neg <- fragments.neg[, names(fragments.zero):= NULL] # remove these from the data set
rm(n, fragments.zero) # clean environment
```
