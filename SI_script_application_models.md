Application of structural alert models to NTS data
================

This is an R Markdown document with the code to classify non-target
screening feature groups for the presence of an aromatic amine and an
organophosphorus structural alert, based on their MS2 spectra.

``` r
# Load required packages
library(patRoon) # for NTS feature finding and processing

# Load required functions
source("structuralalertFunctions.R")

# define seed
seed <- 234
set.seed(seed)
```

### Initialization

``` r
# Convert C18 positive mode files from .raw to .mzML (only run once)
convertMSFiles(
  files = paste0("raw data\\", list.files("raw data\\")),
  algorithm = "pwiz",
  dirs = FALSE,
  from = "thermo",
  to = "mzML"
)

# read file with analysis info
anaInfo.c18pos <- readxl::read_excel(path = "raw data\\analysisinfo.xlsx", sheet = "patroon_final", col_names = TRUE)
```

### Feature finding, grouping and filtering

``` r
set.seed(seed) 

# Find and group all features with OpenMS
fList <- findFeatures(anaInfo.c18pos, "openms", noiseThrInt = 10000, chromSNR = 5, chromFWHM = 10, minFWHM = 1, maxFWHM = 33.154, mzPPM = 5.6)

# Group and align features between analyses
fGroups <- groupFeaturesOpenMS(fList, rtalign = TRUE, maxAlignRT = 30, maxGroupRT = 6)

# Calculate peak qualities
fGroups.pq <- calculatePeakQualities(fGroups) # calculate peak qualities for all features and groups

# Remove feature groups with a GaussianSimilarityScore < 0.3
fGroups <- filter(fGroups.pq, featQualityRange = list(GaussianSimilarityScore = c(0.3, Inf)))

# Basic rule based filtering
fGroups <- filter(fGroups,
  absMinIntensity = 10000, relMinReplicateAbundance = 0.66,
  maxReplicateIntRSD = NULL, blankThreshold = 10, removeBlanks = TRUE,
  retentionRange = c(162, 1620), mzRange = NULL
)
```

### Obtain MS peak lists

``` r
set.seed(seed)

avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.0001)

mslists <- generateMSPeakLists(fGroups, "mzr",
  maxMSRtWindow = 5, precursorMzWindow = 4,
  topMost = 5, # only extract MS peak lists from a maximum of analyses with highest intensity
  avgFeatParams = avgMSListParams,
  avgFGroupParams = avgMSListParams
)
```

### Calculate formula candidates and fingerprints using SIRIUS and CSI:FingerID, generate components

``` r
# Obtain formulas and fingerprints using SIRIUS
formulas.sirius <- generateFormulas(fGroups, mslists, "sirius",
                                    relMzDev = 5, adduct = "[M+H]+", elements = "CHNOPSClBrF",
                                    calculateFeatures = FALSE, # has to be set to FALSE in order to be able to calculate fingerprints
                                    profile = "orbitrap",
                                    getFingerprints = TRUE, # calculate fingerprints using SIRIUS-CSI:FingerID
                                    cores = 7, # specify number of cores, leave one free for other stuff
                                    token = token.sirius,
                                    projectPath = "sirius"
)

# Generate components
componCAM <- generateComponents(fGroups, "camera", ionization = "positive")
fGroupsSel <- selectIons(fGroups, componCAM, "[M+H]+")
```

### Load MassBank and do suspect screening

``` r
# Load and filter MassBank
mslibrary <- loadMSLibrary("MassBank_NIST_2024_06.msp", "msp") # load MassBank.eu MSP library (downloaded from https://github.com/MassBank/MassBank-data/releases, release version 2024.06 (june 5), downloaded on 24th of july 2024)
unique(records(mslibrary)[["Instrument_type"]]) # get all unique instrument types
mslibraryF <- patRoon::filter(mslibrary, properties = list(Instrument_type = c("LC-ESI-ITFT", "ESI-ITFT", "ESI-QTOF", "LC-ESI-QTOF"))) # only select orbitrap and QTOF data

# calculate compounds
compounds.mb <- generateCompounds(fGroupsSel, mslists, "library", MSLibrary = mslibraryF, minSim = 0.5
)
compounds.mb <- addFormulaScoring(compounds.mb, formulas.sirius, updateScore = TRUE)
```

### Application of structural alert models

``` r
# obtain raw MS data from mslists object
ms.data <- unlist(mslists@peakLists, recursive = FALSE)

# retrieve fragments and neutral losses 
ms.data.frag <- pbapply::pbsapply(ms.data, FUN = getFragments, simplify = FALSE, USE.NAMES = TRUE) # get fragments
ms.data.frag <- ms.data.frag[lapply(ms.data.frag, length) > 0] # remove feature groups without fragments
ms.data.nl <- pbapply::pbsapply(ms.data, FUN = getNeutralLosses, simplify = FALSE, USE.NAMES = TRUE) # get neutral losses
ms.data.nl <- ms.data.nl[lapply(ms.data.nl, length) > 0] # remove feature groups without neutral losses

# Get predictions for the organophosphorus structural alerts
set.seed(seed)
result.organophosphorus <- getPredictionAlert(frags = ms.data.frag, nls = ms.data.nl, alert = "organophosphorus")

# Get predictions for the aromatic amine structural alerts
set.seed(seed)
result.aromaticamine <- getPredictionAlert(frags = ms.data.frag, nls = ms.data.nl, alert = "aromaticamine")

# Filter fGroups with organophosphorus alert
library(dplyr) # for filtering
fGroups.organophosphorus <- result.organophosphorus %>% filter(alert > 0.5)
fGroups.aromaticamine <- result.aromaticamine %>% filter(alert > 0.5)
```

### Report results

``` r
report(fGroupsSel,
       MSPeakLists = mslists, formulas = formulas.sirius, compounds = compounds.mb,
       components = NULL, openReport = FALSE, path = paste0("report//", Sys.Date(), "_report_application.html")
)
```
