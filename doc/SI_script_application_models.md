Application of structural alert models to NTS data
================

This is an R Markdown document with the code to classify non-target
screening feature groups for the presence of an aromatic amine and an
organophosphorus structural alert, based on their MS2 spectra.

``` r
# Load required packages
library(patRoon) # for NTS feature finding and processing, using version 4.3.3

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

### Componentization

``` r
# cluster feature groups with similar chromatographic elution profiles and annotate by known chemical rules
components <- generateComponents(
  fGroups = fGroups,
  algorithm = "camera",
  ionization = "positive",
  onlyIsotopes = FALSE,
  minSize = 1,
  relMinReplicates = 0.3
)

# select feature groups with adducts [M+H]+, and assign [M+H]+ to those where the adduct is unknown
fGroups <- selectIons(fGroups, components, prefAdduct = "[M+H]+")
```

### Obtain MS peak lists

``` r
set.seed(seed)
avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.0001)
mslists <- generateMSPeakLists(
  fGroups = fGroups,
  algorithm = "mzr",
  maxMSRtWindow = 5,
  precursorMzWindow = 0.8,
  avgFeatParams = avgMSListParams,
  avgFGroupParams = avgMSListParams,
  topMost = NULL
)

# Rule based filtering of MS peak lists (for screening MassBank etc.)
mslistsF <- patRoon::filter(
  mslists,
  withMSMS = TRUE,
  absMSIntThr = NULL,
  absMSMSIntThr = NULL,
  relMSIntThr = NULL,
  relMSMSIntThr = NULL,
  topMSPeaks = NULL,
  topMSMSPeaks = 25,
  deIsotopeMS = FALSE,
  deIsotopeMSMS = FALSE
)
```

### Calculate formula candidates and fingerprints using SIRIUS and CSI:FingerID

``` r
# Calculate formula candidates and fingerprints using SIRIUS (version 5.8.2)
formulas.sirius <- generateFormulas( # use formulas for fingerprint predictions
  fGroups = fGroups,
  MSPeakLists = mslists,
  algorithm = "sirius",
  relMzDev = 5, # in ppm
  adduct = "[M+H]+",
  elements = "CHNOPSClBrF",
  profile = "orbitrap",
  database = "pubchem", # use PubChem as database for retrieval of formula candidates
  cores = 7,
  getFingerprints = TRUE,
  topMost = 5,
  calculateFeatures = FALSE,
  absAlignMzDev = 0.0015, # in Da
  extraOptsFormula = "--ppm-max-ms2=50",
  verbose = TRUE,
  login = "interactive"
)
```

### Load MassBank and do suspect screening

``` r
# Load and filter MassBank
mslibrary <- loadMSLibrary("MassBank_NIST_2024_06.msp", "msp") # load MassBank.eu MSP library (downloaded from https://github.com/MassBank/MassBank-data/releases, release version 2024.06 (june 5), downloaded on 24th of july 2024)
unique(records(mslibrary)[["Instrument_type"]]) # get all unique instrument types
mslibraryF <- patRoon::filter(mslibrary, properties = list(Instrument_type = c("LC-ESI-ITFT", "ESI-ITFT", "ESI-QTOF", "LC-ESI-QTOF"))) # only select orbitrap and QTOF data

# Find candidate structures using MassBank Europe (only look into [M+H]+ adducts)
compounds.mb <- generateCompounds(fGroups,
  MSPeakLists = mslistsF,
  algorithm = "library",
  MSLibrary = mb,
  minSim = 0.5,
  absMzDev = 0.05, # in Da
  spectrumType = "MS2",
  adduct = "[M+H]+",
  checkIons = "adduct"
)
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

# get predictions for the organophosphorus structural alerts
set.seed(seed)
result.organophosphorus <- getPredictionAlert(frags = ms.data.frag, nls = ms.data.nl, alert = "organophosphorus")
result.organophosphorus$name <- row.names(result.organophosphorus)
result.organophosphorus$group <- stringr::str_split_fixed(string = result.organophosphorus$name, pattern = "\\.", n = 2)[, 2]

# get predictions for the aromatic amine structural alerts
set.seed(seed)
result.aromaticamine <- getPredictionAlert(frags = ms.data.frag, nls = ms.data.nl, alert = "aromaticamine")
result.aromaticamine$name <- row.names(result.aromaticamine)
result.aromaticamine$group <- stringr::str_split_fixed(string = result.aromaticamine$name, pattern = "\\.", n = 2)[, 2]

# select fGroups with prediction > 0.5
fGroups.aromaticamine <- result.aromaticamine %>% # filter(alert > 0.5) %>%
  select(alert, group)
fGroups.organophosphorus <- result.organophosphorus %>% # filter(alert > 0.5) %>%
  select(alert, group)

# find max predicted probability for every group - aromatic amine
fGroups.aromaticamine.clean <- data.frame(
  fGroup = unique(fGroups.aromaticamine$group),
  alert.max = NA
)

for (i in 1:nrow(fGroups.aromaticamine.clean)) {
  fGroup.aa <- fGroups.aromaticamine.clean$fGroup[i]
  alert.max <- max(fGroups.aromaticamine$alert[which(fGroups.aromaticamine$group == fGroup.aa)])
  fGroups.aromaticamine.clean$alert.max[i] <- alert.max
  rm(fGroup.aa, alert.max)
}

# find max predicted probability for every group - organophosphorus
fGroups.organophosphorus.clean <- data.frame(
  fGroup = unique(fGroups.organophosphorus$group),
  alert.max = NA
)

for (j in 1:nrow(fGroups.organophosphorus.clean)) {
  fGroup.op <- fGroups.organophosphorus.clean$fGroup[j]
  alert.max <- max(fGroups.organophosphorus$alert[which(fGroups.organophosphorus$group == fGroup.op)])
  fGroups.organophosphorus.clean$alert.max[j] <- alert.max
  rm(fGroup.op, alert.max)
}
```
