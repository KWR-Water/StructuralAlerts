#------------------------------------------------------------------------------#
# Functions for data pre-processing and prediction of organophosphorus         #
# structural alert or aromatic amine structural alert                          #
#                                                                              #
# Developed by Nienke Meekel                                                   #
# KWR Water Research Institute, Vrije Universiteit Amsterdam                   #
#                                                                              #
# Current version: 9 August 2024                                               #
# developed on R version 4.3.3                                                 #
# ---------------------------------------------------------------------------- #

# Load required model
model.organophosphorus <- readRDS("model_organophosphorus.RDS")
model.aromaticamine <- readRDS("model_aromaticamine.RDS")

# Function to get fragments
getFragments <- function(msdata) {
  # get data frame
  ms2 <- msdata[["MSMS"]]

  if (is.null(ms2)) {
    output <- NULL
  } else {
    # calculate relative intensities
    ms2$rel.intensity <- (ms2$intensity / max(ms2$intensity)) * 999

    # retrieve MS2 fragments, get corresponding parent mass
    precursor <- ms2$mz[which(ms2$precursor == TRUE)]
    fragments <- ms2$mz[which(ms2$rel.intensity > 50)] # only include fragments with a relative intensity > 50

    # remove fragments larger than precursor
    fragments <- fragments[which(fragments < precursor)]

    # round to 1 digits, make sure unique fragments are present
    output <- unique(round(fragments, digits = 1))
  }
  return(output)
}

# Function to get fragments
getNeutralLosses <- function(msdata) {
  # get data frame
  ms2 <- msdata[["MSMS"]]

  if (is.null(ms2)) {
    output <- NULL
  } else {

    # calculate relative intensities
    ms2$rel.intensity <- (ms2$intensity / max(ms2$intensity)) * 999

    # retrieve MS2 fragments, get corresponding parent mass
    precursor <- ms2$mz[which(ms2$precursor == TRUE)]
    fragments <- ms2$mz[which(ms2$rel.intensity > 50)] # only include fragments with a relative intensity > 50

    # remove fragments larger than precursor
    fragments <- fragments[which(fragments < precursor)]

    # calculate neutral losses
    neutral.losses <- precursor - fragments

    # round to 1 digits, make sure unique fragments are present
    output <- unique(round(neutral.losses, digits = 1))
  }

  return(output)
}

# Function to create a matrix that is suitable for prediction
getPredictionAlert <- function(frags, nls, alert) {
  if (alert == "organophosphorus") {
    ml.model <- model.organophosphorus
  } else if (alert == "aromaticamine") {
    ml.model <- model.aromaticamine
  }

  if (length(frags) != length(nls)) {
    stop("Different number of spectra for fragments and neutral losses.")
  }

  # get required input variables
  variables.model <- ml.model[["finalModel"]][["xNames"]]

  # define matrix
  variables <- matrix(nrow = length(frags), ncol = length(variables.model))
  colnames(variables) <- as.character(variables.model)
  rownames(variables) <- names(frags)

  # fill matrix
  for (n in 1:length(frags)) {
    if (names(frags)[n] != names(nls)[n]) {
      stop("spectrum is not the same")
    }

    # retrieve data frame and spectrum name from fragments
    temp <- paste0("f", as.character(frags[[n]]))
    tokeep <- temp[which(temp %in% colnames(variables))]
    name <- names(frags[n])

    # assign '1' to cells that correspond to the fragments
    variables[name, tokeep] <- 1

    # clean environment
    rm(temp, name, tokeep)

    # retrieve data frame and spectrum name from neutral losses
    temp <- paste0("nl", as.character(nls[[n]]))
    tokeep <- temp[which(temp %in% colnames(variables))]
    name <- names(nls[n])

    # assign '1' to cells that correspond to the neutral losses
    variables[name, tokeep] <- 1

    # clean up environment
    rm(temp, name, tokeep)
  }

  # replace all NA values with zero
  variables[is.na(variables)] <- 0

  # predict presence of alert
  outcome <- predict(ml.model, newdata = variables, type = "prob")
  outcome$alert <- round(outcome$alert, digits = 3)
  outcome$none <- round(outcome$none, digits = 3)

  return(outcome)
}
