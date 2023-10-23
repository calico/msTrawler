# Main file for analyzing isobaric tag proteomics data
# All other code should be called from this function.

#' Analyze Proteomics Data
#'
#' @param DF  A dataframe containing identifiers and
#'  intensities that matches the structure in the sample
#'  data.
#' @param sampleFile A csv file containing all the covariate
#'  values that map to individual samples in the study. An
#'  example file called sampleFile can be found in the data
#'  folder.
#' @param covariateFile A csv file containing all covariate
#'  metadata that is used to determine the type of analyses
#'  that will be performed.  See covariateFile in the data
#'  folder.
#' @param scaleSN The technical variance of an observation is modeled
#'  as the inverse of the square root of each SNR * some constant.
#'  scaleSN is that constant.  This value will be used for both
#'  imputation and modeling.
#' @param lod A value that determines the limit of reliability.
#'  If the value is between zero and one, then
#'  the threshold will be a proportion of the total scan
#'  intensity.  If lod >= 1, then the threshold will
#'  be set at a SNR = lod for all scans.
#' @imputePenalty A value that diminishes the impact of LOD
#'  imputations.  Each value below the LOD will be replaced
#'  by a random variate centered at LOD / 2.  The corresponding
#'  SNR, which is used to create weights for each observation,
#'  is also imputed to be LOD / (2 * imputePenalty).
#' @param minAbove It is possible that all of the scans will
#'  be below the limit of reliability.  Imputing an entire scan
#'  is senseless.  For small numbers of values above lod the
#'  decision is questionable.  If fewer than minAbove observations
#'  are above the lod, then the scan will be flagged for
#'  removal.
#' @param ssnFilter A cutoff for each scan based on the sum of
#'  signal-to-noise ratios across all channels.  This filter
#'  precedes all others, so if the a scan has less total signal
#'  than "ssnFilter", then no LOD imputation will occur.  If
#'  NULL then no SSN filtering will be done.
#' @param outlierCutoff Outliers defined by consensus ratios
#'  within each protein grouping will be flagged whenever
#'  the jacknife residuals exceed this value.  Flagged scans
#'  will be removed from the protein grouping, but will
#'  remain in the analysis by replacing the protein label
#'  with a peptide sequence.  Set this parameter to NULL to
#'  skip this form of outlier removal.
#' @param N_SUM If the number of observations from a protein
#'  within a sample, is less than N_SUM, then both the fluxes
#'  and corresponding ion counts will be aggregated into a single
#'  data point.  Setting this to a high value will activate the
#'  single-level msTrawler model for all proteins.
#' @param swapProtein Boolean variable that determines the
#'  handling of outlier scans.  When a scan includes an observation
#'  that exceeds the outlierCutoff, we will either remove the scan
#'  from the analysis, or replace the protein label with the
#'  corresponding peptide label.  This lets us look for relationships
#'  that may be indicative of interesting post translational modifications.
#'  When swapProtein == TRUE, we perform the label swap.  Otherwise,
#'  we remove the scan.
#' @param maxPep An integer denoting the maximum number of scans
#'  to use from a single sample.  In plasma many proteins have
#'  hundreds of scans that can result in an out of memory error
#'  when working with large experiments.  Only the scans with the
#'  highest summed-signal-to-noise ratios will be used in the modeling.
#' @param colAdjust This parameter specifies how global column
#'  adjustments will be performed.  The parameter can
#'  take 3 types of input.  A value between 0 and 1
#'  will be interpreted as a standard deviation percentile.
#'  Only rows with a standard deviation percentile less than
#'  the entered value will be used to calculate column
#'  adjustment factors.  If colAdjust = NULL, then no
#'  column adjustment will be done at all.  The last
#'  type of acceptable parameter entry is a boolean
#'  vector with length equal to the number of rows in DF.
#'  Only the rows corresponding to 1s in this vector will
#'  be used to calculate column adjustment factors.
#' @param colRatios A numeric vector with length equal to the
#'  number of normalization factors.  The ratios in this
#'  vector determine the expected relationship between
#'  all of the columns.  Typically it will be a vector of 1's,
#'  but occasionally experiments will include columns at known
#'  dilutions, e.g. a 2x bridge channel.
#' @param dropContam Boolean variable that determines whether or
#'  not scans containing the string "contaminant" within
#'  the protein name will be removed prior to analysis.
#' @param dropReverse Boolean variable that determines whether or
#'  not scans containing the string "##" within
#'  the protein name will be removed prior to analysis.
#' @param peptideAnalysis A boolean parameter set to false by
#'  default.  When TRUE, protein labels will be replaced with
#'  peptide the peptide labels prior to calling any other functions.
#'  In words, we will estimate relationships for peptides instead of
#'  proteins.  This also triggers the single-sample msTrawler model
#'  regardless of the number of scans observed per peptide.
#' @param minRE The number of levels required to fit a random intercept.
#' @param timeDiff A boolean parameter that determines whether or
#'  not hypothesis tests for differential time trends will occur.
#'  The only reason to set this to false is if computational time
#'  is a concern and trends across categories are not of interest.
#' @export
msTrawl <- function(DF,
                    sampleFile = NULL,
                    covariateFile = NULL,
                    scaleSN = 1,
                    lod = 0.01,
                    imputePenalty = 1,
                    minAbove = 4,
                    ssnFilter = 20,
                    outlierCutoff = 3,
                    N_SUM = 3,
                    swapProtein = FALSE,
                    maxPep = 25,
                    colAdjust = 0.5,
                    colRatios = NULL,
                    dropContam = TRUE,
                    dropReverse = TRUE,
                    peptideAnalysis = FALSE,
                    minRE = 5,
                    timeDiff = TRUE) {

  # Step1.  Make sure the covariate file and sample file are consistent
  if (!is.null(covariateFile)) {
    # FORCES ROW ORDERING OF VARIABLES TO MATCH SAMPLE FILE COLUMNS (EXCLUDING BRIDGE AND SAMPLEID)
    covariateFile <- covariateFile[order(match(covariateFile$Covariate, colnames(sampleFile))), ]

    coVector <- covariateFile$Covariate
    sampVector <- colnames(sampleFile)
    nMatches <- sum(coVector %in% sampVector)
    if (nMatches < length(coVector)) {
      stop("Error: At least one of the covariate names does not match across files.")
    }
  }

    # Make sure there are no missing covariates
  if (!is.null(covariateFile)) {
    #ignore bridge channel
    bridgeI <- grep("BRIDGE", toupper(colnames(sampleFile)))
    if(length(bridgeI) > 0){
      nMissing <- sum(is.na(sampleFile[-which(sampleFile[ , bridgeI] == 1) , ]))
    }else{
      nMissing <- sum(is.na(sampleFile))
    }
    
    if (nMissing > 0) {
      stop("Missing values are not allowed in covariates")
    }
  }
  
  # If colAdjust is a vector, add it to DF so it will be included in the filtering steps
  if (length(colAdjust) > 1) {
    DF$colAdjust <- colAdjust
  }

  # Remove any junk
  if (dropReverse == TRUE) {
    revI <- grep("##", DF$Protein.ID)
    if (length(revI) > 0) {
      DF <- DF[-revI, ]
    }
  }
  if (dropContam == TRUE) {
    contamI <- grep("contam", DF$Protein.ID)
    if (length(contamI) > 0) {
      DF <- DF[-contamI, ]
    }
  }

  # Order the rows
  DF <- DF[order(DF$Protein.ID, DF$Plex, DF$Peptide), ]

  # Swap label columns if peptideAnalysis == TRUE
  if (peptideAnalysis) {
    DF$PA.Gene.Symbol <- paste0(DF$PA.Gene.Symbol, "___", DF$Protein.ID)
    DF$Protein.ID <- DF$Peptide
  }

  # Now create the two primary quantification matrices
  snMat <- DF[, grep("\\.Sn", colnames(DF))]
  iMat <- DF[, grep("Adjusted", colnames(DF))]

  # Make sure these have the same number of columns
  if (ncol(snMat) != ncol(iMat)) {
    stop("Error: The number of SNR columns
    does not equal the number of intensity columns. Often this occurs
    because of an unwanted column containing the string \".Sn\" in the
                                     column name.")
  }

  # Apply ssnFilter
  if (!is.null(ssnFilter)) {
    lowI <- which(apply(snMat, 1, sum) < ssnFilter)
    if (length(lowI) > 0) {
      DF <- DF[-lowI, ]
      snMat <- snMat[-lowI, ]
      iMat <- iMat[-lowI, ]
    }
  }

  # Now get rid of any rows with less than 2 non-zero intensities
  emptyRow <- which(apply(iMat, 1, function(x) {
    sum(x > 0, na.rm = T)
  }) < 2)
  if (length(emptyRow > 0)) {
    DF <- DF[-emptyRow, ]
    iMat <- iMat[-emptyRow, ]
    snMat <- snMat[-emptyRow, ]
  }


  ###################### Deal with the LOD#################
  lodRES <- tmtLOD(snMat, iMat, lod, minAbove, scaleSN, imputePenalty)
  # Remove rows with too few entries above LOD
  tooFew <- lodRES[[3]]
  iMat <- lodRES[[1]]
  snMat <- lodRES[[2]]

  if (sum(tooFew) > 0) {
    DF <- DF[-which(tooFew == TRUE), ]
    iMat <- iMat[-which(tooFew == TRUE), ]
    snMat <- snMat[-which(tooFew == TRUE), ]
  }


  ################## Remove outliers from protein groupings##########

  # Make sure the relevant ID columns are character vectors
  # so that replacement with peptide names does not result in NA
  DF$Protein.ID <- as.character(DF$Protein.ID)
  DF$Peptide <- as.character(DF$Peptide)
  uPlex <- unique(DF$Plex)
  uProt <- unique(DF$Protein.ID)
  for (j in 1:length(uPlex)) {
    for (i in 1:length(uProt)) {
      protIndex <- which(DF$Protein.ID == uProt[i] &
        DF$Plex == uPlex[j])
      
      #If the protein is missing in this plex, move on
      if(length(protIndex) == 0){next}
      
      subDat <- iMat[protIndex, , drop = F]
      subSn <- snMat[protIndex, , drop = F]
      subPep <- DF[protIndex, "Peptide"]
      subSSN <- apply(subSn, 1, sum)
      
      
      # Create outlier indicator matrix
      if (!is.null(outlierCutoff)) {
        outlierMat <- findOutliers(subDat, subSn, outlierCutoff, scaleSN)
      } else {
        outlierMat <- matrix(0, nrow = nrow(subDat), ncol = ncol(subDat))
      }

      # Figure out which rows have outliers
      outlierRows <- apply(outlierMat, 1, sum, na.rm = T)
      outlierIndex <- protIndex[which(outlierRows > 0)]

      # Now change the protein labels of each outlier
      if (length(outlierIndex) > 0) {
        if (swapProtein == TRUE) {
          DF[outlierIndex, "Protein.ID"] <- paste0(DF[outlierIndex, "Protein.ID"], "___", DF[outlierIndex, "Peptide"])
        } else {
          DF[outlierIndex, "Protein.ID"] <- "OUTLIER_REMOVE_AT_ONCE!"
        }
      }

      #Check our aggregation criteria.  
      n_after <- length(which(outlierRows == 0))
      if ((n_after < N_SUM & n_after > 1) |  (peptideAnalysis == TRUE & n_after > 1)) {
        
      newSnVec <- apply(subSn[which(outlierRows == 0), , drop=FALSE], 2, sum)
      newFluxVec <- apply(subDat[which(outlierRows == 0), , drop=FALSE], 2, sum)
      #Update the peptide name to indicate aggregation
      #Keep the first non-outlier row
      keepIndex <- protIndex[which(outlierRows == 0)][1]
      DF[keepIndex, "Peptide"] <- paste0(DF[protIndex[1], "Peptide"], "_SUM")
      #Update all the quant rows
      DF[keepIndex, grep("\\.Sn", colnames(DF))] <- newSnVec
      DF[keepIndex, grep("Adjusted", colnames(DF))] <- newFluxVec
      snMat[keepIndex, ] <- newSnVec
      iMat[keepIndex, ] <- newFluxVec
      #Flag the rest of the rows for removal
      DF[setdiff(protIndex, keepIndex), "Protein.ID"] <- "OUTLIER_REMOVE_AT_ONCE!"
      next
      }
      #End peptide aggregation step
      
      # Limit the remaining observations to maxPep # of scans
      if (nrow(subDat) - length(outlierIndex) > maxPep) {

        # Rank order each scan, smallest to largest, within each peptide
        # Data was ordered by protein, plex, peptide.  We are in a protein
        # and plex loop.  So peptides should be consecutive
        # also the unique() function returns strings in the order observed
        listOranks <- lapply(unique(subPep), function(x) {
          order(subSSN[which(subPep == x)], decreasing = TRUE)
        })
        rankVector <- do.call(c, listOranks)

        # Assign outliers an unrealistically high rank
        if (length(outlierIndex) > 0) {
          rankVector[which(outlierRows > 0)] <- 999999
        }

        snOrder <- order(rankVector, (-1 * subSSN)) # small rank, large SSN

        # Find the sub-indices of the top maxPep signals
        bigIndex <- head(snOrder, n = maxPep)
        # Find the indices of everything but the top maxPep
        smallIndex <- protIndex[setdiff(snOrder, bigIndex)]
        # Change the protein name of the small scans
        DF[smallIndex, "Protein.ID"] <- "TOO_MANY!!!"
      }
    } # End protein loop
  } # End plex loop

  outlierIndex <- grep("OUTLIER_REMOVE_AT_ONCE!", DF$Protein.ID)
  tooManyIndex <- grep("TOO_MANY!!!", DF$Protein.ID)
  removeIndex <- union(outlierIndex, tooManyIndex)
  if (length(removeIndex > 0)) {
    DF <- DF[-removeIndex, ]
    iMat <- iMat[-removeIndex, ]
    snMat <- snMat[-removeIndex, ]
  }



  ############## Implement a global column adjustment############

  if (is.null(colAdjust)) {
    normI <- iMat
  } else {
    if (length(colAdjust) == 1) {
      # Intialize boolean
      normBool <- rep(0, nrow(iMat))
      for (p_ in 1:length(uPlex)) {
        # Loop through each plex to find the median
        plexIndex <- which(DF$Plex == uPlex[p_])
        if (length(plexIndex) <= 1) {
          stop("You have a plex with <= 1 observation")
        } # An error occurs, in small subsets of data,
        # where only one observation exists in a plex.  This removes these cases
        # from the normalization

        pepSD <- apply(
          iMat[plexIndex, ],
          1, function(x) sd(log(x), na.rm = TRUE)
        )
        medSD <- quantile(pepSD, probs = colAdjust, na.rm = T)

        sdI <- which(pepSD < medSD)


        normBool[plexIndex[sdI]] <- 1
      }
    } else {
      normBool <- DF$colAdjust
    }

    if (length(colRatios) < 2) {
      ratios <- rep(1, length(uPlex) * ncol(iMat))
    }
    normed <- geoNorm(iMat, normBool, Plex = DF$Plex, ratios)
    normI <- normed[[1]]
    colnames(normI) <- paste0("norm_", colnames(normI))
    normFacs <- normed[[2]] # Return normalization factors
    # Write the column normalization table
    write.csv(normFacs,
      file = "ColAdjustmentFactors.csv",
      row.names = FALSE
    )
  } # End global column adjustment

  ################# Pre-processing complete. Read the covariate data##########

  # Figure out what we are dealing with.
  # First reduce the number of plexes
  usedPlexes <- intersect(uPlex, substring(sampleFile$SampleID, 1, regexpr("_", sampleFile$SampleID) - 1))

  # Is there a bridge channel?
  nBridges <- sum(sampleFile$Bridge)
  if (nBridges > 0) {
    bridgeMod <- TRUE

    # Quick sanity check on the number of bridge channels
    if (nBridges != length(usedPlexes)) {
      stop("Error: The number of plexes must match the number of bridge
         samples.")
    }
  } else {
    bridgeMod <- FALSE
  }


  # Initialize table names
  tableNames <- "Simple.csv"

  # Figure out the number of categorical covariates
  # Is there a bridge column?
  bridgeCol <- grep("BRIDGE", toupper(colnames(sampleFile)))
  if (length(bridgeCol) > 0) {
    coVars <- colnames(sampleFile)[-c(1, bridgeCol)]
  } else {
    coVars <- colnames(sampleFile)[-1]
  }

  # If sampleFile was NULL or had no covariates then set an intercept model
  if (length(coVars) == 0) {
    coVars <- "1"
  }

  covType <- covariateFile$Type
  covLevel <- covariateFile$Levels
  catIndex <- grep("FACTOR", toupper(covType))
  n_cat <- length(catIndex)


  # Initialize parameter for defining time table workflow
  tParm <- ""

  # Figure out the number of levels of each factor and the names
  if (n_cat > 0) {
    n_levels <- rep(0, length(catIndex)) # initialize vector
    factorNames <- as.character(covariateFile$Covariate[catIndex])

    # Make sure the covariate names match across files
    if (sum(factorNames %in% colnames(sampleFile)) != length(factorNames)) {
      stop("Error: The covariate names in the sampleFile do not match the
         column names in the sampleFile")
    }

    levelNames <- list()
    for (i in 1:n_cat) {
      if (length(bridgeCol) > 0) {
        realIndex <- which(sampleFile[, bridgeCol] == 0)
      } else {
        realIndex <- 1:nrow(sampleFile)
      }

      levelNames[[i]] <- as.character(unique(sampleFile[realIndex, factorNames[i]]))
      # Remove bridge, if it was entered as a factor level

      # bridgeIndex <- grep("BRIDGE", toupper(tempNames))
      # if(length(bridgeIndex) > 0){
      #   levelNames[[i]] <- tempNames[-bridgeIndex]
      # }else{
      #   levelNames[[i]] <- tempNames
      # }
      n_levels[i] <- length(levelNames[[i]])

      # Update table names
      tableNames <- c(tableNames, paste0(
        "Factor_",
        factorNames[i], "_",
        levelNames[[i]], ".csv"
      ))
    } # End level loop
  } else { # End "if" factors are present
    n_levels <- 0
    factorNames <- ""
  }

  # Figure out the number of continuous covariates and center them
  contIndex <- grep("CONTINUOUS", toupper(covariateFile$Type))
  n_cont <- length(contIndex)
  if (n_cont > 0) {
    contNames <- as.character(covariateFile$Covariate[contIndex])
    tableNames <- c(tableNames, paste0("Continuous_", contNames, ".csv"))
    for (i in 1:n_cont) {
      sampleFile[, contNames[i]] <- sampleFile[, contNames[i]] -
        mean(sampleFile[, contNames[i]], na.rm = T)
    }
  } else { # End "if" continuous
    contNames <- ""
  }

  # Figure out time related parameters
  timeIndex <- grep("TIME", toupper(covariateFile$Type))
  if (length(timeIndex) > 0) {
    # First rename the time variable for use in formulas
    colnames(sampleFile)[grep(
      covariateFile$Covariate[timeIndex],
      colnames(sampleFile)
    )] <- "Time"
    covariateFile$Covariate[timeIndex] <- "Time"

    # Make sure this variable is continuous
    if (!is.numeric(sampleFile$Time)) {
      stop("Error:  Your time variable is non-numeric")
    }

    # Remove Time variable from the covariate vector.
    coVars <- coVars[-timeIndex]
    # Start building a separate vector for time variables
    timeVars <- NULL

    # Now set the degree and circadian parameters
    timeDegree <- covariateFile$TimeDegree[timeIndex]
    if (timeDegree >= 1) {
      timeVars <- c(timeVars, "Time")
    }
    if (timeDegree >= 2) {
      timeVars <- c(timeVars, "Time2")
      sampleFile$Time2 <- sampleFile$Time^2
    }
    if (timeDegree == 3) {
      timeVars <- c(timeVars, "Time3")
      sampleFile$Time3 <- sampleFile$Time^3
    }
    # Look for a factor that separates time trends
    # Find position of time category in covariate file
    tCatIndex <- which(covariateFile$TimeCategory == 1)
    # Find position of time category in factor list
    tCatFactorIndex <- match(covariateFile$Covariate[tCatIndex], factorNames)
    if (length(tCatIndex) > 1) {
      stop("Only one time category is allowed")
    }
    if(length(tCatIndex) == 0){
      #Force a value greater than zero to always indicate the presence of a time category
      tCatIndex <- 0
      tCatFactorIndex <- 0
    } 
    
    if (tCatIndex > 0) {
      tCatName <- covariateFile$Covariate[tCatIndex]
      timeLevelN <- n_levels[tCatFactorIndex]
      tCatLevels <- levelNames[[tCatFactorIndex]]
      timeTableNames <- paste0("Time_", tCatLevels, ".csv")
      # There are three possible relevant parameterizations for time
      tParm <- "Category"
    } else { # This means there is no time category but there are time covariates
      timeLevelN <- 1
      timeTableNames <- "Time.csv"
      if (n_cont > 0) {
        tParm <- "Continuous"
      } else {
        tParm <- "Time"
      }
    }
    circadian <- covariateFile$Circadian[timeIndex]
  } else {
    timeDegree <- 0
    timeVars <- NULL
    timeLevelN <- 0
    tCatIndex <- 0
    tCatFactorIndex <- 0
    circadian <- 0
    tParm <- ""
  }

  # Create circadian covariates
  if (circadian != 0) {
    sampleFile$Sine <- sin((2 * pi / 24) * sampleFile$Time)
    sampleFile$Cosine <- cos((2 * pi / 24) * sampleFile$Time)
    timeVars <- c(timeVars, "Sine", "Cosine")
  }

  # Is this a longitudinal model?
  IDindex <- grep("ID", toupper(covariateFile$Type))
  if (length(IDindex) > 0) {
    IDname <- covariateFile$Covariate[IDindex]
    randID <- IDname
    # Remove ID variable from covariate vector
    coVars <- coVars[-grep(IDname, coVars)]
  } else {
    randID <- "SampleID"
  }

  # Now create the full model formulas
  # Start with the fixed effects
  fixedStr <- paste0("lIntensity ~ ", paste(c(coVars, timeVars), collapse = " + "))

  # If we have a time category, add interaction terms
  if (tCatIndex > 0) {
    fixedStr <- paste0(
      fixedStr, " + ",
      paste0(c(rep(paste0(tCatName, ":"), length(timeVars))),
        timeVars,
        collapse = " + "
      )
    )
  }

  fixedForm <- as.formula(fixedStr)




  ################# Restructure the data for modeling##############
  # Recap:  The relevant data tables are normI, snMat, and DF
  # These must be melted and merged into a data frame with one
  # observation per row.  Then we need to merge the table with the
  # information in sampleFile

  idVars <- DF[, c("PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex")]
  meltI <- reshape2::melt(data.frame(Scan = 1:nrow(DF), idVars, normI),
    id.vars = c("Scan", "PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex")
  )
  meltSn <- reshape2::melt(data.frame(Scan = 1:nrow(DF), idVars, snMat),
    id.vars = c("Scan", "PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex"),
    value.name = "SNR", variable.name = "Channel"
  )
  # Extract and transform outcomes
  meltSn$lIntensity <- log2(meltI$value)
  meltSn$SampleID <- paste0(meltSn$Plex, "_", meltSn$Channel)
  # This should match The sample ID in the sampleFile.

  covNames <- colnames(sampleFile)[-1]
  if (length(covNames) > 0) {

    # Now map the sampleFile information into emptyCov
    # make sure the sample ID's match
    dataTosampFile <- match(meltSn$SampleID, sampleFile$SampleID)
    undescribed <- which(is.na(dataTosampFile))
    if (length(undescribed) > 0) {
      meltSn <- meltSn[-undescribed, ]
    }

    idMatch <- match(meltSn$SampleID, sampleFile$SampleID)

    sampFileToData <- match(sampleFile$SampleID, meltSn$SampleID)
    if (sum(is.na(sampFileToData)) > 0) {
      stop("Error: Samples listed in the sampleFile are not present in the data")
    }
    # now remove any samples that were not in the sampleFile

    emptyCov <- as.data.frame(matrix(NA, nrow = nrow(meltSn), ncol = length(covNames)))
    colnames(emptyCov) <- covNames

    for (i in 1:length(covNames)) {
      emptyCov[, i] <- sampleFile[idMatch, i + 1]
    }

    readyDf <- data.frame(meltSn, emptyCov)
  } else { # If there were no covariates
    readyDf <- meltSn
  }

  # Anything set to NA should be removed at this point
  naIndex <- which(is.na(readyDf$lIntensity))
  if (length(naIndex) > 0) {
    readyDf <- readyDf[-naIndex, ]
  }


  # In the absence of a bridge channel (not recommended) we
  # will align all of the scan means
  if (bridgeMod == FALSE) {
    # Get scan means
    scanMean <- tapply(readyDf$lIntensity,
      readyDf$Scan,
      FUN = mean
    )
    scanFactors <- scanMean
    readyDf$lIntensity <- readyDf$lIntensity -
      scanFactors[match(readyDf$Scan, names(scanFactors))]
    bridgeDat <- NULL # This gets passed into the modeling function
  } else {
    # Make a boolean for bridge observations
    bridgeIDs <- sampleFile$SampleID[which(sampleFile$Bridge == 1)]
    readyDf$Bridge <- 0
    readyDf$Bridge[which(readyDf$SampleID %in% bridgeIDs)] <- 1
  }

  # Create names of samples to populate the results table
  sampleNames <- levels(factor(readyDf$SampleID))
  readyDf$Plex <- factor(readyDf$Plex)
  # if(bridgeMod == TRUE){
  # Bridge channels don't get estimated
  # sampleNames <- sampleNames[-grep(bridgeChannel, sampleNames)]
  # }#
  # Prepare the weighting column
  readyDf$techVar <- 1 / (scaleSN * readyDf$SNR)



  ################# Now generate the empty results tables################
  # The general strategy is to create a separate table for each
  # fixed reference and then let the rest of the variable references
  # float.  We repeat, iterating through all possible references.
  uProt <- unique(readyDf$Protein.ID)
  uGene <- readyDf$PA.Gene.Symbol[match(uProt, readyDf$Protein.ID)]
  n_prot <- length(uProt)

  n_tables <- 1 + sum(n_levels) + n_cont

  # Loop through the tables


  tableList <- list()

  # Generate the sample estimate table
  resCols <- 3 * length(sampleNames)
  tempTab <- data.frame(
    uGene, uProt, rep(NA, length(uGene)),
    matrix(NA, nrow = length(uProt), ncol = resCols)
  )
  colnames(tempTab) <- c(
    "Gene", "Protein", "modelFit",
    paste0(
      rep(c("lower_", "est_", "upper_"), length(sampleNames)),
      rep(sampleNames, each = 3)
    )
  )
  tableList[[1]] <- tempTab


  # Generate any factor tables
  tabIndex <- 2

  if (n_cat > 0) {
    for (i in 1:length(factorNames)) {
      levelN <- length(levelNames[[i]])
      for (j in 1:levelN) {
        # Define Columns: est, pval, qval for (levels-1)
        resCols <- 3 * (levelN - 1)
        tempTab <- data.frame(
          uGene, uProt, rep(NA, length(uGene)),
          matrix(NA, nrow = length(uProt), ncol = resCols)
        )
        colnames(tempTab) <- c(
          "Gene", "Protein", "modelFit",
          paste0(
            rep(c("Est_", "Pval_", "Qval_"), length(levelN - 1)),
            rep(levelNames[[i]][-j], each = 3)
          )
        )
        tableList[[tabIndex]] <- tempTab
        tabIndex <- tabIndex + 1
      } # End level loop
    } # End factor loop
  } # End "if" factors exist

  # Generate tables for continuous covariates
  if (n_cont > 0) {
    for (i in 1:length(contNames)) {
      # Define Columns: est, pval, qval
      tempTab <- data.frame(
        uGene, uProt, rep(NA, length(uGene)),
        matrix(NA, nrow = length(uProt), ncol = 3)
      )
      colnames(tempTab) <- c(
        "Gene", "Protein", "modelFit",
        paste0(
          c("Est_", "Pval_", "Qval_"),
          contNames[i]
        )
      )
      tableList[[tabIndex]] <- tempTab
      tabIndex <- tabIndex + 1
    } # End continuous covariate loop
  } # End "if" continuous

  # Generate time tables
  # We need a separate list for these tables since they may be
  # populated at different parts of the code, depending on the
  # appropriate parameterization (timeParm)
  timeTables <- list()
  timeTabIndex <- 1
  if (timeLevelN > 0) { # This means that there are time parameters
    nParams <- length(timeVars)

    # First figure out how many we need to make
    if (timeLevelN == 1) {
      # Just one table
      nTests <- 1

      tempTab <- data.frame(
        uGene, uProt, rep(NA, length(uGene)),
        matrix(NA,
          nrow = length(uProt),
          ncol = nParams + nTests * 2
        )
      )

      colStr <- c("Gene", "Protein", "modelFit", timeVars, "Pval_Time", "Qval_Time")

      colnames(tempTab) <- colStr

      timeTables[[timeTabIndex]] <- tempTab
      timeTabIndex <- timeTabIndex + 1
    } else { # There must be a time category
      nTests <- timeLevelN
      # make timeLevelN number of tables
      for (l_ in 1:timeLevelN) {
        tempTab <- data.frame(
          uGene, uProt, rep(NA, length(uGene)),
          matrix(NA,
            nrow = length(uProt),
            ncol = nParams + nTests * 2
          )
        )
        colStr <- c(
          "Gene", "Protein", "modelFit", timeVars,
          paste0(
            rep(c("Pval_Time_", "Qval_Time_"), length(timeLevelN)),
            rep(c("REF", tCatLevels[-l_]), each = 2)
          )
        )

        colnames(tempTab) <- colStr

        timeTables[[timeTabIndex]] <- tempTab
        timeTabIndex <- timeTabIndex + 1
      } # End time level loop
    } # End else statement
  } # End creation of time tables





  ##################### Table Creation Complete######################

  # Now loop through references and fit the full model for each protein
  for (prot in 1:length(uProt)) {
    protDat <- readyDf[which(readyDf$Protein.ID == uProt[prot]), ]

    # refactor scan numbers or convert to a standard intercept
    protDat$Scan <- factor(protDat$Scan)
    if (length(levels(protDat$Scan)) == 1) {
      # If there is only one scan, then it should be a single
      # vector of 1's in the design matrix
      protDat$Scan <- 1
    }
    # Split Data and refactor SampleIDs

    if (bridgeMod == TRUE) {
      # This was set to NULL above and will remain NULL
      # if bridgeMod ==FALSE
      bridgeDat <- protDat[which(protDat$Bridge == 1), ]
      notBridge <- protDat[which(protDat$Bridge == 0), ]
      notBridge$SampleID <- factor(notBridge$SampleID)
    } else {
      bridgeDat <- NULL
      notBridge <- protDat
      notBridge$SampleID <- factor(notBridge$SampleID)
    }


    # First Estimate the protein in each sample
    tempRes <- try(adaptModel(notBridge,
      covType = "None",
      fixedForm, bridgeDat, tParm,
      fullColumns = sampleNames, minRE, randID,
      reducedMod = FALSE, multiLevel = FALSE
    ))
    tableList[[1]][prot, 4:ncol(tableList[[1]])] <- tempRes[[1]]
    tableList[[1]][prot, 3] <- tempRes[[3]]

    # Now loop through possible references
    # Mimic the conditional loops used for table generation
    # Start with factor covariates

    # Reset the table index
    tabIndex <- 2

    # Populate any factor tables
    if (n_cat > 0) {
      factorCols <- which(colnames(protDat) %in% factorNames)

     # Now create a table counting the number of observed
      # samples in each factor cell
      
      #Re-factor to get rid of any levels that won't appear in this model
      for(f in 1:length(factorCols)){
        notBridge[, factorCols[f]] <- factor(notBridge[, factorCols[f]])
      }
      
      countTable <- tapply(notBridge$SampleID, notBridge[, factorCols],
                           FUN = function(x) length(unique(x))
      )
      
      N_PerLevel <- lapply(factorCols, function(x)
        tapply(notBridge$SampleID, notBridge[, x],
               FUN = function(x) length(unique(x))))
      minPerLevel <- min(unlist(lapply(N_PerLevel, min)))
      
      nCells <- sum(!is.na(countTable))
      nFactors <- length(factorCols)
      nLevels <- sum(unlist(lapply(N_PerLevel, length)))
        
        
      # If we have continuous covariates, make sure they have
      # support in every level (> 2 + n_cont)
      if (n_cont > 0 & minPerLevel < n_cont + 2) {
        reducedMod <- TRUE
        # This will first be used to used to remove continous
        # covariates from tempCov
        # Then it will be passed into adapt model to return
        # an indicator of the reduction
        coVars_reduced <- coVars[-which(coVars %in% contNames)]
        p_ <- nLevels - nFactors + length(timeVars)
      } else {
        reducedMod <- FALSE
        coVars_reduced <- coVars
        p_ <- nLevels - nFactors + length(timeVars) + n_cont
      }
      
      # Make sure we have at least minRE "extra" samples to estimate between
      # sample variance 
      if (sum(countTable, na.rm = T) - p_ < minRE) {
        multiLevel <- FALSE
      } else {
        multiLevel <- TRUE
      }


      for (i in 1:length(factorNames)) {
        allLevels <- levelNames[[i]]
        levelN <- length(levelNames[[i]])
        factorIndex <- which(colnames(protDat) == factorNames[i])
        for (j in 1:levelN) {
          # See if the reference is present
          obsLevels <- levels(factor(notBridge[, factorNames[i]]))
          # Create a boolean vector that tells us which contrasts to expect
          # Used in the creation of time related contrasts
          # defines tests that are not the reference, hence the -1
          obsBool <- rep(FALSE, length(allLevels) - 1)
          obsBool[allLevels[-j] %in% obsLevels] <- TRUE


          refIndex <- which(obsLevels == allLevels[j])
          # If the reference is missing, move on
          if (length(refIndex) == 0) {
            tabIndex <- tabIndex + 1
            next
          }


          # refactor target variable
          notBridge[, factorNames[i]] <- factor(notBridge[, factorNames[i]],
            levels = c(obsLevels[refIndex], obsLevels[-refIndex])
          )

          # refactor the other factor variables and adjust model formula
          dropInteract <- FALSE # Initialize boolean for dropping interactions
          tempCov <- coVars_reduced
          for (z in 1:length(factorCols)) {
            if (colnames(notBridge)[factorCols[z]] != factorNames[i]) {
              notBridge[, factorCols[z]] <- factor(notBridge[, factorCols[z]])
            }
            # If there is only one level, change the model string
            if (length(levels(notBridge[, factorCols[z]])) < 2) {
              # Create a reduced covariate vector
              tempCov <- tempCov[-which(tempCov == factorNames[z])]
              # Set a boolean to drop interactions if this is a time category
              if (tCatFactorIndex == z) {
                dropInteract <- TRUE
              }
            } # End if ("one-level")
          } # End adjust factor loop
          # tempCov has been reduced for missing data.
          # Reset the model string
          tempStr <- paste0("lIntensity ~ ", paste(c(tempCov, timeVars), collapse = " + "))

          # If we have a time category
          if (tCatIndex > 0 & dropInteract == FALSE) {
            tempStr <- paste0(
              tempStr, " + ",
              paste0(c(rep(paste0(tCatName, ":"), length(timeVars))),
                timeVars,
                collapse = " + "
              )
            )
          }



          fullColumns <- paste0(factorNames[i], allLevels[-j])
          lhtList <- list() # Initialize empty list
          if (tCatFactorIndex == i) { # If the factor "i" is the time category
            # Set the joint hypothesis tests

            for (tNumber in 1:length(allLevels)) {
              if (tNumber == 1) {
                testStr <- paste0(paste0("X_", timeVars)) # ,
                # rep(" = 0", length(timeVars)))
              } else {
                if (obsBool[tNumber - 1] == FALSE) { # Set the test to NA if level is missing
                  lhtList[[tNumber]] <- NA
                  next
                }
                testStr <- paste0(paste0(
                  "X_",
                  factorNames[i], allLevels[-j][tNumber - 1],
                  ":", timeVars
                ))
                #  rep(" = 0", length(timeVars)))
              }
              lhtList[[tNumber]] <- testStr
            }
            # If we don't want to test time trends across categories, then shrink lhtList
            if (timeDiff == FALSE) {
              lhtList <- lhtList[1]
            }
          }


          # If only one level, move on
          if (length(obsLevels) < 2) {
            tabIndex <- tabIndex + 1
            next
          } # Note this escape needs to occure before turning
          # tempStr into a formula but after lhtList is initialized

          tempForm <- as.formula(tempStr)

          # If there is only 1 level, we move on (unless there is a time trend)
          # Note this was placed after the creation of tempForm so that
          # tempForm will be correct in the continuous covariate section.
          if (length(obsLevels) < 2 & (tCatFactorIndex != i)) {
            tabIndex <- tabIndex + 1
            next
          }

          tempRes <- try(adaptModel(
            protDat = notBridge,
            covType = "Factor", fixedForm = tempForm,
            bridgeDat, tParm, lhtList,
            fullColumns, minRE, randID,
            reducedMod, multiLevel
          ))

          tableList[[tabIndex]][prot, 4:ncol(tableList[[tabIndex]])] <- tempRes[[1]]
          tableList[[tabIndex]][prot, 3] <- tempRes[[3]]
          # Extract time results if they exist
          if (length(lhtList) > 0) {
            timeTables[[j]][prot, 4:ncol(timeTables[[j]])] <- tempRes[[2]]
            timeTables[[j]][prot, 3] <- tempRes[[3]]
          }

          # Increment to the next reference
          tabIndex <- tabIndex + 1
        } # End level loop
      } # End factor loop
    } else { # If there are no categorical variables, the model does not need
      # any adjustments.  We still pass "tempForm".
      tempForm <- fixedForm
    } # End "if" categorical
    

    # Populate tables for continuous covariates
    if (n_cont > 0) {
      # For continuous covariates, only one model fit is necessary (The
      # model that allows references to float across proteins).
      # So first we will fit the model and then loop through the
      # covariates to populate each table

      # If there are any categorical covariates, refactor them
      if (n_cat > 0) {
        for (z in 1:length(factorCols)) {
          notBridge[, factorCols[z]] <- factor(notBridge[, factorCols[z]])
        }
        # If the continuous covariate was removed, put it back
        if (reducedMod == TRUE) {
          tempCov <- c(tempCov, contNames)
          tempStr <- paste0("lIntensity ~ ", paste(c(tempCov), collapse = " + "))
          tempForm <- as.formula(tempStr)
        }
      } else {
        multiLevel <- TRUE # needs to be initialized
      }

      tempRes <- try(adaptModel(
        protDat = notBridge,
        covType = "Continuous", fixedForm = tempForm,
        bridgeDat, tParm, lhtList = NULL,
        fullColumns = contNames, minRE, randID,
        reducedMod = FALSE, multiLevel
      ))

      for (i in 1:length(contNames)) {
        tableList[[tabIndex]][prot, 4:6] <- tempRes[[1]][i, ]
        tableList[[tabIndex]][prot, 3] <- tempRes[[3]]
        tabIndex <- tabIndex + 1
      } # End populating continuous tables
    } # End "if" continuous

    
    # If it hasn't been done yet, populate the time tables
    if (tParm == "Time") {
      #If we are here, then lhtList was not initialized in the factor loop
      lhtList <- list()
      lhtList[[1]] <- paste0(paste0("X_", timeVars))
      
      # If there are any categorical covariates, refactor them
      if (n_cat > 0) {
        for (z in 1:length(factorCols)) {
          notBridge[, factorCols[z]] <- factor(notBridge[, factorCols[z]])
        }
      }

      tempRes <- try(adaptModel(
        protDat = notBridge,
        covType = "Time", fixedForm = tempForm,
        bridgeDat, tParm, lhtList,
        fullColumns = timeVars, minRE, randID,
        reducedMod = FALSE, multiLevel = TRUE
      ))

    #There is only 1 time table
      timeTables[[1]][prot, 4:ncol(timeTables[[1]])] <- tempRes[[2]]
      timeTables[[1]][prot, 3] <- tempRes[[3]]
      
  
  } # End Time tests
    
    
    # Progress message
    print(paste0("Protein ", prot, " Complete"))
  } # End Protein Loop



  # We're almost done.  Calculate FDR adjustments and save the tables
  #Also replace p-values that rounded to zeroe with the minimum observed value

  for (z_ in 1:length(tableList)) {
    tempTable <- tableList[[z_]]
    qIndex <- grep("Qval", colnames(tempTable))
    if (length(qIndex) > 0) {
      for (q in 1:length(qIndex)) {
        #Replace 0's
        indexZero <- which(tempTable[ , qIndex[q] - 1] == 0)
        if(length(indexZero) > 0){
          tempTable[indexZero , qIndex[q] - 1] <- 10^-300
        }
        tempTable[, qIndex[q]] <- p.adjust(tempTable[, qIndex[q] - 1],
          method = "fdr"
        )
      }
    }
    write.csv(tempTable, file = tableNames[z_])
  } # End writing of results tables

  # Repeat for time tables
  if (length(timeTables) > 0) {
    for (z_ in 1:length(timeTables)) {
      tempTable <- timeTables[[z_]]
      qIndex <- grep("Qval", colnames(tempTable))
      if (length(qIndex) > 0) {
        for (q in 1:length(qIndex)) {
            #Replace 0's
            indexZero <- which(tempTable[ , qIndex[q] - 1] == 0)
            if(length(indexZero) > 0){
              tempTable[indexZero , qIndex[q] - 1] <- 10^-300
            }
          tempTable[, qIndex[q]] <- p.adjust(tempTable[, qIndex[q] - 1],
            method = "fdr"
          )
        }
      }
      write.csv(tempTable,
        file = timeTableNames[z_],
        row.names = FALSE
      )
    }
  } # End writing time tables
} # End msTrawl function



# Create a function to perform global adjustments and make
# protein subsets for parallel analysis

#' Prepare protein subsets for parallel analysis
#' @param setSize An integer value denoting the approximate
#'  number of proteins that will be included in each
#'  subset of data.
#'  if setSize == 50 and the number of proteins is 99,
#'      there will be one set of size 99.
#'  If setSize == 50 and the number of proteins is 100,
#'      there will be two sets of size 50.
#'  If set Size == 50 and the number of proteins is 101,
#'      there will be two set, one of size 50, the other of size 51.
#'
#' @inheritParams msTrawl
#' @export
protPrep <- function(DF,
                     setSize = 50,
                     sampleFile = NULL,
                     covariateFile = NULL,
                     scaleSN = 1,
                     lod = 0.01,
                     imputePenalty = 1,
                     minAbove = 3,
                     ssnFilter = 20,
                     outlierCutoff = 3,
                     N_SUM = 3,
                     swapProtein = FALSE,
                     maxPep = 25,
                     colAdjust = 0.5,
                     colRatios = NULL,
                     dropContam = TRUE,
                     dropReverse = TRUE,
                     peptideAnalysis = FALSE,
                     minRE = 5,
                     timeDiff = TRUE) {

  # Step1.  Make sure the covariate file and sample file are consistent
  if (!is.null(covariateFile)) {
    # FORCES ROW ORDERING OF VARIABLES TO MATCH SAMPLE FILE COLUMNS (EXCLUDING BRIDGE AND SAMPLEID)
    covariateFile <- covariateFile[order(match(covariateFile$Covariate, colnames(sampleFile))), ]

    coVector <- covariateFile$Covariate
    sampVector <- colnames(sampleFile)
    nMatches <- sum(coVector %in% sampVector)
    if (nMatches < length(coVector)) {
      stop("Error: At least one of the covariate names does not match across files.")
    }
  }

    # Make sure there are no missing covariates
  if (!is.null(covariateFile)) {
    #ignore bridge channel
    bridgeI <- grep("BRIDGE", toupper(colnames(sampleFile)))
    if(length(bridgeI) > 0){
      nMissing <- sum(is.na(sampleFile[-which(sampleFile[ , bridgeI] == 1) , ]))
    }else{
      nMissing <- sum(is.na(sampleFile))
    }
    
    if (nMissing > 0) {
      stop("Missing values are not allowed in covariates")
    }
  }
  
  # If colAdjust is a vector, add it to DF so it will be included in the filtering steps
  if (length(colAdjust) > 1) {
    DF$colAdjust <- colAdjust
  }

  # Remove any junk
  if (dropReverse == TRUE) {
    revI <- grep("##", DF$Protein.ID)
    if (length(revI) > 0) {
      DF <- DF[-revI, ]
    }
  }
  if (dropContam == TRUE) {
    contamI <- grep("contam", DF$Protein.ID)
    if (length(contamI) > 0) {
      DF <- DF[-contamI, ]
    }
  }
  
  # MSTR-137: TODO: At this point, sometimes DF$Protein.ID is an empty string.
  # These should be be filtered out or otherwise handled at this point or before.

  # Order the rows
  DF <- DF[order(DF$Protein.ID, DF$Plex, DF$Peptide), ]

  # Swap label columns if peptideAnalysis == TRUE
  if (peptideAnalysis) {
    DF$PA.Gene.Symbol <- paste0(DF$PA.Gene.Symbol, "___", DF$Protein.ID)
    DF$Protein.ID <- DF$Peptide
  }

  # Now create the two primary quantification matrices
  snMat <- DF[, grep("\\.Sn", colnames(DF))]
  iMat <- DF[, grep("Adjusted", colnames(DF))]

  # Make sure these have the same number of columns
  if (ncol(snMat) != ncol(iMat)) {
    stop("Error: The number of SNR columns
    does not equal the number of intensity columns. Often this occurs
    because of an unwanted column containing the string \".Sn\" in the
                                     column name.")
  }

  # Apply ssnFilter
  if (!is.null(ssnFilter)) {
    lowI <- which(apply(snMat, 1, sum) < ssnFilter)
    if (length(lowI) > 0) {
      DF <- DF[-lowI, ]
      snMat <- snMat[-lowI, ]
      iMat <- iMat[-lowI, ]
    }
  }

  # Now get rid of any rows with less than 2 non-zero intensities
  emptyRow <- which(apply(iMat, 1, function(x) {
    sum(x > 0, na.rm = T)
  }) < 2)
  if (length(emptyRow > 0)) {
    DF <- DF[-emptyRow, ]
    iMat <- iMat[-emptyRow, ]
    snMat <- snMat[-emptyRow, ]
  }


  ###################### Deal with the LOD#################
  lodRES <- tmtLOD(snMat, iMat, lod, minAbove, scaleSN, imputePenalty)
  # Remove rows with too few entries above LOD
  tooFew <- lodRES[[3]]
  iMat <- lodRES[[1]]
  snMat <- lodRES[[2]]

  if (sum(tooFew) > 0) {
    DF <- DF[-which(tooFew == TRUE), ]
    iMat <- iMat[-which(tooFew == TRUE), ]
    snMat <- snMat[-which(tooFew == TRUE), ]
  }


  ################## Remove outliers from protein groupings##########
  
  # Make sure the relevant ID columns are character vectors
  # so that replacement with peptide names does not result in NA
  DF$Protein.ID <- as.character(DF$Protein.ID)
  DF$Peptide <- as.character(DF$Peptide)
  uPlex <- unique(DF$Plex)
  uProt <- unique(DF$Protein.ID)
  for (j in 1:length(uPlex)) {
    for (i in 1:length(uProt)) {
      protIndex <- which(DF$Protein.ID == uProt[i] &
                           DF$Plex == uPlex[j])
      
      #If the protein is missing in this plex, move on
      if(length(protIndex) == 0){next}
      
      subDat <- iMat[protIndex, , drop = F]
      subSn <- snMat[protIndex, , drop = F]
      subPep <- DF[protIndex, "Peptide"]
      subSSN <- apply(subSn, 1, sum)
      
      
      # Create outlier indicator matrix
      if (!is.null(outlierCutoff)) {
        outlierMat <- findOutliers(subDat, subSn, outlierCutoff, scaleSN)
      } else {
        outlierMat <- matrix(0, nrow = nrow(subDat), ncol = ncol(subDat))
      }
      
      # Figure out which rows have outliers
      outlierRows <- apply(outlierMat, 1, sum, na.rm = T)
      outlierIndex <- protIndex[which(outlierRows > 0)]
      
      # Now change the protein labels of each outlier
      if (length(outlierIndex) > 0) {
        if (swapProtein == TRUE) {
          DF[outlierIndex, "Protein.ID"] <- paste0(DF[outlierIndex, "Protein.ID"], "___", DF[outlierIndex, "Peptide"])
        } else {
          DF[outlierIndex, "Protein.ID"] <- "OUTLIER_REMOVE_AT_ONCE!"
        }
      }
      
      #Check our aggregation criteria.  
      n_after <- length(which(outlierRows == 0))
      if ((n_after < N_SUM & n_after > 1) | (peptideAnalysis == TRUE & n_after > 1)) {
        
        newSnVec <- apply(subSn[which(outlierRows == 0), , drop=FALSE], 2, sum)
        newFluxVec <- apply(subDat[which(outlierRows == 0), , drop=FALSE], 2, sum)
        #Update the peptide name to indicate aggregation
        #Keep the first non-outlier row
        keepIndex <- protIndex[which(outlierRows == 0)][1]
        DF[keepIndex, "Peptide"] <- paste0(DF[protIndex[1], "Peptide"], "_SUM")
        #Update all the quant rows
        DF[keepIndex, grep("\\.Sn", colnames(DF))] <- newSnVec
        DF[keepIndex, grep("Adjusted", colnames(DF))] <- newFluxVec
        snMat[keepIndex, ] <- newSnVec
        iMat[keepIndex, ] <- newFluxVec
        #Flag the rest of the rows for removal
        DF[setdiff(protIndex, keepIndex), "Protein.ID"] <- "OUTLIER_REMOVE_AT_ONCE!"
        next
      }
      #End peptide aggregation step
      
      # Limit the remaining observations to maxPep # of scans
      if (nrow(subDat) - length(outlierIndex) > maxPep) {
        
        # Rank order each scan, smallest to largest, within each peptide
        # Data was ordered by protein, plex, peptide.  We are in a protein
        # and plex loop.  So peptides should be consecutive
        # also the unique() function returns strings in the order observed
        listOranks <- lapply(unique(subPep), function(x) {
          order(subSSN[which(subPep == x)], decreasing = TRUE)
        })
        rankVector <- do.call(c, listOranks)
        
        # Assign outliers an unrealistically high rank
        if (length(outlierIndex) > 0) {
          rankVector[which(outlierRows > 0)] <- 999999
        }
        
        snOrder <- order(rankVector, (-1 * subSSN)) # small rank, large SSN
        
        # Find the sub-indices of the top maxPep signals
        bigIndex <- head(snOrder, n = maxPep)
        # Find the indices of everything but the top maxPep
        smallIndex <- protIndex[setdiff(snOrder, bigIndex)]
        # Change the protein name of the small scans
        DF[smallIndex, "Protein.ID"] <- "TOO_MANY!!!"
      }
    } # End protein loop
  } # End plex loop
  

  outlierIndex <- grep("OUTLIER_REMOVE_AT_ONCE!", DF$Protein.ID)
  tooManyIndex <- grep("TOO_MANY!!!", DF$Protein.ID)
  removeIndex <- union(outlierIndex, tooManyIndex)
  if (length(removeIndex > 0)) {
    DF <- DF[-removeIndex, ]
    iMat <- iMat[-removeIndex, ]
    snMat <- snMat[-removeIndex, ]
  }




  ############## Implement a global column adjustment############

  if (is.null(colAdjust)) {
    normI <- iMat
  } else {
    if (length(colAdjust) == 1) {
      # Intialize boolean
      normBool <- rep(0, nrow(iMat))
      for (p_ in 1:length(uPlex)) {
        # Loop through each plex to find the median
        plexIndex <- which(DF$Plex == uPlex[p_])
        if (length(plexIndex) <= 1) {
          stop("You have a plex with <= 1 observation")
        } # An error occurs, in small subsets of data,
        # where only one observation exists in a plex.  This removes these cases
        # from the normalization

        pepSD <- apply(
          iMat[plexIndex, ],
          1, function(x) sd(log(x), na.rm = TRUE)
        )
        medSD <- quantile(pepSD, probs = colAdjust, na.rm = T)

        sdI <- which(pepSD < medSD)


        normBool[plexIndex[sdI]] <- 1
      } # End plex for loop
    } else {
      normBool <- DF$colAdjust
    }

    if (length(colRatios) < 2) {
      ratios <- rep(1, length(uPlex) * ncol(iMat))
    }
    normed <- geoNorm(iMat, normBool, Plex = DF$Plex, ratios)
    normI <- normed[[1]]
    colnames(normI) <- paste0("norm_", colnames(normI))
    normFacs <- normed[[2]] # Return normalization factors
    # Write the column normalization table
    write.csv(normFacs,
      file = "ColAdjustmentFactors.csv",
      row.names = FALSE
    )
  } # End global column adjustment

  ################# Pre-processing complete. Read the covariate data##########

  # Figure out what we are dealing with.

  # First reduce the number of plexes
  usedPlexes <- intersect(uPlex, substring(sampleFile$SampleID, 1, regexpr("_", sampleFile$SampleID) - 1))

  # Is there a bridge channel?
  nBridges <- sum(sampleFile$Bridge)
  if (nBridges > 0) {
    bridgeMod <- TRUE

    # Quick sanity check on the number of bridge channels
    if (nBridges != length(usedPlexes)) {
      stop("Error: The number of plexes must match the number of bridge
         samples.")
    }
  } else {
    bridgeMod <- FALSE
  }

  # Initialize table names
  tableNames <- "Simple.csv"

  # Figure out the number of categorical covariates
  # Is there a bridge column?
  bridgeCol <- grep("BRIDGE", toupper(colnames(sampleFile)))
  if (length(bridgeCol) > 0) {
    coVars <- colnames(sampleFile)[-c(1, bridgeCol)]
  } else {
    coVars <- colnames(sampleFile)[-1]
  }

  # If sampleFile was NULL or had no covariates then set an intercept model
  if (length(coVars) == 0) {
    coVars <- "1"
  }

  covType <- covariateFile$Type
  covLevel <- covariateFile$Levels
  catIndex <- grep("FACTOR", toupper(covType))
  n_cat <- length(catIndex)


  # Initialize parameter for defining time table workflow
  tParm <- ""

  # Figure out the number of levels of each factor and the names
  if (n_cat > 0) {
    n_levels <- rep(0, length(catIndex)) # initialize vector
    factorNames <- as.character(covariateFile$Covariate[catIndex])

    # Make sure the covariate names match across files
    if (sum(factorNames %in% colnames(sampleFile)) != length(factorNames)) {
      stop("Error: The covariate names in the sampleFile do not match the
         column names in the sampleFile")
    }

    levelNames <- list()
    for (i in 1:n_cat) {
      if (length(bridgeCol) > 0) {
        realIndex <- which(sampleFile[, bridgeCol] == 0)
      } else {
        realIndex <- 1:nrow(sampleFile)
      }

      levelNames[[i]] <- as.character(unique(sampleFile[realIndex, factorNames[i]]))
      n_levels[i] <- length(levelNames[[i]])

      # Update table names
      tableNames <- c(tableNames, paste0(
        "Factor_",
        factorNames[i], "_",
        levelNames[[i]], ".csv"
      ))
    } # End level loop
  } else { # End "if" factors are present
    n_levels <- 0
    factorNames <- ""
    levelNames <- NULL
  }

  # Figure out the number of continuous covariates and center them
  contIndex <- grep("CONTINUOUS", toupper(covariateFile$Type))
  n_cont <- length(contIndex)
  if (n_cont > 0) {
    contNames <- as.character(covariateFile$Covariate[contIndex])
    tableNames <- c(tableNames, paste0("Continuous_", contNames, ".csv"))
    for (i in 1:n_cont) {
      sampleFile[, contNames[i]] <- sampleFile[, contNames[i]] -
        mean(sampleFile[, contNames[i]], na.rm = T)
    }
  } else { # End "if" continuous
    contNames <- ""
  }
  # Figure out time related parameters
  timeIndex <- grep("TIME", toupper(covariateFile$Type))
  if (length(timeIndex) > 0) {
    # First rename the time variable for use in formulas
    colnames(sampleFile)[grep(
      covariateFile$Covariate[timeIndex],
      colnames(sampleFile)
    )] <- "Time"
    covariateFile$Covariate[timeIndex] <- "Time"

    # Make sure this variable is continuous
    if (!is.numeric(sampleFile$Time)) {
      stop("Error:  Your time variable is non-numeric")
    }

    # Remove Time variable from the covariate vector.
    coVars <- coVars[-timeIndex]
    # Start building a separate vector for time variables
    timeVars <- NULL

    # Now set the degree and circadian parameters
    timeDegree <- covariateFile$TimeDegree[timeIndex]
    if (timeDegree >= 1) {
      timeVars <- c(timeVars, "Time")
    }
    if (timeDegree >= 2) {
      timeVars <- c(timeVars, "Time2")
      sampleFile$Time2 <- sampleFile$Time^2
    }
    if (timeDegree >= 3) {
      timeVars <- c(timeVars, "Time3")
      sampleFile$Time3 <- sampleFile$Time^3
    }
    # Look for a factor that separates time trends
    # Find position of time category in covariate file
    tCatIndex <- which(covariateFile$TimeCategory == 1)
    # Find position of time category in factor list
    tCatFactorIndex <- match(covariateFile$Covariate[tCatIndex], factorNames)
    if (length(tCatIndex) > 1) {
      stop("Only one time category is allowed")
    }
    
    if(length(tCatIndex) == 0){
      #Force a value greater than zero to always indicate the presence of a time category
      tCatIndex <- 0
      tCatFactorIndex <- 0
    } 
    
    if (tCatIndex > 0) {
      tCatName <- covariateFile$Covariate[tCatIndex]
      timeLevelN <- n_levels[tCatFactorIndex]
      tCatLevels <- levelNames[[tCatFactorIndex]]
      timeTableNames <- paste0("Time_", tCatLevels, ".csv")
      # There are three possible relevant parameterizations for time
      tParm <- "Category"
    } else { # This means there is no time category but there are time covariates
      timeLevelN <- 1
      timeTableNames <- "Time.csv"
      if (n_cont > 0) {
        tParm <- "Continuous"
      } else {
        tParm <- "Time"
      }
    tCatName <- NULL
    tCatLevels <- NULL
    }
    circadian <- covariateFile$Circadian[timeIndex]
  } else { # No time related covariates
    timeDegree <- 0
    timeVars <- NULL
    timeLevelN <- 0
    tCatIndex <- 0
    tCatFactorIndex <- 0
    circadian <- 0
    tParm <- ""
    timeTableNames <- NULL
    tCatName <- NULL
    tCatLevels <- NULL
  }

  # Create circadian covariates
  if (circadian != 0) {
    sampleFile$Sine <- sin((2 * pi / 24) * sampleFile$Time)
    sampleFile$Cosine <- cos((2 * pi / 24) * sampleFile$Time)
    timeVars <- c(timeVars, "Sine", "Cosine")
  }

  # Is this a longitudinal model?
  IDindex <- grep("ID", toupper(covariateFile$Type))
  if (length(IDindex) > 0) {
    IDname <- covariateFile$Covariate[IDindex]
    randID <- IDname
    # Remove ID variable from covariate vector
    coVars <- coVars[-grep(IDname, coVars)]
  } else {
    randID <- "SampleID"
  }

  # Now create the full model formulas
  # Start with the fixed effects
  fixedStr <- paste0("lIntensity ~ ", paste(c(coVars, timeVars), collapse = " + "))

  # If we have a time category, add interaction terms
  if (tCatIndex > 0) {
    fixedStr <- paste0(
      fixedStr, " + ",
      paste0(c(rep(paste0(tCatName, ":"), length(timeVars))),
        timeVars,
        collapse = " + "
      )
    )
  }

  fixedForm <- as.formula(fixedStr)




  ################# Restructure the data for modeling##############
  # Recap:  The relevant data tables are normI, snMat, and DF
  # These must be melted and merged into a data frame with one
  # observation per row.  Then we need to merge the table with the
  # information in sampleFile

  idVars <- DF[, c("PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex")]
  meltI <- reshape2::melt(data.frame(Scan = 1:nrow(DF), idVars, normI),
    id.vars = c("Scan", "PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex")
  )

  meltSn <- reshape2::melt(data.frame(Scan = 1:nrow(DF), idVars, snMat),
    id.vars = c("Scan", "PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex"),
    value.name = "SNR", variable.name = "Channel"
  )
  # Extract and transform outcomes
  meltSn$lIntensity <- log2(meltI$value)
  meltSn$SampleID <- paste0(meltSn$Plex, "_", meltSn$Channel)
  # This should match The sample ID in the sampleFile.

  covNames <- colnames(sampleFile)[-1]
  if (length(covNames) > 0) {



    # Now map the sampleFile information into emptyCov
    # make sure the sample ID's match
    dataTosampFile <- match(meltSn$SampleID, sampleFile$SampleID)
    undescribed <- which(is.na(dataTosampFile))
    if (length(undescribed) > 0) {
      meltSn <- meltSn[-undescribed, ]
    }

    idMatch <- match(meltSn$SampleID, sampleFile$SampleID)

    sampFileToData <- match(sampleFile$SampleID, meltSn$SampleID)
    if (sum(is.na(sampFileToData)) > 0) {
      stop("Error: Samples listed in the sampleFile are not present in the data")
    }
    # now remove any samples that were not in the sampleFile

    emptyCov <- as.data.frame(matrix(NA, nrow = nrow(meltSn), ncol = length(covNames)))
    colnames(emptyCov) <- covNames

    for (i in 1:length(covNames)) {
      emptyCov[, i] <- sampleFile[idMatch, i + 1]
    }

    readyDf <- data.frame(meltSn, emptyCov)
  } else { # If there were no covariates
    readyDf <- meltSn
  }

  # Anything set to NA should be removed at this point
  naIndex <- which(is.na(readyDf$lIntensity))
  if (length(naIndex) > 0) {
    readyDf <- readyDf[-naIndex, ]
  }


  # In the absence of a bridge channel (not recommended) we
  # will align all of the scan means
  if (bridgeMod == FALSE) {
    # Get scan means
    scanMean <- tapply(readyDf$lIntensity,
      readyDf$Scan,
      FUN = mean
    )
    scanFactors <- scanMean
    readyDf$lIntensity <- readyDf$lIntensity -
      scanFactors[match(readyDf$Scan, names(scanFactors))]
    bridgeDat <- NULL # This gets passed into the modeling function
  } else {
    # Make a boolean for bridge observations
    bridgeIDs <- sampleFile$SampleID[which(sampleFile$Bridge == 1)]
    readyDf$Bridge <- 0
    readyDf$Bridge[which(readyDf$SampleID %in% bridgeIDs)] <- 1
  }

  # Create names of samples to populate the results table
  sampleNames <- levels(factor(readyDf$SampleID))
  readyDf$Plex <- factor(readyDf$Plex)
  # if(bridgeMod == TRUE){
  # Bridge channels don't get estimated
  # sampleNames <- sampleNames[-grep(bridgeChannel, sampleNames)]
  # }
  # Prepare the weighting column
  readyDf$techVar <- 1 / (scaleSN * readyDf$SNR)



  # Now split up readyDf to create subsets for parallel analysis
  uProt <- unique(readyDf$Protein.ID)
  uGene <- readyDf$PA.Gene.Symbol[match(uProt, readyDf$Protein.ID)]
  n_prot <- length(uProt)

  # Create set label for each protein
  
  # # gives quotient operator, e.g. 6 %/% 2 = 3, 10 %/% 3 = 3
  nSets <- n_prot %/% setSize 
  setLabel <- sample(1:nSets, n_prot, replace = TRUE)
  readyDf$SetLabel <- setLabel[match(readyDf$Protein.ID, uProt)]
  levelN <- length(levelNames)
  dfList <- split(readyDf, readyDf$SetLabel)

  neededParams <- list()
  for (s_ in 1:length(dfList)) {
    tempName <- paste0("Subset_", s_, ".rda")
    tempNameCSV <- paste0("normalized___", s_, ".csv")
    readyDf <- dfList[[s_]]

    write.csv(readyDf, file = tempNameCSV, row.names = FALSE)

    # TODO: MSTR-136: Switch from save() and load() to saveRDS() and readRDS()
    save(readyDf, bridgeMod, coVars, tableNames, timeTableNames, tCatName,
      n_levels, n_cont, sampleNames, n_cat, levelNames, factorNames,
      levelN, contNames, timeLevelN, timeVars, tCatLevels, fixedStr,
      tParm, minRE, randID, tCatFactorIndex, tCatIndex, timeDiff,
      file = tempName, compress = "xz"
    )
  }
} # End protPrep

#' Function for modeling rdata created by the protPrep function.
#' The function is intended to be used as part of a parallel
#' processing procedure.
#' @param rdaSubset An rda datafile created by protPrep.  This is a
#'  subset of the data and all the relevant parameters needed to
#'  fit models and generate results.
#' @export
miniTrawl <- function(rdaSubset) {
  load(rdaSubset)
  suffix <- gsub("Subset_", "", gsub(".rda", "", rdaSubset))

  fixedForm <- as.formula(fixedStr)

  ################# Now generate the empty results tables################
  # The general strategy is to create a separate table for each
  # fixed reference and then let the rest of the variable references
  # float.  We repeat, iterating through all possible references.
  uProt <- unique(readyDf$Protein.ID)
  uGene <- readyDf$PA.Gene.Symbol[match(uProt, readyDf$Protein.ID)]
  n_prot <- length(uProt)

  n_tables <- 1 + sum(n_levels) + n_cont

  # Loop through the tables


  tableList <- list()

  # Generate the sample estimate table
  resCols <- 3 * length(sampleNames)
  tempTab <- data.frame(
    uGene, uProt, rep(NA, length(uGene)),
    matrix(NA, nrow = length(uProt), ncol = resCols)
  )
  colnames(tempTab) <- c(
    "Gene", "Protein", "modelFit",
    paste0(
      rep(c("lower_", "est_", "upper_"), length(sampleNames)),
      rep(sampleNames, each = 3)
    )
  )
  tableList[[1]] <- tempTab


  # Generate any factor tables
  tabIndex <- 2

  if (n_cat > 0) {
    for (i in 1:length(factorNames)) {
      levelN <- length(levelNames[[i]])
      for (j in 1:levelN) {
        # Define Columns: est, pval, qval for (levels-1)
        resCols <- 3 * (levelN - 1)
        tempTab <- data.frame(
          uGene, uProt, rep(NA, length(uGene)),
          matrix(NA, nrow = length(uProt), ncol = resCols)
        )
        colnames(tempTab) <- c(
          "Gene", "Protein", "modelFit",
          paste0(
            rep(c("Est_", "Pval_", "Qval_"), length(levelN - 1)),
            rep(levelNames[[i]][-j], each = 3)
          )
        )
        tableList[[tabIndex]] <- tempTab
        tabIndex <- tabIndex + 1
      } # End level loop
    } # End factor loop
  } # End "if" factors exist

  # Generate tables for continuous covariates
  if (n_cont > 0) {
    for (i in 1:length(contNames)) {
      # Define Columns: est, pval, qval
      tempTab <- data.frame(
        uGene, uProt, rep(NA, length(uGene)),
        matrix(NA, nrow = length(uProt), ncol = 3)
      )
      colnames(tempTab) <- c(
        "Gene", "Protein", "modelFit",
        paste0(
          c("Est_", "Pval_", "Qval_"),
          contNames[i]
        )
      )
      tableList[[tabIndex]] <- tempTab
      tabIndex <- tabIndex + 1
    } # End continuous covariate loop
  } # End "if" continuous

  # Generate time tables
  # We need a separate list for these tables since they may be
  # populated at different parts of the code, depending on the
  # appropriate parameterization (timeParm)
  timeTables <- list()
  timeTabIndex <- 1
  if (timeLevelN > 0) { # This means that there are time parameters
    nParams <- length(timeVars)

    # First figure out how many we need to make
    if (timeLevelN == 1) {
      # Just one table
      nTests <- 1

      tempTab <- data.frame(
        uGene, uProt, rep(NA, length(uGene)),
        matrix(NA,
          nrow = length(uProt),
          ncol = nParams + nTests * 2
        )
      )

      colStr <- c("Gene", "Protein", "modelFit", timeVars, "Pval_Time", "Qval_Time")

      colnames(tempTab) <- colStr

      timeTables[[timeTabIndex]] <- tempTab
      timeTabIndex <- timeTabIndex + 1
    } else { # There must be a time category
      nTests <- timeLevelN
      # make timeLevelN number of tables
      for (l_ in 1:timeLevelN) {
        tempTab <- data.frame(
          uGene, uProt, rep(NA, length(uGene)),
          matrix(NA,
            nrow = length(uProt),
            ncol = nParams + nTests * 2
          )
        )
        colStr <- c(
          "Gene", "Protein", "modelFit", timeVars,
          paste0(
            rep(c("Pval_Time_", "Qval_Time_"), length(timeLevelN)),
            rep(c("REF", tCatLevels[-l_]), each = 2)
          )
        )

        colnames(tempTab) <- colStr

        timeTables[[timeTabIndex]] <- tempTab
        timeTabIndex <- timeTabIndex + 1
      } # End time level loop
    } # End else statement
  } # End creation of time tables





  ##################### Table Creation Complete######################

  # Now loop through references and fit the full model for each protein
  for (prot in 1:length(uProt)) {
    protDat <- readyDf[which(readyDf$Protein.ID == uProt[prot]), ]

    # refactor scan numbers or convert to a standard intercept
    protDat$Scan <- factor(protDat$Scan)
    if (length(levels(protDat$Scan)) == 1) {
      # If there is only one scan, then it should be a single
      # vector of 1's in the design matrix
      protDat$Scan <- 1
    }
    # Split Data and refactor SampleIDs

    if (bridgeMod == TRUE) {
      # This was set to NULL above and will remain NULL
      # if bridgeMod ==FALSE
      bridgeDat <- protDat[which(protDat$Bridge == 1), ]
      notBridge <- protDat[which(protDat$Bridge == 0), ]
      notBridge$SampleID <- factor(notBridge$SampleID)
    } else {
      bridgeDat <- NULL
      notBridge <- protDat
      notBridge$SampleID <- factor(notBridge$SampleID)
    }


    # First Estimate the protein in each sample
    tempRes <- try(adaptModel(notBridge,
      covType = "None",
      fixedForm, bridgeDat, tParm,
      fullColumns = sampleNames, minRE, randID,
      reducedMod = FALSE, multiLevel = FALSE
    ))
    tableList[[1]][prot, 4:ncol(tableList[[1]])] <- tempRes[[1]]
    tableList[[1]][prot, 3] <- tempRes[[3]]

    # Now loop through possible references
    # Mimic the conditional loops used for table generation
    # Start with factor covariates

    # Reset the table index
    tabIndex <- 2

    # Populate any factor tables
    if (n_cat > 0) {
      factorCols <- which(colnames(protDat) %in% factorNames)

        # Now create a table counting the number of observed
      # samples in each factor cell
      
      #Re-factor to get rid of any levels that won't appear in this model
      for(f in 1:length(factorCols)){
        notBridge[, factorCols[f]] <- factor(notBridge[, factorCols[f]])
      }
      
      countTable <- tapply(notBridge$SampleID, notBridge[, factorCols],
                           FUN = function(x) length(unique(x))
      )
      
      N_PerLevel <- lapply(factorCols, function(x)
        tapply(notBridge$SampleID, notBridge[, x],
               FUN = function(x) length(unique(x))))
      minPerLevel <- min(unlist(lapply(N_PerLevel, min)))
      
      nCells <- sum(!is.na(countTable))
      nFactors <- length(factorCols)
      nLevels <- sum(unlist(lapply(N_PerLevel, length)))
        
        
      # If we have continuous covariates, make sure they have
      # support in every level (> 2 + n_cont)
      if (n_cont > 0 & minPerLevel < n_cont + 2) {
        reducedMod <- TRUE
        # This will first be used to used to remove continous
        # covariates from tempCov
        # Then it will be passed into adapt model to return
        # an indicator of the reduction
        coVars_reduced <- coVars[-which(coVars %in% contNames)]
        p_ <- nLevels - nFactors + length(timeVars)
      } else {
        reducedMod <- FALSE
        coVars_reduced <- coVars
        p_ <- nLevels - nFactors + length(timeVars) + n_cont
      }
      
      # Make sure we have at least minRE "extra" samples to estimate between
      # sample variance 
      if (sum(countTable, na.rm = T) - p_ < minRE) {
        multiLevel <- FALSE
      } else {
        multiLevel <- TRUE
      }
      
      #End establish multi-level rules
      
      for (i in 1:length(factorNames)) {
        allLevels <- levelNames[[i]]
        levelN <- length(levelNames[[i]])
        factorIndex <- which(colnames(protDat) == factorNames[i])
        for (j in 1:levelN) {
          # See if the reference is present
          obsLevels <- levels(factor(notBridge[, factorNames[i]]))

          # Create a boolean vector that tells us which contrasts to expect
          # Used in the creation of time related contrasts
          # defines tests that are not the reference, hence the -1
          obsBool <- rep(FALSE, length(allLevels) - 1)
          obsBool[allLevels[-j] %in% obsLevels] <- TRUE

          refIndex <- which(obsLevels == allLevels[j])
          # If the reference is missing, move on
          if (length(refIndex) == 0) {
            tabIndex <- tabIndex + 1
            next
          }



          # refactor target variable
          notBridge[, factorNames[i]] <- factor(notBridge[, factorNames[i]],
            levels = c(obsLevels[refIndex], obsLevels[-refIndex])
          )

          # refactor the other factor variables and adjust model formula
          dropInteract <- FALSE # Initialize boolean for dropping interactions
          tempCov <- coVars_reduced
          for (z in 1:length(factorCols)) {
            if (colnames(notBridge)[factorCols[z]] != factorNames[i]) {
              notBridge[, factorCols[z]] <- factor(notBridge[, factorCols[z]])
            }
            # If there is only one level, change the model string
            if (length(levels(notBridge[, factorCols[z]])) < 2) {
              # Create a reduced covariate vector
              tempCov <- tempCov[-which(tempCov == factorNames[z])]
              # Set a boolean to drop interactions if this is a time category
              if (tCatFactorIndex == z) {
                dropInteract <- TRUE
              }
            } # End if ("one-level")
          } # End adjust factor loop
          # tempCov has been reduced for missing data.
          # Reset the model string
          tempStr <- paste0("lIntensity ~ ", paste(c(tempCov, timeVars), collapse = " + "))

          # If we have a time category
          if (tCatIndex > 0 & dropInteract == FALSE) {
            tempStr <- paste0(
              tempStr, " + ",
              paste0(c(rep(paste0(tCatName, ":"), length(timeVars))),
                timeVars,
                collapse = " + "
              )
            )
          }





          fullColumns <- paste0(factorNames[i], allLevels[-j])
          lhtList <- list() # Initialize empty list
          if (tCatFactorIndex == i) { # If the factor "i" is the time category
            # Set the joint hypothesis tests

            for (tNumber in 1:length(allLevels)) {
              if (tNumber == 1) {
                testStr <- paste0("X_", timeVars) # ),
                # rep(" = 0", length(timeVars)))
              } else {
                if (obsBool[tNumber - 1] == FALSE) { # Set the test to NA if level is missing
                  lhtList[[tNumber]] <- NA
                  next
                }
                testStr <- paste0(
                  "X_",
                  factorNames[i], allLevels[-j][tNumber - 1],
                  ":", timeVars
                ) # ,
                # rep(" = 0", length(timeVars)))
              }
              lhtList[[tNumber]] <- testStr
            }
            # If we don't want to test time trends across category levels, shrink lhtList
            if (timeDiff == FALSE) {
              lhtList <- lhtList[1]
            }
          }

          # Now that the tempStr has been updated,
          # If there is nothing to compare, move on
          if (length(obsLevels) < 2) {
            tabIndex <- tabIndex + 1
            next
          } # Note this escape needs to occure before turning
          # tempStr into a formula but after lhtList is initialized

          tempForm <- as.formula(tempStr)


          # If there is only 1 level, we move on (unless there is a time trend)
          # Note this was placed after the creation of tempForm so that
          # tempForm will be correct in the continuous covariate section.
          if (length(obsLevels) < 2 & (tCatFactorIndex != i)) {
            tabIndex <- tabIndex + 1
            next
          }


          tempRes <- try(adaptModel(
            protDat = notBridge,
            covType = "Factor", fixedForm = tempForm,
            bridgeDat, tParm, lhtList,
            fullColumns, minRE, randID,
            reducedMod, multiLevel
          ))

          tableList[[tabIndex]][prot, 4:ncol(tableList[[tabIndex]])] <- tempRes[[1]]
          tableList[[tabIndex]][prot, 3] <- tempRes[[3]]
          # Extract time results if they exist
          if (length(lhtList) > 0) {
            timeTables[[j]][prot, 4:ncol(timeTables[[j]])] <- tempRes[[2]]
            timeTables[[j]][prot, 3] <- tempRes[[3]]
          }

          # Increment to the next reference
          tabIndex <- tabIndex + 1
        } # End level loop
      } # End factor loop
    } else { # If there are no categorical variables, the model does not need
      # any adjustments.  We still pass "tempForm".
      tempForm <- fixedForm
    } # End "if" categorical


    # Populate tables for continuous covariates
    if (n_cont > 0) {
      # For continuous covariates, only one model fit is necessary (The
      # model that allows references to float across proteins).
      # So first we will fit the model and then loop through the
      # covariates to populate each table

      # If there are any categorical covariates, refactor them
      if (n_cat > 0) {
        for (z in 1:length(factorCols)) {
          notBridge[, factorCols[z]] <- factor(notBridge[, factorCols[z]])
        }
        # If the continuous covariate was removed, put it back
        if (reducedMod == TRUE) {
          tempCov <- c(tempCov, contNames)
          tempStr <- paste0("lIntensity ~ ", paste(c(tempCov), collapse = " + "))
          tempForm <- as.formula(tempStr)
        }
      } else {
        multiLevel <- TRUE # Needs to be initialized if no factors are present
      }

      tempRes <- try(adaptModel(
        protDat = notBridge,
        covType = "Continuous", fixedForm = tempForm,
        bridgeDat, tParm, lhtList = NULL,
        fullColumns = contNames, minRE, randID,
        reducedMod = FALSE, multiLevel
      ))
      if(is.null(tempRes)){next}

      for (i in 1:length(contNames)) {
        tableList[[tabIndex]][prot, 4:6] <- tempRes[[1]][i, ]
        tableList[[tabIndex]][prot, 3] <- tempRes[[3]]
        tabIndex <- tabIndex + 1
      } # End populating continuous tables
    } # End "if" continuous


    # If it hasn't been done yet, populate the time tables
    if (tParm == "Time") {
      #If we are here, then lhtList was not initialized in the factor loop
      lhtList <- list()
      lhtList[[1]] <- paste0(paste0("X_", timeVars))
      
      # If there are any categorical covariates, refactor them
      if (n_cat > 0) {
        for (z in 1:length(factorCols)) {
          notBridge[, factorCols[z]] <- factor(notBridge[, factorCols[z]])
        }
      }

      tempRes <- try(adaptModel(
        protDat = notBridge,
        covType = "Time", fixedForm = tempForm,
        bridgeDat, tParm, lhtList,
        fullColumns = timeVars, minRE, randID,
        reducedMod = FALSE, multiLevel = TRUE
      ))

    #There is only 1 time table
      timeTables[[1]][prot, 4:ncol(timeTables[[1]])] <- tempRes[[2]]
      timeTables[[1]][prot, 3] <- tempRes[[3]]
      
  
  } # End Time tests
    
    
    # Progress message
    print(paste0("Protein ", prot, " Complete"))
  } # End Protein Loop



  # We're almost done.  Save the tables

  for (z_ in 1:length(tableList)) {
    tempTable <- tableList[[z_]]
    qIndex <- grep("Qval", colnames(tempTable))
    if (length(qIndex) > 0) {
      for (q in 1:length(qIndex)) {
        #Replace 0's
        indexZero <- which(tempTable[ , qIndex[q] - 1] == 0)
        if(length(indexZero) > 0){
          tempTable[indexZero , qIndex[q] - 1] <- 10^-300
        }
        #No FDR adjustment for subsets
        #tempTable[, qIndex[q]] <- p.adjust(tempTable[, qIndex[q] - 1],
         #                                  method = "fdr"
        #)
      }
    }
    
    # TODO MSTR-138: This will fail if the input parameter rdaSubset
    # refers to an absolute path, or a relative path different from the current
    # working directory.
    filename <- gsub(".csv", paste0("___", suffix, ".csv"), tableNames[z_])
    write.csv(tempTable, file = filename)
  } # End writing of results tables

  # Repeat for time tables
  if (length(timeTables) > 0) {
    for (z_ in 1:length(timeTables)) {
      tempTable <- timeTables[[z_]]
      qIndex <- grep("Qval", colnames(tempTable))
      if (length(qIndex) > 0) {
        for (q in 1:length(qIndex)) {
          #Replace 0's
          indexZero <- which(tempTable[ , qIndex[q] - 1] == 0)
          if(length(indexZero) > 0){
            tempTable[indexZero , qIndex[q] - 1] <- 10^-300
          }
          #No FDR adjustment for subsets
          # tempTable[, qIndex[q]] <- p.adjust(tempTable[, qIndex[q] - 1],
          #                                    method = "fdr"
          # )
        }
      }
      filename <- gsub(".csv", paste0("___", suffix, ".csv"), timeTableNames[z_])
      write.csv(tempTable,
        file = filename,
        row.names = FALSE
      )
    }
  } # End writing time tables
} # End miniTrawl function





