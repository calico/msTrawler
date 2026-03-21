# Internal helper functions for msTrawler (not exported)
# These eliminate code duplication between msTrawl(), protPrep(), and miniTrawl().

#' Validate covariate and sample file consistency.
#' Returns the (possibly reordered) covariateFile.
#' @keywords internal
.validate_covariates <- function(covariateFile, sampleFile) {
  # FORCES ROW ORDERING OF VARIABLES TO MATCH SAMPLE FILE COLUMNS (EXCLUDING BRIDGE AND SAMPLEID)
  covariateFile <- covariateFile[order(match(covariateFile$Covariate, colnames(sampleFile))), ]

  coVector <- covariateFile$Covariate
  sampVector <- colnames(sampleFile)
  nMatches <- sum(coVector %in% sampVector)
  if (nMatches < length(coVector)) {
    stop("Error: At least one of the covariate names does not match across files.")
  }

  # Make sure there are no missing covariates (ignore bridge channel)
  bridgeI <- grep("BRIDGE", toupper(colnames(sampleFile)))
  if(length(bridgeI) > 0){
    nMissing <- sum(is.na(sampleFile[-which(sampleFile[ , bridgeI] == 1) , ]))
  }else{
    nMissing <- sum(is.na(sampleFile))
  }

  if (nMissing > 0) {
    stop("Missing values are not allowed in covariates")
  }

  covariateFile
}


#' Clean input data: filter contaminants/reverse, order rows, extract matrices,
#' apply SSN filter, remove empty rows.
#' @return list(DF, snMat, iMat)
#' @keywords internal
.clean_input_data <- function(DF, colAdjust, dropReverse, dropContam,
                              peptideAnalysis, ssnFilter) {
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

  list(DF = DF, snMat = snMat, iMat = iMat)
}


#' Run LOD imputation, outlier detection/removal, and global normalization.
#' @return list(DF, snMat, normI, uPlex)
#' @keywords internal
.preprocess_pipeline <- function(DF, snMat, iMat, lod, minAbove, scaleSN,
                                 imputePenalty, outlierCutoff, N_SUM,
                                 swapProtein, maxPep, peptideAnalysis,
                                 colAdjust, colRatios) {
  ###################### Deal with the LOD #################
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


  ################## Remove outliers from protein groupings ##########

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
          DF[outlierIndex, "Protein.ID"] <- .MARKER_OUTLIER
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
      DF[setdiff(protIndex, keepIndex), "Protein.ID"] <- .MARKER_OUTLIER
      next
      }
      #End peptide aggregation step

      # Limit the remaining observations to maxPep # of scans
      if (nrow(subDat) - length(outlierIndex) > maxPep) {

        # Rank order each scan, smallest to largest, within each peptide
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
        DF[smallIndex, "Protein.ID"] <- .MARKER_TOO_MANY
      }
    } # End protein loop
  } # End plex loop

  outlierIndex <- grep(.MARKER_OUTLIER, DF$Protein.ID)
  tooManyIndex <- grep(.MARKER_TOO_MANY, DF$Protein.ID)
  removeIndex <- union(outlierIndex, tooManyIndex)
  if (length(removeIndex > 0)) {
    DF <- DF[-removeIndex, ]
    iMat <- iMat[-removeIndex, ]
    snMat <- snMat[-removeIndex, ]
  }


  ############## Implement a global column adjustment ############

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
        }

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

  list(DF = DF, snMat = snMat, normI = normI, uPlex = uPlex)
}


#' Set up covariates, factors, time variables, and model formula from sampleFile
#' and covariateFile. Returns a named list of all derived configuration.
#' @return named list with all covariate-related variables
#' @keywords internal
.setup_covariates <- function(sampleFile, covariateFile, uPlex) {

  # First reduce the number of plexes
  usedPlexes <- intersect(uPlex, substring(sampleFile$SampleID, 1, regexpr("_", sampleFile$SampleID) - 1))

  # Is there a bridge channel?
  nBridges <- sum(sampleFile$Bridge)
  if (nBridges > 0) {
    bridgeMod <- TRUE
    if (nBridges != length(usedPlexes)) {
      stop("Error: The number of plexes must match the number of bridge
         samples.")
    }
  } else {
    bridgeMod <- FALSE
  }

  # Initialize table names
  tableNames <- "Simple.csv"

  # Is there a bridge column?
  bridgeCol <- grep("BRIDGE", toupper(colnames(sampleFile)))
  if (length(bridgeCol) > 0) {
    coVars <- colnames(sampleFile)[-c(1, bridgeCol)]
  } else {
    coVars <- colnames(sampleFile)[-1]
  }

  if (length(coVars) == 0) {
    coVars <- "1"
  }

  covType <- covariateFile$Type
  covLevel <- covariateFile$Levels
  catIndex <- grep("FACTOR", toupper(covType))
  n_cat <- length(catIndex)

  tParm <- ""

  if (n_cat > 0) {
    n_levels <- rep(0, length(catIndex))
    factorNames <- as.character(covariateFile$Covariate[catIndex])

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

      tableNames <- c(tableNames, paste0(
        "Factor_",
        factorNames[i], "_",
        levelNames[[i]], ".csv"
      ))
    }
  } else {
    n_levels <- 0
    factorNames <- ""
    levelNames <- NULL
  }

  contIndex <- grep("CONTINUOUS", toupper(covariateFile$Type))
  n_cont <- length(contIndex)
  if (n_cont > 0) {
    contNames <- as.character(covariateFile$Covariate[contIndex])
    tableNames <- c(tableNames, paste0("Continuous_", contNames, ".csv"))
    for (i in 1:n_cont) {
      sampleFile[, contNames[i]] <- sampleFile[, contNames[i]] -
        mean(sampleFile[, contNames[i]], na.rm = T)
    }
  } else {
    contNames <- ""
  }

  timeIndex <- grep("TIME", toupper(covariateFile$Type))
  if (length(timeIndex) > 0) {
    colnames(sampleFile)[grep(
      covariateFile$Covariate[timeIndex],
      colnames(sampleFile)
    )] <- "Time"
    covariateFile$Covariate[timeIndex] <- "Time"

    if (!is.numeric(sampleFile$Time)) {
      stop("Error:  Your time variable is non-numeric")
    }

    coVars <- coVars[-timeIndex]
    timeVars <- NULL

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

    tCatIndex <- which(covariateFile$TimeCategory == 1)
    tCatFactorIndex <- match(covariateFile$Covariate[tCatIndex], factorNames)
    if (length(tCatIndex) > 1) {
      stop("Only one time category is allowed")
    }
    if(length(tCatIndex) == 0){
      tCatIndex <- 0
      tCatFactorIndex <- 0
    }

    if (tCatIndex > 0) {
      tCatName <- covariateFile$Covariate[tCatIndex]
      timeLevelN <- n_levels[tCatFactorIndex]
      tCatLevels <- levelNames[[tCatFactorIndex]]
      timeTableNames <- paste0("Time_", tCatLevels, ".csv")
      tParm <- "Category"
    } else {
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
  } else {
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

  if (circadian != 0) {
    sampleFile$Sine <- sin((2 * pi / 24) * sampleFile$Time)
    sampleFile$Cosine <- cos((2 * pi / 24) * sampleFile$Time)
    timeVars <- c(timeVars, "Sine", "Cosine")
  }

  IDindex <- grep("ID", toupper(covariateFile$Type))
  if (length(IDindex) > 0) {
    IDname <- covariateFile$Covariate[IDindex]
    randID <- IDname
    coVars <- coVars[-grep(IDname, coVars)]
  } else {
    randID <- "SampleID"
  }

  fixedStr <- paste0("lIntensity ~ ", paste(c(coVars, timeVars), collapse = " + "))

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

  list(
    sampleFile = sampleFile,
    bridgeMod = bridgeMod,
    bridgeCol = bridgeCol,
    tableNames = tableNames,
    coVars = coVars,
    n_cat = n_cat,
    factorNames = factorNames,
    levelNames = levelNames,
    n_levels = n_levels,
    n_cont = n_cont,
    contNames = contNames,
    timeVars = timeVars,
    timeDegree = timeDegree,
    timeLevelN = timeLevelN,
    timeTableNames = timeTableNames,
    tCatIndex = tCatIndex,
    tCatFactorIndex = tCatFactorIndex,
    tCatName = tCatName,
    tCatLevels = tCatLevels,
    tParm = tParm,
    circadian = circadian,
    randID = randID,
    fixedStr = fixedStr,
    fixedForm = fixedForm
  )
}


#' Restructure data from wide to long format for modeling.
#' Melts intensity/SNR matrices, merges with sample covariates.
#' @return list(readyDf, sampleNames)
#' @keywords internal
.restructure_for_modeling <- function(DF, normI, snMat, sampleFile, bridgeMod, scaleSN) {

  idVars <- DF[, c("PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex")]
  meltI <- reshape2::melt(data.frame(Scan = 1:nrow(DF), idVars, normI),
    id.vars = c("Scan", "PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex")
  )
  meltSn <- reshape2::melt(data.frame(Scan = 1:nrow(DF), idVars, snMat),
    id.vars = c("Scan", "PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex"),
    value.name = "SNR", variable.name = "Channel"
  )

  meltSn$lIntensity <- log2(meltI$value)
  meltSn$SampleID <- paste0(meltSn$Plex, "_", meltSn$Channel)

  covNames <- colnames(sampleFile)[-1]
  if (length(covNames) > 0) {
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

    emptyCov <- as.data.frame(matrix(NA, nrow = nrow(meltSn), ncol = length(covNames)))
    colnames(emptyCov) <- covNames

    for (i in 1:length(covNames)) {
      emptyCov[, i] <- sampleFile[idMatch, i + 1]
    }

    readyDf <- data.frame(meltSn, emptyCov)
  } else {
    readyDf <- meltSn
  }

  # Remove NA lIntensity
  naIndex <- which(is.na(readyDf$lIntensity))
  if (length(naIndex) > 0) {
    readyDf <- readyDf[-naIndex, ]
  }

  # Bridge or scan-mean normalization
  if (bridgeMod == FALSE) {
    scanMean <- tapply(readyDf$lIntensity,
      readyDf$Scan,
      FUN = mean
    )
    scanFactors <- scanMean
    readyDf$lIntensity <- readyDf$lIntensity -
      scanFactors[match(readyDf$Scan, names(scanFactors))]
    bridgeDat <- NULL
  } else {
    bridgeIDs <- sampleFile$SampleID[which(sampleFile$Bridge == 1)]
    readyDf$Bridge <- 0
    readyDf$Bridge[which(readyDf$SampleID %in% bridgeIDs)] <- 1
  }

  sampleNames <- levels(factor(readyDf$SampleID))
  readyDf$Plex <- factor(readyDf$Plex)
  readyDf$techVar <- 1 / (scaleSN * readyDf$SNR)

  list(readyDf = readyDf, sampleNames = sampleNames)
}
