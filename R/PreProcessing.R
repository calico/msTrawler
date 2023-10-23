# A file containing key functions for data pre-processing




#' Impute LOD
#'
#' @param snMat  A matrix of signal-to-noise ratios with
#'  channels in the columns and scans in each row
#' @param iMat A matrix of intensities with
#'  channels in the columns and scans in each row
#' @param lod A value that determines the limit of reliability.
#'  If the value is between zero and one, then
#'  the threshold will be a proportion of the total scan
#'  intensity.  If lod >= 1, then the threshold will
#'  be set at a SNR = lod for all scans.
#' @inheritParams msTrawl
#' @return A length three list. This first component has
#'  a matrix of post imputation intensities and the second
#'  has the post imputation SNRs.  The final component contains
#'  an indicator variable for flagging scans with too few
#'  observations above the lod (as determined by minAbove).
#' @export
tmtLOD <- function(snMat, iMat, lod = 0.01, minAbove = 3, scaleSN, imputePenalty) {
  # Calculate aggregate measures that will be used to
  # find a threshold
  SSN <- apply(snMat, 1, sum, na.rm = T)
  SSI <- apply(iMat, 1, sum, na.rm = T)

  PImat <- iMat / SSI

  # Figure out which values will be censored
  if (lod < 1) {
    censoredB <- (PImat < lod | is.na(PImat)) # boolean matrix
    # Added an na check keep NAs out of the subsetting
  } else {
    censoredB <- (snMat < lod & !is.na(PImat))
  }
  censoredI <- which(censoredB == TRUE)

  tooFew <- apply(censoredB, 1, function(x) {
    sum(!(x), na.rm = T) < minAbove
  })

  # If lod is relative then thresholds can be found independently for
  # SNRs and intensities.  Otherwise the SN limit can be converted
  # to an intensity by matching rank ordered observations.

  if (lod < 1) {
    Emiss <- SSI * lod / 2
    Emat <- matrix(rep(Emiss, ncol(iMat)), ncol = ncol(iMat))
    sdMiss <- SSN * lod / (2 * imputePenalty)
    sdMat <- matrix(rep(sdMiss, ncol(iMat)), ncol = ncol(iMat))
  } else {
    # Find the rank of the SNR nearest to LOD
    rankedSN <- c(snMat)
    rankedSN <- rankedSN[order(rankedSN)]
    snIndex <- max(which(rankedSN > lod))
    # Find the matching rank intensity
    rankedI <- c(iMat)
    rankedI <- rankedI[order(rankedI)]
    lodI <- rankedI[snIndex]

    # Now calculate the relevant thresholds
    Emiss <- lodI / 2
    Emat <- matrix(rep(Emiss, ncol(iMat)), ncol = ncol(iMat))
    sdMiss <- lod / 2
    sdMat <- matrix(rep(sdMiss, ncol(iMat)), ncol = ncol(iMat))
  }


  # Generate imputed values
  sdVec <- unlist(sdMat)[censoredI]
  eVec <- unlist(Emat)[censoredI]
  impVec <- 2^(rnorm(length(sdVec), log2(eVec), 1 / sqrt(scaleSN * sdVec)))

  newSN <- unlist(snMat)
  newSN[censoredB] <- sdVec
  newSN <- matrix(newSN, ncol = ncol(snMat))

  newI <- unlist(iMat)
  newI[censoredB] <- impVec
  newI <- matrix(newI, ncol = ncol(iMat))

  colnames(newI) <- colnames(iMat)
  colnames(newSN) <- colnames(snMat)
  list(newI, newSN, tooFew)
}




#' Perform global column adjustments
#'
#' @param mat  A matrix of reporter ion intensities that need
#'  to be adjusted
#' @param normIndex A boolean vector with length equal to the
#'  number of rows in "mat". Each "1" indicates that the
#'  corresponding row will be used to calculate the normalization
#'  factors.
#' @param Plex A vector denoting plex membership.  Normalization
#'  factors are computed for every channel*plex combination.
#'  Accordingly, this variable determines how many adjustments
#'  will be made.
#' @param ratios A numeric vector with length equal to the
#'  number of normalization factors.  The ratios in this
#'  vector determine the expected relationship between
#'  all of the columns.  Typically it will be a vector of 1's,
#'  but occasionally experiments will include columns at known
#'  dilutions, e.g. a 2x bridge channel.
#' @return A length two list.  the first component includes a
#'  matrix that has been adjusted so that the geometric mean
#'  of the subset of rows specified by "normIndex" will be
#'  equivalent across all columns.  The second list component
#'  contains a vector with all the log additive values used
#'  to adjust each column.
#' @export
geoNorm <- function(mat, normIndex, Plex, ratios) {
  mat[mat == 0] <- NA
  lMat <- log2(mat)

  uPlex <- unique(Plex)
  chMeans <- sapply(1:length(unique(Plex)), function(x) {
    apply(lMat[which(normIndex == 1 & Plex == uPlex[x]), ], 2, mean, na.rm = T)
  })

  grandMean <- mean(chMeans - log2(ratios), na.rm = T) # NA removed in case a whole channel
  # is missing
  normFactors <- chMeans - log2(ratios) - grandMean
  if (length(uPlex) > 1) {
    colnames(normFactors) <- uPlex
  }

  # Now adjust each channel*plex combination and then exponentiate
  normed <- lapply(1:length(uPlex), function(x) {
    2^(t(t(lMat[which(Plex == uPlex[x]), ]) - normFactors[, x]))
  })

  newMat <- mat
  for (i in 1:length(uPlex)) {
    newMat[which(Plex == uPlex[i]), ] <- normed[[i]]
  }

  list(newMat, normFactors)
}



#' Remove "consensus outliers" from protein groupings
#'
#' @param protMat  A matrix of reporter ion intensities from a
#'  single protein.
#' @inheritParams msTrawl
#' @return A matrix of boolean values indicating which entries
#'  met the definition of an outlier.
#' @export
findOutliers <- function(protMat, snMat, outlierCutoff = 3,
                         scaleSN) {
  
  # manipulate the data for model fitting
  logMat <- log2(protMat)
  logMat <- logMat - apply(logMat, 1, mean)
  lDf <- data.frame(Scan = 1:nrow(logMat), logMat)
  melted <- reshape2::melt(lDf, id.vars = "Scan")

  # minSn <- apply(snMat[ , -refCol], 2, function(x) pmin(x, snMat[ , refCol]))
  snDf <- data.frame(Scan = 1:nrow(logMat), snMat)
  meltSn <- reshape2::melt(snDf, id.vars = "Scan")
  melted$minSN <- meltSn$value

  # Now fit the average profile
  mod <- lm(value ~ variable + 0,
    weights = I(scaleSN * minSN),
    data = melted
  )
  mod2 <- lm(value ~ 1,
                    weights = I(scaleSN * minSN),
                    data = melted
  )

  # Get standardized residuals
  resid <- rstandard(mod)
  fullResid <- rep(NA, nrow(melted))
  fullResid[which(!is.na(melted$value))] <- resid
  residMat <- reshape2::dcast(data.frame(cbind(melted, fullResid)),
    Scan ~ variable,
    value.var = "fullResid"
  )
  outlierMat <- 1 * (abs(residMat[, -1]) > outlierCutoff)
  
  #repeat for model2
  resid <- rstandard(mod2)
  fullResid <- rep(NA, nrow(melted))
  fullResid[which(!is.na(melted$value))] <- resid
  residMat <- reshape2::dcast(data.frame(cbind(melted, fullResid)),
                              Scan ~ variable,
                              value.var = "fullResid"
  )
  outlierMat2 <- 1 * (abs(residMat[, -1]) > outlierCutoff)
  
  outlierMat <- pmax(outlierMat, outlierMat2, na.rm = TRUE)

  # Put the reference column back
  # fullOutlier <- data.frame(outlierMat, REF = 0)
  # colnames(fullOutlier)[ncol(fullOutlier)] <- colnames(protMat)[refCol]
  # fullOutlier <- fullOutlier[ , match(colnames(protMat), colnames(fullOutlier))]

  as.matrix(outlierMat)
}


#' Fit the appropriate model
#'
#' @param protDat  A dataframe with all the data needed to model
#'   a single protein.
#' @param covType A string defining the type of
#'    reference covariate.  Acceptable types are "Factor",
#'   "Continuous", "Time" or "None.  "None" specifies the
#'   basic model for estimating each sample
#' @param bridgeDat The subset of data that comes from bridge
#'  channels
#' @param tCatIndex An integer denoting the position in
#'  covString of a categorical variable that should have
#'  unique time trajectories for each level.
#' @param fullColumns The name of every variable that we
#'  would estimate if there were no missing values
#' @param minRE The number of levels of a random variable required to fit a 
#'  random intercept
#' @export
adaptModel <- function(protDat, covType, fixedForm,
                       bridgeDat, timeParm, lhtList,
                       fullColumns, minRE, randID,
                       reducedMod, multiLevel) {

  # This function returns 3 list components.  List[[1]] has
  # estimates from any non-time related model component
  # List[[2]] has the time related entries
  # list[[3]] provides the type of model fit ("lm" or "lme")

  protDat$SampleID <- factor(protDat$SampleID)
  protDat$Scan <- factor(protDat$Scan)
  protDat$Plex <- factor(protDat$Plex)
  if (length(grep("nestedID", colnames(protDat))) > 0) {
    protDat$nestedID <- factor(protDat$nestedID)
  }
  resList <- list() # initialize results

  bridgeMod <- TRUE
  if (is.null(bridgeDat)) {
    bridgeMod <- FALSE
  } else {
    # Refactor bridge variables
    bridgeDat$Plex <- factor(bridgeDat$Plex)
  }

  if (length(levels(protDat$Scan)) == 1) {
    # If there is only one scan, then it should be a single
    # vector of 1's in the design matrix
    protDat$Scan <- 1
  }

  # First deal with the sample estimates.
  if (covType == "None") {
    # We are estimating relative abundance in each sample

    # We can do this per plex to speed up processing of large datasets
    nPlexes <- length(unique(protDat$Plex))
    intList <- list()
    for (p_ in 1:nPlexes) {
      subProt <- protDat[which(protDat$Plex == unique(protDat$Plex)[p_]), ]

      if (bridgeMod == FALSE) {
        tempMod <- try(lm(lIntensity ~ 0 + SampleID,
          weights = 1 / techVar,
          data = subProt
        ), silent = TRUE)
        if (class(tempMod) == "try-error") {
          next
        }

        tempInt <- suppressWarnings(confint(tempMod))
        ints <- as.matrix(data.frame(tempInt[, 1], tempMod$coefficients, tempInt[, 2]))
        intList[[p_]] <- ints
        resList[[3]] <- "lm"
      } else { # Else Fit the bridge model
        subBridge <- bridgeDat[which(bridgeDat$Plex == unique(protDat$Plex)[p_]), ]

        allDat <- rbind(subProt, subBridge) 
        # Refactor for the subset

        allDat$Scan <- factor(allDat$Scan)
        subProt$SampleID <- factor(subProt$SampleID)
        if (length(levels(allDat$Scan)) == 1) {
          scanX <- matrix(1, nrow = nrow(allDat), ncol = 1)
        } else {
          scanX <- model.matrix(~ 0 + Scan, data = allDat)
        }

        if (length(levels(subProt$SampleID)) > 1) { #Note: subProt does not include the bridge
          baseX <- model.matrix(~ 0 + SampleID, data = subProt)
        } else {
          baseX <- model.matrix(~1, data = subProt)
          colnames(baseX) <- levels(subProt$SampleID)
        }

        basePlus <- rbind(baseX, matrix(0, nrow = nrow(subBridge), ncol = ncol(baseX)))
        X_ <- cbind(scanX, basePlus)
        y_ <- allDat$lIntensity
        w_ <- allDat$techVar

        tempMod <- try(lm(y_ ~ 0 + X_, weights = 1 / w_), silent = TRUE)
        if (class(tempMod) == "try-error") {
          next
        }

        tempInt <- suppressWarnings(confint(tempMod))
        ints <- as.matrix(data.frame(tempInt[, 1], tempMod$coefficients, tempInt[, 2]))
        intList[[p_]] <- ints
        resList[[3]] <- "lm"
      }
    } # End plex loop

    # Now reassemble the ints
    ints <- do.call(rbind, intList)

    resVec <- rep(NA, 3 * length(fullColumns))
    for (j in 1:length(fullColumns)) {
      startCol <- 1 + (j - 1) * 3
      targetRow <- grep(paste0("*", fullColumns[j]), rownames(ints))
      if (length(targetRow) > 0) {
        resVec[startCol:(startCol + 2)] <- ints[targetRow, ]
      }
    }

    resList[[1]] <- resVec
    return(resList)
  } # End basic sample estimation

  # Now fit the full model for the given reference
  # We will figure out what to report from the model later
  if (is.null(bridgeDat)) {
    # If we only have 1 scan or if we are fitting random scans, we need a
    # 1 column intercept for scan
    nBatches <- length(unique(protDat$Plex))
    nScans <- length(levels(protDat$Scan))
    # Gate multi-level model on the number of scans/batch
    if (nScans - nBatches < 3) {
      multiLevel <- FALSE
    }

    X_ <- model.matrix(fixedForm, data = protDat)


    y_ <- protDat$lIntensity
    w_ <- protDat$techVar

    subjectID <- protDat[, as.character(randID)]
    sampleID <- protDat[, "SampleID"]
    # Create an intercept with a name to match the bridge channel indicator
    sampleI <- rep(1, length(sampleID))
  } else { # Prepare the bridge model
    allDat <- suppressWarnings(rbind(protDat, bridgeDat))
    # Factors for bridge channels are forced to be NA.  Warning supressed.

    # If we only have 1 scan or if we are fitting random scans, we need a
    # 1 column intercept for scan
    nBatches <- length(unique(allDat$Plex))
    nScans <- length(levels(allDat$Scan))
    # Gate multi-level model on the number of scans/batch
    # Require >= 3 scans in a sample to estimate the within sample variance
    # Require >= 3 samples in a cell to estimate sample variance
    if (nScans - nBatches < 3) {
      multiLevel <- FALSE
    }

    if (nScans == 1) {
      # Directly create intercept when relevant factors has 1 level
      scanX <- matrix(1, nrow = nrow(allDat), ncol = 1)
    } else {      
        # Bridge intercept for each scan (ideal if data permits)
        scanX <- model.matrix(~ 0 + Scan, data = allDat)
    }

    baseX <- model.matrix(fixedForm, data = protDat)




    basePlus <- rbind(baseX, matrix(0, nrow = nrow(bridgeDat), ncol = ncol(baseX)))
    X_ <- cbind(scanX, basePlus)
    y_ <- allDat$lIntensity
    w_ <- allDat$techVar
    Scan <- factor(allDat$Scan)

    # Give bridge channels some subject ID (will be set to zero later)
    subjectID <- protDat[, as.character(randID)]
    subjectID <- factor(c(subjectID, rep(subjectID[1], nrow(bridgeDat))))
    sampleID <- allDat[, "SampleID"]

    # remove bridge samples from the sample IDs (we do not want a random effect
    # for each bridge.  We only want the bridge to estimate scan effects)
    replaceID <- allDat$SampleID[which(allDat$Bridge == 1)[1]]
    sampleID[which(allDat$Bridge == 1)] <- replaceID
    sampleID <- factor(sampleID) # Refactor with fewer levels
    # Now create an indicator variable so that we know which samples
    # truly are not bridge samples
    sampleI <- -1 * (allDat$Bridge - 1)
  } # End bridge model prep
  
  # Issue MSTR-133: dimension mismatch
  if (nrow(X_) != length(y_)) {
    err_msg <- paste0("msTrawler::adaptModel() (in file PreProcessing.R) has encountered an illegal state, and cannot continue.\n",
                      "\ty_ has a length of ", length(y_), " elements.\n",
                      "\tX_ contains ", nrow(X_), " rows (and ", ncol(X_), " columns).\n",
                      "In order to fit the linear model, the length of y_ must equal the number of rows in X_.\n",
                      "This problem may be due to an issue with the initial model conditions, or it may be due to a programming error.\n",
                      "Please contact the application maintainers for more information.\n")
    stop(err_msg)
  }
  
  # Look for any aliased variables, remove them from X_ and
  # flag inestimable parameters for removal at the testing stage
  aliasTest <- alias(y_ ~ 0 + X_)$Complete
  
  if (!is.null(aliasTest)) {
    badPars <- NULL
    for (a_ in 1:nrow(aliasTest)) {
      badPars <- union(badPars, colnames(aliasTest)[which(aliasTest[a_, ] != "0")])
    }
    # Removing the rows from aliasTest will give us a linearly independent X_
    # However the variables flagged in badPars are nonsense
    X_ <- X_[, -which(paste0("X_", colnames(X_)) %in% rownames(aliasTest))]
    # The model has been reduced.  Update the parameter to denote this fact
    reducedMod <- TRUE
  } else {
    badPars <- ""
  }
  
  #Now check to see if we have enough samples to levels to fit a subjectID
  #If we have many measurements but too few subjects, we will convert to fixed
  #effects for each subject.  If almost all the observations are unique 
  #subjects then we drop subject altogether.  
  if(randID != "SampleID"){
    #Do we have enough levels to estimate subject variance?
    subjectVarOK <- (length(unique(subjectID)) >= minRE)
    #Do we have at least 3 repeat measurements?
    repeatOK <- (length(unique(sampleID)) - length(unique(subjectID)) >= 3)
    
    #if either condition failed then reduce the model by changing randID 
    #and converting random to fixed effects
    if(subjectVarOK * repeatOK == 0){
      randID <- "noRE"  
      reducedMod <- TRUE
      
      if(repeatOK){
        uSubjects <- unique(subjectID)
        for(subN in 1: length(uSubjects)){
          subjectVec <- 1 * (subjectID == uSubjects[subN])
          tempX_ <- cbind(X_, subjectVec)
          aliasTest1 <- alias(y_ ~ 0 + tempX_)$Complete
          if(is.null(aliasTest1)){
            X_ <- tempX_
          }else{
            next
          }
        }#End convert to fixed effects
      }

    }#End if something is wrong with the repeat structure
  }#End if SampleID 

    if (randID == "SampleID") {
      lmerForm <- as.formula("y_ ~ X_ + 0 + (0 + sampleI | sampleID)")
      lmerForm2 <- NULL
    } else {
      lmerStr <- paste0("y_ ~ X_ + 0 + (0 + sampleI | subjectID) + (0 + sampleI | subjectID:sampleID)")
      lmerForm <- as.formula(lmerStr)
      lmerStr2 <- paste0("y_ ~ X_ + 0 + (0 + sampleI | sampleID)")
      lmerForm2 <- as.formula(lmerStr2)
    }
  

  if (multiLevel == TRUE) {
    tempMod <- try(lme4::lmer(lmerForm,
      weights = 1 / w_,
      control = lme4::lmerControl(
        check.nobs.vs.nlev = "stop",
        check.nobs.vs.nRE = "stop"
      )
    ),
    silent = TRUE
    )

    if (class(tempMod) == "try-error") {
      # Try the second lmer model, if we have longitudinal data
      if (!is.null(lmerForm2)) {
        tempMod <- try(lme4::lmer(lmerForm2,
          weights = 1 / w_,
          control = lme4::lmerControl(
            check.nobs.vs.nlev = "stop",
            check.nobs.vs.nRE = "stop"
          )
        ),
        silent = TRUE
        )
        # Report lm results
        modelType <- "lmer2"
        resList[[3]] <- "lmer2"
      } # End option B

      if (class(tempMod) == "try-error") {
        # If both failed then use lm (option C)
        tempMod <- try(lm(y_ ~ 0 + X_, weights = 1 / w_), silent = T)
        # If this didn't work, move on
        if (class(tempMod) == "try-error") {
          resList[[3]] <- "none"
          next
        }
        # Report lm results

        modelType <- "lm"
        if (reducedMod == FALSE) {
          resList[[3]] <- "lm"
        } else {
          resList[[3]] <- "lm reduced"
        }
      } # End option C
    } else { # Make note that option A worked
      modelType <- "lmer"
      if (reducedMod == FALSE) {
        resList[[3]] <- "lmer"
      } else {
        resList[[3]] <- "lmer reduced"
      }
    }
  } else { # End attempt multi-level
    if (randID != "SampleID" & randID != "noRE") {
      lmerForm <- as.formula("y_ ~ X_ + 0 + (1 | subjectID)")
      tempMod <- try(lme4::lmer(lmerForm,
                                weights = 1 / w_,
                                control = lme4::lmerControl(
                                  check.nobs.vs.nlev = "stop",
                                  check.nobs.vs.nRE = "stop"
                                )
      ),
      silent = TRUE
      )
      modelType <- "lmer_1_Level"
      if (reducedMod == FALSE) {
        resList[[3]] <- "lmer_1_Level"
      } else {
        resList[[3]] <- "lmer_1_Level reduced"
      }
    }else{
      tempMod <- try(lm(y_ ~ 0 + X_, weights = 1 / w_), silent = T)
      # Report lm results
      modelType <- "lm"
      if (reducedMod == FALSE) {
        resList[[3]] <- "lm"
      } else {
        resList[[3]] <- "lm reduced"
      }
    }
    
    
    # If this didn't work, move on
    if (class(tempMod) == "try-error") {
      resList[[3]] <- "none"
      return(NULL)
    }

  } # End option C


  
################ Model fitting is done###############

  # Extract results depending on type of reference covariate
  if (covType == "Factor") {
    if (length(grep("lm$", modelType)) > 0 | length(grep("lm ", modelType)) > 0) {
      tTable <- summary(tempMod)$coefficients
      estPars <- tTable[, 1] # Analogous to the fixef() vector
      tTable <- tTable[, -c(2:3)]
      # Add an empty column to later populate with q-values
      tTable <- cbind(tTable, rep(NA, nrow(tTable)))

      # Set inestimable parameters to NA
      if (badPars[1] != "") {
        badI <- match(badPars, rownames(tTable))
        tTable[badI, ] <- NA
      }
    }
    if (length(grep("lmer", modelType)) > 0) {
      # If there was no error, report results from the lme model
      # tTable  <- summary(tempMod)$tTable
      # tTable <- tTable[ , -c(2:4)]
      # #Add an empty column to later populate with q-values
      # tTable <- cbind(tTable, rep(NA, nrow(tTable)))

      # Set tTable for an lmer model
      tTable <- matrix(NA, nrow = length(fullColumns), ncol = 3)
      row.names(tTable) <- fullColumns

      # Get model estimates
      estPars <- lme4::fixef(tempMod)

      for (j in 1:length(fullColumns)) {
        # If this parameter is inestimable move on
        badMatch <- grep(paste0(fullColumns[j], "$"), badPars)
        if (length(badMatch) > 0) {
          next
        }


        targetRow <- grep(paste0(fullColumns[j], "$"), names(estPars))
        if (length(targetRow) == 0) {
          next
        } # If we didnt' get an estimate
        # then move on

        tTable[j, 1] <- estPars[targetRow]
        # Do the Kenward Roger Test
        testMat <- matrix(0, nrow = 1, ncol = length(estPars))
        testMat[, targetRow] <- 1
        krT <- try(pbkrtest::KRmodcomp(tempMod, testMat))
        if (!(class(krT) == "try-error")) {
          if(all(!is.na(krT$test$p.value))){
            if(krT$test$p.value[1]  == 1){
              tTable[j, 2] <- 10^(-300)  #There is a bug in pbkrtest
            }else{
              tTable[j, 2] <- krT$test$p.value[1]  
            }#record p-val
          }#End if not NA
        }#End if not a try error
          
        }#End column loop
      }# End lme extraction
   

    # Populate results table based on observed variable names
    # Note: this should work for both lmer and lm model fits
    resVec <- rep(NA, 3 * length(fullColumns))
    for (j in 1:length(fullColumns)) {
      startCol <- 1 + (j - 1) * 3
      targetRow <- grep(paste0(fullColumns[j], "$"), rownames(tTable))
      if (length(targetRow) > 0) {
        resVec[startCol:(startCol + 2)] <- tTable[targetRow, ]
      }
    }

    resList[[1]] <- resVec

    # Results for regular factors are done.  But we might have time
    # parameters to extract from this model as well
    if (length(lhtList) > 0) {
      n_timeParams <- length(lhtList[[1]])
      timeRes <- rep(NA, n_timeParams + length(lhtList) * 2)

      # Populate point estimates of reference parameters
      timePars <- lhtList[[1]] # substring(lhtList[[1]], 1, regexpr("=", lhtList[[1]]) - 2)
      for (par_ in 1:n_timeParams) {
        targetRow <- grep(paste0(timePars[par_], "$"), names(estPars))
        if (length(targetRow) == 0) {
          next
        }
        timeRes[par_] <- estPars[targetRow]
      }
      # Now do all the hypothesis tests
      for (test_ in 1:length(lhtList)) {
        if (is.na(lhtList[[test_]][1])) {
          # Missing values might have led to dropping some tests
          timeRes[n_timeParams + 2 * test_ - 1] <- NA
        } else {
          testPars <- lhtList[[test_]]
          # make sure these parameters are estimable
          inestimable <- 0
          for (b_ in 1:length(testPars)) {
            badMatch <- grep(paste0(testPars[b_], "$"), badPars)
            if (length(badMatch) > 0) {
              inestimable <- inestimable + 1
            }
          }
          if (inestimable > 0) {
            next
          }

          # Now proceed with the testing
          testMat <- matrix(0, nrow = length(testPars), ncol = length(estPars))
          for (tPar in 1:length(testPars)) {
            parIndex <- grep(paste0(testPars[tPar], "$"), names(estPars))
            if (length(parIndex) > 0) {
              testMat[tPar, parIndex] <- 1
            } else {
              testMat[tPar, 1] <- NA
            }
          }

          if (sum(is.na(testMat))) {
            # If a testing parameter was missing move on
            next
          }

          # Do the test
          if (modelType == "lm") {
            tempTest <- try(car::lht(tempMod,
              hypothesis.matrix = testMat, singular.ok = T,
              test = "F"
            )[2, "Pr(>F)"], silent = T)
          } else {
            tempTest <- try(pbkrtest::KRmodcomp(tempMod, testMat)$test$p.value[1], silent = T)
          }



          if (class(tempTest) != "try-error") {
            timeRes[n_timeParams + 2 * test_ - 1] <- tempTest
          }
        } # End if test is !missing
      } # End lht loop

      resList[[2]] <- timeRes
    }#End time testing withing Factor loop

    # Return a list object.  The first entry is for non-time related
    # Parameters.  The second is for the time table.
    return(resList)
  } # End extract factor parameters

  # Now extract continuous covariate
  if (covType == "Continuous") {
    if  (length(grep("lm$", modelType)) > 0 | length(grep("lm ", modelType)) > 0){
      tTable <- summary(tempMod)$coefficients
      tTable <- tTable[, -c(2:3)]
      # Add an empty column to later populate with q-values
      tTable <- cbind(tTable, rep(NA, nrow(tTable)))
    }



    if (length(grep("lmer", modelType)) > 0) {

      # Set tTable for an lmer model
      tTable <- matrix(NA, nrow = length(fullColumns), ncol = 3)
      row.names(tTable) <- fullColumns

      # Get model estimates
      estPars <- lme4::fixef(tempMod)

      for (j in 1:length(fullColumns)) {
        # If this parameter is inestimable move on
        badMatch <- grep(paste0(fullColumns[j], "$"), badPars)
        if (length(badMatch) > 0) {
          next
        }


        targetRow <- grep(paste0(fullColumns[j], "$"), names(estPars))
        if (length(targetRow) == 0) {
          next
        } # If we didnt' get an estimate
        # then move on

        tTable[j, 1] <- estPars[targetRow]
        # Do the Kenward Roger Test
        testMat <- matrix(0, nrow = 1, ncol = length(estPars))
        testMat[, targetRow] <- 1
        krT <- try(pbkrtest::KRmodcomp(tempMod, testMat))
        if (!(class(krT) == "try-error")) {
          if(!is.na(krT$test$p.value)){
            if(krT$test$p.value[1]  == 1){
              tTable[j, 2] <- 10^(-300)  #There is a bug in pbkrtest
            }else{
              tTable[j, 2] <- krT$test$p.value[1]  
            }#record p-val
          }#End if not NA
        }#End if not a try error
        
        
      }
    } # End lme extraction





    # Populate results table based on observed variable names
    resVec <- matrix(NA, nrow = length(fullColumns), ncol = 3)
    for (j in 1:length(fullColumns)) {
      targetRow <- grep(paste0(fullColumns[j], "$"), rownames(tTable))
      if (length(targetRow) > 0) {
        resVec[j, ] <- tTable[targetRow, ]
      }
    }
    # Note that this returns a matrix that must be read
    # differently than the categorical vector
    resList[[1]] <- resVec
  } # End extract continuous parameter

  # Results from factors and the simple model were returned early.
  # If we are here, then we factored to use all the data for either
  # a time covariate.  Check for time tests.
  if (length(lhtList) > 0) {
    #Extract point estimates from an lm model
    if (length(grep("lm$", modelType)) > 0 | length(grep("lm ", modelType)) > 0) {
      tTable <- summary(tempMod)$coefficients
      estPars <- tTable[, 1] # Analagous to the fixef() vector
      tTable <- tTable[, -c(2:3)]
      # Add an empty column to later populate with q-values
      tTable <- cbind(tTable, rep(NA, nrow(tTable)))

      # Set inestimable parameters to NA
      if (badPars[1] != "") {
        badI <- match(badPars, rownames(tTable))
        tTable[badI, ] <- NA
      }
    }
    #Extract from a lmer model
    if (length(grep("lmer", modelType)) > 0) {

      tTable <- matrix(NA, nrow = length(fullColumns), ncol = 3)
      row.names(tTable) <- fullColumns

      # Get model estimates
      estPars <- lme4::fixef(tempMod)

      for (j in 1:length(fullColumns)) {
        # If this parameter is inestimable move on
        badMatch <- grep(paste0(fullColumns[j], "$"), badPars)
        if (length(badMatch) > 0) {
          next
        }


        targetRow <- grep(paste0(fullColumns[j], "$"), names(estPars))
        if (length(targetRow) == 0) {
          next
        } # If we didnt' get an estimate
        # then move on

        tTable[j, 1] <- estPars[targetRow]
      }
    } # End lme extraction


    # 
    # # Results for regular factors are done.  But we might have time
    # # parameters to extract from this model as well
 
      n_timeParams <- length(lhtList[[1]])
      timeRes <- rep(NA, n_timeParams + length(lhtList) * 2)

      # Populate point estimates of reference parameters
      timePars <- lhtList[[1]] # substring(lhtList[[1]], 1, regexpr("=", lhtList[[1]]) - 2)
      for (par_ in 1:n_timeParams) {
        targetRow <- grep(paste0(timePars[par_], "$"), names(estPars))
        if (length(targetRow) == 0) {
          next
        }
        timeRes[par_] <- estPars[targetRow]
      }
      # Now do all the hypothesis tests
      for (test_ in 1:length(lhtList)) {
        if (is.na(lhtList[[test_]][1])) {
          # Missing values might have led to dropping some tests
          timeRes[n_timeParams + 2 * test_ - 1] <- NA
        } else {
          testPars <- lhtList[[test_]]
          # make sure these parameters are estimable
          inestimable <- 0
          for (b_ in 1:length(testPars)) {
            badMatch <- grep(paste0(testPars[b_], "$"), badPars)
            if (length(badMatch) > 0) {
              inestimable <- inestimable + 1
            }
          }
          if (inestimable > 0) {
            next
          }

          # Now proceed with the testing
          testMat <- matrix(0, nrow = length(testPars), ncol = length(estPars))
          for (tPar in 1:length(testPars)) {
            parIndex <- grep(paste0(testPars[tPar], "$"), names(estPars))
            if (length(parIndex) > 0) {
              testMat[tPar, parIndex] <- 1
            } else {
              testMat[tPar, 1] <- NA
            }
          }

          if (sum(is.na(testMat))) {
            # If a testing parameter was missing move on
            next
          }

          # Do the test
          if (modelType == "lm") {
            tempTest <- try(car::lht(tempMod,
                                     hypothesis.matrix = testMat, singular.ok = T,
                                     test = "F"
            )[2, "Pr(>F)"], silent = T)
          } else {
            tempTest <- try(pbkrtest::KRmodcomp(tempMod, testMat)$test$p.value[1], silent = T)
          }



          if (class(tempTest) != "try-error") {
            timeRes[n_timeParams + 2 * test_ - 1] <- tempTest
          }
        } # End if test is !missing
      } # End lht loop

      resList[[2]] <- timeRes
    
    
    #####End test section######
    
    
  }#End lonely time variable
 

  return(resList)
} # End model fitting function
