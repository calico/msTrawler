# Main file for analyzing isobaric tag proteomics data
# All other code should be called from this function.

# Internal constants for row marking during outlier/overflow filtering
.MARKER_OUTLIER <- "OUTLIER_REMOVE_AT_ONCE!"
.MARKER_TOO_MANY <- "TOO_MANY!!!"

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
    covariateFile <- .validate_covariates(covariateFile, sampleFile)
  }
  
  # Clean input data: filter, order, extract matrices
  cleanRes <- .clean_input_data(DF, colAdjust, dropReverse, dropContam,
                                peptideAnalysis, ssnFilter)
  DF <- cleanRes$DF
  snMat <- cleanRes$snMat
  iMat <- cleanRes$iMat


  # Run LOD imputation, outlier detection/removal, and normalization
  ppRes <- .preprocess_pipeline(DF, snMat, iMat, lod, minAbove, scaleSN,
                                imputePenalty, outlierCutoff, N_SUM,
                                swapProtein, maxPep, peptideAnalysis,
                                colAdjust, colRatios)
  DF <- ppRes$DF
  snMat <- ppRes$snMat
  normI <- ppRes$normI
  uPlex <- ppRes$uPlex

  ################# Pre-processing complete. Set up covariates and restructure ##########

  # Set up covariates, factors, time variables, and model formula
  covSetup <- .setup_covariates(sampleFile, covariateFile, uPlex)
  sampleFile    <- covSetup$sampleFile
  bridgeMod     <- covSetup$bridgeMod
  bridgeCol     <- covSetup$bridgeCol
  tableNames    <- covSetup$tableNames
  coVars        <- covSetup$coVars
  n_cat         <- covSetup$n_cat
  factorNames   <- covSetup$factorNames
  levelNames    <- covSetup$levelNames
  n_levels      <- covSetup$n_levels
  n_cont        <- covSetup$n_cont
  contNames     <- covSetup$contNames
  timeVars      <- covSetup$timeVars
  timeDegree    <- covSetup$timeDegree
  timeLevelN    <- covSetup$timeLevelN
  timeTableNames <- covSetup$timeTableNames
  tCatIndex     <- covSetup$tCatIndex
  tCatFactorIndex <- covSetup$tCatFactorIndex
  tCatName      <- covSetup$tCatName
  tCatLevels    <- covSetup$tCatLevels
  tParm         <- covSetup$tParm
  circadian     <- covSetup$circadian
  randID        <- covSetup$randID
  fixedStr      <- covSetup$fixedStr
  fixedForm     <- covSetup$fixedForm

  # Restructure data for modeling (melt + merge with covariates)
  restructRes <- .restructure_for_modeling(DF, normI, snMat, sampleFile, bridgeMod, scaleSN)
  readyDf     <- restructRes$readyDf
  sampleNames <- restructRes$sampleNames


  # Initialize empty results tables
  tabRes <- .init_result_tables(readyDf, sampleNames, n_cat, factorNames,
                                levelNames, n_levels, n_cont, contNames,
                                timeLevelN, timeVars, tCatLevels)
  tableList   <- tabRes$tableList
  timeTables  <- tabRes$timeTables
  uProt       <- tabRes$uProt
  uGene       <- tabRes$uGene
  n_prot      <- tabRes$n_prot





  ##################### Table Creation Complete######################

  # Run the per-protein modeling loop
  modelRes <- .run_protein_models(readyDf, uProt, tableList, timeTables,
                                  bridgeMod, sampleNames, fixedForm,
                                  tParm, minRE, randID,
                                  n_cat, factorNames, levelNames, n_levels,
                                  n_cont, contNames, coVars,
                                  timeVars, timeDiff,
                                  tCatIndex, tCatFactorIndex, tCatName,
                                  timeLevelN)
  tableList  <- modelRes$tableList
  timeTables <- modelRes$timeTables

  # Calculate FDR adjustments and save the tables
  .write_result_tables(tableList, tableNames, timeTables, timeTableNames,
                       applyFDR = TRUE, suffix = NULL)
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
    covariateFile <- .validate_covariates(covariateFile, sampleFile)
  }
  
  # Clean input data: filter, order, extract matrices
  cleanRes <- .clean_input_data(DF, colAdjust, dropReverse, dropContam,
                                peptideAnalysis, ssnFilter)
  DF <- cleanRes$DF
  snMat <- cleanRes$snMat
  iMat <- cleanRes$iMat


  # Run LOD imputation, outlier detection/removal, and normalization
  ppRes <- .preprocess_pipeline(DF, snMat, iMat, lod, minAbove, scaleSN,
                                imputePenalty, outlierCutoff, N_SUM,
                                swapProtein, maxPep, peptideAnalysis,
                                colAdjust, colRatios)
  DF <- ppRes$DF
  snMat <- ppRes$snMat
  normI <- ppRes$normI
  uPlex <- ppRes$uPlex

  ################# Pre-processing complete. Set up covariates and restructure ##########

  # Set up covariates, factors, time variables, and model formula
  covSetup <- .setup_covariates(sampleFile, covariateFile, uPlex)
  sampleFile    <- covSetup$sampleFile
  bridgeMod     <- covSetup$bridgeMod
  bridgeCol     <- covSetup$bridgeCol
  tableNames    <- covSetup$tableNames
  coVars        <- covSetup$coVars
  n_cat         <- covSetup$n_cat
  factorNames   <- covSetup$factorNames
  levelNames    <- covSetup$levelNames
  n_levels      <- covSetup$n_levels
  n_cont        <- covSetup$n_cont
  contNames     <- covSetup$contNames
  timeVars      <- covSetup$timeVars
  timeDegree    <- covSetup$timeDegree
  timeLevelN    <- covSetup$timeLevelN
  timeTableNames <- covSetup$timeTableNames
  tCatIndex     <- covSetup$tCatIndex
  tCatFactorIndex <- covSetup$tCatFactorIndex
  tCatName      <- covSetup$tCatName
  tCatLevels    <- covSetup$tCatLevels
  tParm         <- covSetup$tParm
  circadian     <- covSetup$circadian
  randID        <- covSetup$randID
  fixedStr      <- covSetup$fixedStr
  fixedForm     <- covSetup$fixedForm

  # Restructure data for modeling (melt + merge with covariates)
  restructRes <- .restructure_for_modeling(DF, normI, snMat, sampleFile, bridgeMod, scaleSN)
  readyDf     <- restructRes$readyDf
  sampleNames <- restructRes$sampleNames


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

  # Initialize empty results tables
  tabRes <- .init_result_tables(readyDf, sampleNames, n_cat, factorNames,
                                levelNames, n_levels, n_cont, contNames,
                                timeLevelN, timeVars, tCatLevels)
  tableList   <- tabRes$tableList
  timeTables  <- tabRes$timeTables
  uProt       <- tabRes$uProt
  uGene       <- tabRes$uGene
  n_prot      <- tabRes$n_prot





  ##################### Table Creation Complete######################

  # Run the per-protein modeling loop
  modelRes <- .run_protein_models(readyDf, uProt, tableList, timeTables,
                                  bridgeMod, sampleNames, fixedForm,
                                  tParm, minRE, randID,
                                  n_cat, factorNames, levelNames, n_levels,
                                  n_cont, contNames, coVars,
                                  timeVars, timeDiff,
                                  tCatIndex, tCatFactorIndex, tCatName,
                                  timeLevelN)
  tableList  <- modelRes$tableList
  timeTables <- modelRes$timeTables

  # Save tables (no FDR for parallel subsets)
  # TODO MSTR-138: This will fail if rdaSubset is an absolute path
  .write_result_tables(tableList, tableNames, timeTables, timeTableNames,
                       applyFDR = FALSE, suffix = suffix)
} # End miniTrawl function





