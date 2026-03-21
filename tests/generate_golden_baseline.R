#!/usr/bin/env Rscript
#
# generate_golden_baseline.R
#
# Runs the full msTrawler pipeline on tutorial data with a fixed seed,
# producing deterministic "golden" output CSVs that serve as the reference
# for E2E regression tests.
#
# Usage (inside Docker):
#   Rscript /workspace/tests/generate_golden_baseline.R
#
# Outputs written to /workspace/tests/golden/

library(msTrawler)

cat("=== msTrawler E2E Golden Baseline Generator ===\n")
cat("Package version:", as.character(packageVersion("msTrawler")), "\n")
cat("R version:", R.version.string, "\n\n")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SEED <- 42
OUTPUT_DIR <- "/workspace/tests/golden"
DATA_DIR <- "/msTrawler/data/tutorial_example"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Test 1: convertPdFile — file conversion
# ---------------------------------------------------------------------------
cat("--- Test 1: convertPdFile ---\n")
set.seed(SEED)

df <- convertPdFile(
  pd_psm_file = file.path(DATA_DIR, "pd_example_export_PSMs.txt"),
  pd_protein_file = file.path(DATA_DIR, "pd_example_export_Proteins.txt")
)

write.csv(df, file = file.path(OUTPUT_DIR, "converted_df.csv"), row.names = FALSE)
cat("  Rows:", nrow(df), " Cols:", ncol(df), "\n")
cat("  Saved: converted_df.csv\n\n")

# ---------------------------------------------------------------------------
# Test 2: tmtLOD — LOD imputation (isolated)
# ---------------------------------------------------------------------------
cat("--- Test 2: tmtLOD ---\n")
set.seed(SEED)

snMat <- df[, grep("\\.Sn", colnames(df))]
iMat <- df[, grep("Adjusted", colnames(df))]

lodResult <- tmtLOD(snMat, iMat, lod = 0.01, minAbove = 4, scaleSN = 1, imputePenalty = 1)

write.csv(lodResult[[1]], file = file.path(OUTPUT_DIR, "lod_intensities.csv"), row.names = FALSE)
write.csv(lodResult[[2]], file = file.path(OUTPUT_DIR, "lod_snr.csv"), row.names = FALSE)
write.csv(data.frame(tooFew = lodResult[[3]]), file = file.path(OUTPUT_DIR, "lod_toofew.csv"), row.names = FALSE)
cat("  Saved: lod_intensities.csv, lod_snr.csv, lod_toofew.csv\n\n")

# ---------------------------------------------------------------------------
# Test 3: geoNorm — normalization (isolated)
# ---------------------------------------------------------------------------
cat("--- Test 3: geoNorm ---\n")
set.seed(SEED)

# Use post-LOD data, remove tooFew rows
tooFew <- lodResult[[3]]
cleanI <- lodResult[[1]]
if (sum(tooFew) > 0) {
  cleanI <- cleanI[-which(tooFew == TRUE), ]
}
cleanPlex <- df$Plex
if (sum(tooFew) > 0) {
  cleanPlex <- cleanPlex[-which(tooFew == TRUE)]
}

# Build normalization index (all rows for simplicity)
normBool <- rep(1, nrow(cleanI))
ratios <- rep(1, length(unique(cleanPlex)) * ncol(cleanI))

normResult <- geoNorm(cleanI, normBool, Plex = cleanPlex, ratios)
write.csv(normResult[[1]], file = file.path(OUTPUT_DIR, "geonorm_intensities.csv"), row.names = FALSE)
write.csv(normResult[[2]], file = file.path(OUTPUT_DIR, "geonorm_factors.csv"), row.names = FALSE)
cat("  Saved: geonorm_intensities.csv, geonorm_factors.csv\n\n")

# ---------------------------------------------------------------------------
# Test 4: findOutliers — outlier detection (isolated, one protein)
# ---------------------------------------------------------------------------
cat("--- Test 4: findOutliers ---\n")
set.seed(SEED)

# Pick first protein, first plex
firstProt <- unique(df$Protein.ID)[1]
firstPlex <- unique(df$Plex)[1]
idx <- which(df$Protein.ID == firstProt & df$Plex == firstPlex)

if (length(idx) >= 2) {
  subI <- iMat[idx, , drop = FALSE]
  subSn <- snMat[idx, , drop = FALSE]
  outlierMat <- findOutliers(subI, subSn, cutoff = 3, scaleSN = 1)
  write.csv(outlierMat, file = file.path(OUTPUT_DIR, "outlier_matrix.csv"), row.names = FALSE)
  cat("  Protein:", firstProt, "Plex:", firstPlex, "Scans:", length(idx), "\n")
  cat("  Saved: outlier_matrix.csv\n\n")
} else {
  cat("  Skipped: too few scans for first protein/plex\n\n")
}

# ---------------------------------------------------------------------------
# Test 5: msTrawl — full pipeline
# ---------------------------------------------------------------------------
cat("--- Test 5: msTrawl (full pipeline) ---\n")
set.seed(SEED)

sampleFile <- read.csv(file.path(DATA_DIR, "sample_file.csv"), stringsAsFactors = FALSE)
covariateFile <- read.csv(file.path(DATA_DIR, "covariate_file.csv"), stringsAsFactors = FALSE)

# Run in a temp dir so output CSVs land there
pipeline_dir <- file.path(OUTPUT_DIR, "pipeline")
dir.create(pipeline_dir, recursive = TRUE, showWarnings = FALSE)

old_wd <- getwd()
setwd(pipeline_dir)

tryCatch({
  msTrawl(
    DF = df,
    sampleFile = sampleFile,
    covariateFile = covariateFile,
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
    dropContam = TRUE,
    dropReverse = TRUE,
    peptideAnalysis = FALSE,
    minRE = 5,
    timeDiff = TRUE
  )
  cat("  Pipeline completed successfully.\n")
}, error = function(e) {
  cat("  Pipeline ERROR:", conditionMessage(e), "\n")
})

setwd(old_wd)

# List all outputs
outputs <- list.files(pipeline_dir, pattern = "\\.csv$")
cat("  Output files:", paste(outputs, collapse = ", "), "\n\n")

# ---------------------------------------------------------------------------
# Test 6: msTrawl — no covariates (simple mode)
# ---------------------------------------------------------------------------
cat("--- Test 6: msTrawl (no covariates) ---\n")
set.seed(SEED)

simple_dir <- file.path(OUTPUT_DIR, "pipeline_simple")
dir.create(simple_dir, recursive = TRUE, showWarnings = FALSE)
setwd(simple_dir)

tryCatch({
  msTrawl(
    DF = df,
    sampleFile = sampleFile,
    covariateFile = NULL,
    scaleSN = 1,
    lod = 0.01,
    imputePenalty = 1,
    minAbove = 4,
    colAdjust = 0.5
  )
  cat("  Simple pipeline completed successfully.\n")
}, error = function(e) {
  cat("  Simple pipeline ERROR:", conditionMessage(e), "\n")
})

setwd(old_wd)
outputs2 <- list.files(simple_dir, pattern = "\\.csv$")
cat("  Output files:", paste(outputs2, collapse = ", "), "\n\n")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
all_golden <- list.files(OUTPUT_DIR, pattern = "\\.csv$", recursive = TRUE)
cat("=== Golden baseline complete ===\n")
cat("Total files:", length(all_golden), "\n")
cat("Files:\n")
for (f in all_golden) {
  cat("  ", f, "\n")
}
