#!/usr/bin/env bash
#
# run_e2e.sh — E2E regression test harness for msTrawler
#
# Usage:
#   ./tests/run_e2e.sh generate    # Generate golden baseline (run once before refactoring)
#   ./tests/run_e2e.sh test        # Run candidate and compare against golden baseline
#   ./tests/run_e2e.sh both        # Generate golden + run test (identity check)
#
# Prerequisites: Docker
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
IMAGE_NAME="mstrawler-e2e"
GOLDEN_DIR="$SCRIPT_DIR/golden"
CANDIDATE_DIR="$SCRIPT_DIR/candidate"

# ---------------------------------------------------------------------------
# Build Docker image
# ---------------------------------------------------------------------------
build_image() {
  echo "==> Building Docker image '$IMAGE_NAME'..."
  docker build -t "$IMAGE_NAME" -f "$SCRIPT_DIR/Dockerfile" "$PROJECT_DIR"
  echo "==> Image built successfully."
}

# ---------------------------------------------------------------------------
# Generate golden baseline
# ---------------------------------------------------------------------------
generate_golden() {
  echo "==> Generating golden baseline..."

  # Clean previous golden outputs
  rm -rf "$GOLDEN_DIR"
  mkdir -p "$GOLDEN_DIR"

  docker run --rm \
    -v "$PROJECT_DIR:/workspace:rw" \
    "$IMAGE_NAME" \
    Rscript /workspace/tests/generate_golden_baseline.R

  echo ""
  echo "==> Golden baseline saved to $GOLDEN_DIR"
  echo "    Files:"
  find "$GOLDEN_DIR" -name "*.csv" | sort | while read -r f; do
    echo "      $(basename "$f") ($(wc -l < "$f") lines)"
  done
}

# ---------------------------------------------------------------------------
# Run candidate (re-runs the same pipeline, producing outputs in candidate/)
# ---------------------------------------------------------------------------
run_candidate() {
  echo "==> Running candidate..."

  rm -rf "$CANDIDATE_DIR"
  mkdir -p "$CANDIDATE_DIR"

  # Run the same generator but output to candidate dir
  docker run --rm \
    -v "$PROJECT_DIR:/workspace:rw" \
    -e MSTRAWLER_OUTPUT_DIR=/workspace/tests/candidate \
    "$IMAGE_NAME" \
    Rscript -e "
      # Re-source the generator but redirect output to candidate dir
      SEED <- 42
      OUTPUT_DIR <- Sys.getenv('MSTRAWLER_OUTPUT_DIR', '/workspace/tests/candidate')
      DATA_DIR <- '/msTrawler/data/tutorial_example'

      library(msTrawler)
      dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

      # --- Test 1: convertPdFile ---
      set.seed(SEED)
      df <- convertPdFile(
        pd_psm_file = file.path(DATA_DIR, 'pd_example_export_PSMs.txt'),
        pd_protein_file = file.path(DATA_DIR, 'pd_example_export_Proteins.txt')
      )
      write.csv(df, file = file.path(OUTPUT_DIR, 'converted_df.csv'), row.names = FALSE)

      # --- Test 2: tmtLOD ---
      set.seed(SEED)
      snMat <- df[, grep('\\\\.Sn', colnames(df))]
      iMat <- df[, grep('Adjusted', colnames(df))]
      lodResult <- tmtLOD(snMat, iMat, lod = 0.01, minAbove = 4, scaleSN = 1, imputePenalty = 1)
      write.csv(lodResult[[1]], file = file.path(OUTPUT_DIR, 'lod_intensities.csv'), row.names = FALSE)
      write.csv(lodResult[[2]], file = file.path(OUTPUT_DIR, 'lod_snr.csv'), row.names = FALSE)
      write.csv(data.frame(tooFew = lodResult[[3]]), file = file.path(OUTPUT_DIR, 'lod_toofew.csv'), row.names = FALSE)

      # --- Test 3: geoNorm ---
      set.seed(SEED)
      tooFew <- lodResult[[3]]
      cleanI <- lodResult[[1]]
      cleanPlex <- df\$Plex
      if (sum(tooFew) > 0) {
        cleanI <- cleanI[-which(tooFew == TRUE), ]
        cleanPlex <- cleanPlex[-which(tooFew == TRUE)]
      }
      normBool <- rep(1, nrow(cleanI))
      ratios <- rep(1, length(unique(cleanPlex)) * ncol(cleanI))
      normResult <- geoNorm(cleanI, normBool, Plex = cleanPlex, ratios)
      write.csv(normResult[[1]], file = file.path(OUTPUT_DIR, 'geonorm_intensities.csv'), row.names = FALSE)
      write.csv(normResult[[2]], file = file.path(OUTPUT_DIR, 'geonorm_factors.csv'), row.names = FALSE)

      # --- Test 4: findOutliers ---
      set.seed(SEED)
      firstProt <- unique(df\$Protein.ID)[1]
      firstPlex <- unique(df\$Plex)[1]
      idx <- which(df\$Protein.ID == firstProt & df\$Plex == firstPlex)
      if (length(idx) >= 2) {
        subI <- iMat[idx, , drop = FALSE]
        subSn <- snMat[idx, , drop = FALSE]
        outlierMat <- findOutliers(subI, subSn, cutoff = 3, scaleSN = 1)
        write.csv(outlierMat, file = file.path(OUTPUT_DIR, 'outlier_matrix.csv'), row.names = FALSE)
      }

      # --- Test 5: msTrawl full pipeline ---
      set.seed(SEED)
      sampleFile <- read.csv(file.path(DATA_DIR, 'sample_file.csv'), stringsAsFactors = FALSE)
      covariateFile <- read.csv(file.path(DATA_DIR, 'covariate_file.csv'), stringsAsFactors = FALSE)
      pipeline_dir <- file.path(OUTPUT_DIR, 'pipeline')
      dir.create(pipeline_dir, recursive = TRUE, showWarnings = FALSE)
      old_wd <- getwd()
      setwd(pipeline_dir)
      tryCatch({
        msTrawl(DF = df, sampleFile = sampleFile, covariateFile = covariateFile,
                scaleSN = 1, lod = 0.01, imputePenalty = 1, minAbove = 4,
                ssnFilter = 20, outlierCutoff = 3, N_SUM = 3, swapProtein = FALSE,
                maxPep = 25, colAdjust = 0.5, dropContam = TRUE, dropReverse = TRUE,
                peptideAnalysis = FALSE, minRE = 5, timeDiff = TRUE)
      }, error = function(e) { cat('Pipeline ERROR:', conditionMessage(e), '\n') })
      setwd(old_wd)

      # --- Test 6: msTrawl simple (no covariates) ---
      set.seed(SEED)
      simple_dir <- file.path(OUTPUT_DIR, 'pipeline_simple')
      dir.create(simple_dir, recursive = TRUE, showWarnings = FALSE)
      setwd(simple_dir)
      tryCatch({
        msTrawl(DF = df, sampleFile = sampleFile, covariateFile = NULL,
                scaleSN = 1, lod = 0.01, imputePenalty = 1, minAbove = 4,
                colAdjust = 0.5)
      }, error = function(e) { cat('Simple pipeline ERROR:', conditionMessage(e), '\n') })
      setwd(old_wd)

      cat('Candidate run complete.\n')
    "

  echo "==> Candidate outputs saved to $CANDIDATE_DIR"
}

# ---------------------------------------------------------------------------
# Compare golden vs candidate
# ---------------------------------------------------------------------------
run_comparison() {
  echo "==> Comparing golden baseline vs candidate..."
  echo ""

  docker run --rm \
    -v "$PROJECT_DIR:/workspace:rw" \
    "$IMAGE_NAME" \
    Rscript /workspace/tests/compare_outputs.R \
      /workspace/tests/golden \
      /workspace/tests/candidate \
      1e-6

  local exit_code=$?
  echo ""
  if [ $exit_code -eq 0 ]; then
    echo "==> E2E TEST PASSED"
  else
    echo "==> E2E TEST FAILED"
  fi
  return $exit_code
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
ACTION="${1:-both}"

case "$ACTION" in
  generate)
    build_image
    generate_golden
    ;;
  test)
    build_image
    run_candidate
    run_comparison
    ;;
  both)
    build_image
    generate_golden
    run_candidate
    run_comparison
    ;;
  *)
    echo "Usage: $0 {generate|test|both}"
    exit 1
    ;;
esac
