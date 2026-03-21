#!/usr/bin/env Rscript
#
# compare_outputs.R
#
# E2E regression test: compares "candidate" outputs against golden baseline.
# Returns exit code 0 if all checks pass, 1 if any fail.
#
# Usage:
#   Rscript /workspace/tests/compare_outputs.R <golden_dir> <candidate_dir> [tolerance]
#
# Arguments:
#   golden_dir     Path to the golden baseline directory
#   candidate_dir  Path to the candidate output directory
#   tolerance      Numeric tolerance for floating point comparison (default: 1e-6)
#
# The script checks:
#   1. Same set of output files exist
#   2. Same dimensions (rows x cols) per file
#   3. Same column names
#   4. Numeric columns match within tolerance
#   5. Character columns match exactly

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript compare_outputs.R <golden_dir> <candidate_dir> [tolerance]\n")
  quit(status = 1)
}

golden_dir <- args[1]
candidate_dir <- args[2]
tolerance <- if (length(args) >= 3) as.numeric(args[3]) else 1e-6

cat("=== msTrawler E2E Regression Test ===\n")
cat("Golden dir:   ", golden_dir, "\n")
cat("Candidate dir:", candidate_dir, "\n")
cat("Tolerance:    ", tolerance, "\n\n")

# ---------------------------------------------------------------------------
# Collect files
# ---------------------------------------------------------------------------
golden_files <- sort(list.files(golden_dir, pattern = "\\.csv$", recursive = TRUE))
candidate_files <- sort(list.files(candidate_dir, pattern = "\\.csv$", recursive = TRUE))

n_pass <- 0
n_fail <- 0
failures <- character(0)

fail <- function(msg) {
  n_fail <<- n_fail + 1
  failures <<- c(failures, msg)
  cat("  FAIL:", msg, "\n")
}

pass <- function(msg) {
  n_pass <<- n_pass + 1
  cat("  PASS:", msg, "\n")
}

# ---------------------------------------------------------------------------
# Check 1: same file set
# ---------------------------------------------------------------------------
cat("--- Check: file set ---\n")
missing_in_candidate <- setdiff(golden_files, candidate_files)
extra_in_candidate <- setdiff(candidate_files, golden_files)

if (length(missing_in_candidate) > 0) {
  fail(paste("Missing files in candidate:", paste(missing_in_candidate, collapse = ", ")))
}
if (length(extra_in_candidate) > 0) {
  fail(paste("Extra files in candidate:", paste(extra_in_candidate, collapse = ", ")))
}
if (length(missing_in_candidate) == 0 && length(extra_in_candidate) == 0) {
  pass(paste("File sets match (", length(golden_files), "files )"))
}

# ---------------------------------------------------------------------------
# Check 2-5: per-file comparison
# ---------------------------------------------------------------------------
common_files <- intersect(golden_files, candidate_files)

for (f in common_files) {
  cat("\n--- Comparing:", f, "---\n")

  g <- read.csv(file.path(golden_dir, f), stringsAsFactors = FALSE, check.names = FALSE)
  c <- read.csv(file.path(candidate_dir, f), stringsAsFactors = FALSE, check.names = FALSE)

  # Dimensions
  if (nrow(g) != nrow(c) || ncol(g) != ncol(c)) {
    fail(paste0(f, ": dimension mismatch. Golden ", nrow(g), "x", ncol(g),
                " vs Candidate ", nrow(c), "x", ncol(c)))
    next
  }
  pass(paste0(f, ": dimensions match (", nrow(g), "x", ncol(g), ")"))

  # Column names
  if (!identical(colnames(g), colnames(c))) {
    fail(paste0(f, ": column name mismatch"))
    cat("    Golden:", paste(head(colnames(g), 10), collapse = ", "), "\n")
    cat("    Candidate:", paste(head(colnames(c), 10), collapse = ", "), "\n")
    next
  }
  pass(paste0(f, ": column names match"))

  # Per-column comparison
  col_ok <- TRUE
  for (col in colnames(g)) {
    gv <- g[[col]]
    cv <- c[[col]]

    if (is.numeric(gv) && is.numeric(cv)) {
      # Handle NAs
      na_match <- identical(is.na(gv), is.na(cv))
      if (!na_match) {
        fail(paste0(f, " col '", col, "': NA pattern mismatch"))
        col_ok <- FALSE
        next
      }

      # Compare non-NA values
      valid <- !is.na(gv)
      if (sum(valid) > 0) {
        max_diff <- max(abs(gv[valid] - cv[valid]))
        # Also check relative difference for large values
        rel_denom <- pmax(abs(gv[valid]), abs(cv[valid]), 1e-10)
        max_rel_diff <- max(abs(gv[valid] - cv[valid]) / rel_denom)

        if (max_diff > tolerance && max_rel_diff > tolerance) {
          fail(paste0(f, " col '", col, "': max abs diff = ", signif(max_diff, 4),
                      ", max rel diff = ", signif(max_rel_diff, 4)))
          col_ok <- FALSE
        }
      }
    } else {
      # Character comparison
      if (!identical(as.character(gv), as.character(cv))) {
        n_diff <- sum(as.character(gv) != as.character(cv), na.rm = TRUE)
        fail(paste0(f, " col '", col, "': ", n_diff, " character mismatches"))
        col_ok <- FALSE
      }
    }
  }

  if (col_ok) {
    pass(paste0(f, ": all values match within tolerance"))
  }
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat("Passed:", n_pass, "\n")
cat("Failed:", n_fail, "\n")

if (n_fail > 0) {
  cat("\nFailures:\n")
  for (f in failures) {
    cat("  -", f, "\n")
  }
  quit(status = 1)
} else {
  cat("\nAll checks passed!\n")
  quit(status = 0)
}
