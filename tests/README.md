# msTrawler E2E Regression Tests

End-to-end test harness that captures deterministic outputs from the current msTrawler code and compares future runs against them. Designed to be run **before and during refactoring** to ensure behavioral equivalence.

## Prerequisites

- Docker

## Quick Start

```bash
# From the msTrawler root directory:

# Generate golden baseline (run once before refactoring)
./tests/run_e2e.sh generate

# After making changes, verify nothing broke
./tests/run_e2e.sh test

# Do both (identity check)
./tests/run_e2e.sh both
```

## How It Works

1. **Docker image** (`Dockerfile`) builds R 4.3.2 with all msTrawler dependencies and installs the package from source.
2. **Golden baseline** (`generate_golden_baseline.R`) runs 6 test scenarios with `set.seed(42)`, producing 14 deterministic reference CSVs in `tests/golden/`.
3. **Candidate run** re-executes the same scenarios, writing outputs to `tests/candidate/`.
4. **Comparison** (`compare_outputs.R`) checks every CSV pair for matching file sets, dimensions, column names, and values within `1e-6` tolerance.

## Test Scenarios

| # | Function | What it validates | Key outputs |
|---|----------|-------------------|-------------|
| 1 | `convertPdFile()` | Proteome Discoverer TSV parsing | `converted_df.csv` |
| 2 | `tmtLOD()` | LOD imputation (intensities, SNRs, flags) | `lod_intensities.csv`, `lod_snr.csv`, `lod_toofew.csv` |
| 3 | `geoNorm()` | Geometric mean normalization | `geonorm_intensities.csv`, `geonorm_factors.csv` |
| 4 | `findOutliers()` | Outlier detection matrix | `outlier_matrix.csv` (skipped if < 2 scans) |
| 5 | `msTrawl()` with covariates | Full pipeline (factors + time) | `pipeline/Simple.csv`, `Factor_Age_*.csv`, `Time_*.csv`, `ColAdjustmentFactors.csv` |
| 6 | `msTrawl()` without covariates | Simple mode | `pipeline_simple/Simple.csv`, `ColAdjustmentFactors.csv` |

All scenarios use the bundled tutorial data in `data/tutorial_example/`.

## Determinism

`tmtLOD()` uses `rnorm()` for below-LOD imputation. `set.seed(42)` is called before each test scenario to ensure identical random draws across runs.

## Comparison Checks

For each CSV file, `compare_outputs.R` verifies:

1. Same set of output files exists in both directories
2. Same dimensions (rows x columns)
3. Same column names
4. Numeric columns match within tolerance (default `1e-6`, both absolute and relative)
5. Character columns match exactly

## Adjusting Tolerance

```bash
# Pass a custom tolerance as the third argument to compare_outputs.R
docker run --rm -v "$(pwd):/workspace:rw" mstrawler-e2e \
  Rscript /workspace/tests/compare_outputs.R \
    /workspace/tests/golden \
    /workspace/tests/candidate \
    1e-4
```

## Regenerating the Baseline

If you make an **intentional** behavior change (e.g., fixing a bug in LOD imputation), regenerate the golden baseline:

```bash
./tests/run_e2e.sh generate
```

Then commit the updated `tests/golden/` files alongside your code changes.

## File Structure

```
tests/
├── README.md                         # This file
├── Dockerfile                        # R 4.3.2 + msTrawler build environment
├── generate_golden_baseline.R        # Produces golden reference outputs
├── compare_outputs.R                 # Compares golden vs candidate
├── run_e2e.sh                        # Shell orchestrator
├── golden/                           # Reference outputs (committed to repo)
│   ├── converted_df.csv
│   ├── lod_intensities.csv
│   ├── lod_snr.csv
│   ├── lod_toofew.csv
│   ├── geonorm_intensities.csv
│   ├── geonorm_factors.csv
│   ├── pipeline/
│   │   ├── Simple.csv
│   │   ├── Factor_Age_Old.csv
│   │   ├── Factor_Age_Young.csv
│   │   ├── Time_Old.csv
│   │   ├── Time_Young.csv
│   │   └── ColAdjustmentFactors.csv
│   └── pipeline_simple/
│       ├── Simple.csv
│       └── ColAdjustmentFactors.csv
└── candidate/                        # Test outputs (gitignored)
```
