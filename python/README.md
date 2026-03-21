# msTrawler (Python)

Python implementation of the msTrawler proteomics analysis pipeline. Analyzes multiplexed (TMT) mass spectrometry data with signal-aware weighted modeling.

Supports two model-fitting backends:
- **`backend="python"`** (default) — pure Python, no R dependency, uses statsmodels WLS + Wald F-tests
- **`backend="r"`** — uses R's lme4 + pbkrtest via rpy2 for mixed-effects models + Kenward-Roger tests

## Installation

```bash
cd python

# Pure Python (default — no R needed)
pip install -e ".[dev]"

# With R backend support (requires R + lme4 + pbkrtest + car installed)
pip install -e ".[r,dev]"
```

**Core dependencies:** numpy, pandas, scipy, statsmodels (all from PyPI).
**Optional R dependency:** rpy2 (plus R runtime with lme4, pbkrtest, car packages).

## Quick Start

```python
from mstrawler.file_converter import convert_pd_file
from mstrawler.pipeline import ms_trawl
import pandas as pd

# Step 1: Convert Proteome Discoverer output
df = convert_pd_file("pd_export_PSMs.txt", "pd_export_Proteins.txt")

# Step 2: Load sample and covariate metadata
sample_file = pd.read_csv("sample_file.csv")
covariate_file = pd.read_csv("covariate_file.csv")

# Step 3: Run the full pipeline (pure Python — default)
ms_trawl(
    df=df,
    sample_file=sample_file,
    covariate_file=covariate_file,
    lod=0.01,
    min_above=4,
    col_adjust=0.5,
    seed=42,
    backend="python",  # default — no R needed
)

# Or: use R backend for mixed-effects models + Kenward-Roger tests
ms_trawl(
    df=df,
    sample_file=sample_file,
    covariate_file=covariate_file,
    lod=0.01,
    min_above=4,
    col_adjust=0.5,
    seed=42,
    backend="r",  # requires R + rpy2 + R packages
)
```

### Check backend availability

```python
from mstrawler.model_fitting import check_backend

print(check_backend("python"))
# {'available': True, 'backend': 'python', 'details': 'statsmodels WLS + Wald F-tests (always available)'}

print(check_backend("r"))
# {'available': True/False, 'backend': 'r', 'details': '...'}
```
```

Output CSVs are written to the current directory:
- `Simple.csv` — per-sample protein abundance estimates
- `Factor_<Name>_<Level>.csv` — factor contrasts (log2 fold-change, p-value, q-value)
- `Time_<Level>.csv` — time trend estimates and tests
- `ColAdjustmentFactors.csv` — normalization factors applied

## Input Files

### PSM file (from Proteome Discoverer)

Tab-delimited export containing columns:
- `Protein Accessions` or `Master Protein Accessions`
- `Annotated Sequence`
- `Modifications`
- `File ID` (plex identifier)
- `Abundance` columns (one per TMT channel)
- `Average Reporter S/N`

### Sample file (CSV)

Maps TMT channels to experimental conditions:

```csv
SampleID,Bridge,MouseID,Age,Time
F1_X126.Sn,1,Bridge,Bridge,NA
F1_X127n.Sn,0,39,Young,2
F1_X127c.Sn,0,39,Young,4
...
```

### Covariate file (CSV)

Defines the type of each covariate:

```csv
Covariate,Type,Levels,TimeDegree,TimeCategory,Circadian
MouseID,id,0,0,0,0
Age,factor,2,0,1,0
Time,time,0,0,0,1
```

| Type | Meaning |
|------|---------|
| `id` | Subject identifier for random effects grouping |
| `factor` | Categorical variable (e.g., Young vs Old) |
| `continuous` | Numeric covariate (centered automatically) |
| `time` | Time variable for longitudinal modeling |

## API Reference

### `convert_pd_file(pd_psm_file, pd_protein_file="", delimiter="\t")`

Converts Proteome Discoverer exports to an msTrawler DataFrame.

### `ms_trawl(df, sample_file, covariate_file, ...)`

Full analysis pipeline. Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `scale_sn` | 1.0 | SNR scaling factor for weights |
| `lod` | 0.01 | Limit of detection (< 1 = proportion, >= 1 = SNR) |
| `impute_penalty` | 1.0 | Penalty for imputed values |
| `min_above` | 4 | Minimum channels above LOD per scan |
| `ssn_filter` | 20.0 | Minimum summed SNR per scan |
| `outlier_cutoff` | 3.0 | Standardized residual threshold for outliers |
| `col_adjust` | 0.5 | SD percentile for normalization protein selection |
| `seed` | None | Random seed for reproducible LOD imputation |
| `backend` | `"python"` | Model fitting backend: `"python"` (statsmodels) or `"r"` (lme4 via rpy2) |

### Individual functions

Each pipeline step is also available as a standalone function:

```python
from mstrawler.preprocessing import tmt_lod, geo_norm, find_outliers
from mstrawler.covariate_setup import setup_covariates
from mstrawler.restructure import restructure_for_modeling
from mstrawler.model_fitting import adapt_model
from mstrawler.results import init_result_tables, write_result_tables
```

## Pipeline Steps

```
Raw PSM data
  │
  ▼
1. convert_pd_file()     Parse Proteome Discoverer exports
  │
  ▼
2. tmt_lod()             Impute below-LOD values (log-normal, SNR-penalized)
  │
  ▼
3. geo_norm()            Geometric mean normalization across channels/plexes
  │
  ▼
4. find_outliers()       Consensus outlier detection (dual WLS models)
  │
  ▼
5. setup_covariates()    Parse factors, time terms, circadian, build formula
  │
  ▼
6. restructure()         Melt wide→long, merge covariates, compute weights
  │
  ▼
7. adapt_model()         Per-protein weighted model fitting + hypothesis tests
  │
  ▼
8. write_results()       FDR correction (Benjamini-Hochberg) + CSV output
```

## Testing

```bash
cd python

# Run all tests (13 E2E tests against R golden baseline)
pytest tests/ -v

# Run specific test suites
pytest tests/e2e/test_file_converter.py -v     # File conversion (4 tests)
pytest tests/e2e/test_preprocessing.py -v       # LOD + normalization (5 tests)
pytest tests/e2e/test_pipeline.py -v            # Full pipeline (4 tests)
```

Tests compare Python outputs against golden baseline CSVs generated by the original R package, ensuring behavioral equivalence.

## Choosing a Backend

The `backend` parameter controls which statistical engine is used for model fitting. Everything else (preprocessing, normalization, outlier detection, data restructuring, FDR correction) runs in pure Python regardless of backend.

### `backend="python"` (default)

| Aspect | Details |
|--------|---------|
| **Linear models** | `statsmodels.WLS` with SNR-based weights |
| **Mixed-effects models** | Not available (WLS fallback) |
| **Hypothesis tests** | Wald F-test |
| **Dependencies** | numpy, pandas, scipy, statsmodels |
| **Best for** | Production deployment, large experiments, pipelines |

### `backend="r"`

| Aspect | Details |
|--------|---------|
| **Linear models** | R's `lm()` with weights |
| **Mixed-effects models** | `lme4::lmer` (nested random intercepts + weights) |
| **Hypothesis tests** | Kenward-Roger via `pbkrtest` |
| **Dependencies** | rpy2 + R runtime + R packages (lme4, pbkrtest, car, msTrawler) |
| **Best for** | Small-sample experiments, publication-quality statistics |

### When does the choice matter?

| Scenario | Recommended backend |
|----------|-------------------|
| Large experiments (> 10 subjects per group) | `"python"` — WLS converges to lmer as n grows |
| Small experiments (3-5 subjects per group) | `"r"` — KR correction and mixed effects matter here |
| Experiments with few scans per protein | Either — both fall back to fixed-effects with 1-2 scans |
| Experiments with many scans per protein | `"r"` — mixed effects properly partition technical vs biological variance |
| Integration into a web service or pipeline | `"python"` — no R runtime in Docker |
| Publication where reviewers expect mixed models | `"r"` — lmer + KR is the established standard |
| Exploratory analysis / rapid iteration | `"python"` — faster, simpler setup |

### Architecture

```
                    ms_trawl(backend="python" | "r")
                           │
    ┌──────────────────────┼──────────────────────┐
    │                      │                      │
    │  Preprocessing       │  Model Fitting       │  Results
    │  (always Python)     │  (backend choice)    │  (always Python)
    │                      │                      │
    │  - LOD imputation    │  "python":           │  - Table init
    │  - Normalization     │    statsmodels WLS   │  - FDR correction
    │  - Outlier detect    │    Wald F-tests      │  - CSV output
    │  - Data reshape      │                      │
    │                      │  "r":                │
    │                      │    lme4::lmer        │
    │                      │    pbkrtest::KR      │
    │                      │    car::lht          │
    └──────────────────────┴──────────────────────┘
```

## Roadmap

Future enhancements to the model-fitting layer:

1. **Bayesian mixed models via Bambi/PyMC** — a third `backend="bayesian"` option with no R dependency, better than KR for small samples, and native uncertainty quantification.

2. **`statsmodels.MixedLM` with weights** — once [statsmodels issue #9342](https://github.com/statsmodels/statsmodels/issues/9342) is resolved, the `"python"` backend will support observation-level weights in mixed models.

3. **Bootstrap p-values** — resample-based inference as a backend-agnostic alternative to both KR and Wald tests.

## Docker

Two Docker images are provided — build both from the **repo root** (not from `python/`).

### Pure Python image (no R)

```bash
# Build
docker build -t mstrawler-python -f python/Dockerfile .

# Run tests
docker run --rm mstrawler-python

# Run analysis on your data (mount data directory)
docker run --rm \
  -v /path/to/your/data:/data \
  -v /path/to/output:/output \
  -w /output \
  mstrawler-python \
  python -c "
from mstrawler.file_converter import convert_pd_file
from mstrawler.pipeline import ms_trawl
import pandas as pd

df = convert_pd_file('/data/PSMs.txt', '/data/Proteins.txt')
sample_file = pd.read_csv('/data/sample_file.csv')
covariate_file = pd.read_csv('/data/covariate_file.csv')

ms_trawl(df, sample_file, covariate_file, seed=42, backend='python')
"

# Interactive Python shell
docker run --rm -it \
  -v /path/to/your/data:/data \
  mstrawler-python \
  python
```

| Property | Value |
|----------|-------|
| Base image | `python:3.12-slim` |
| Size | ~300 MB |
| Backend | `"python"` only |
| R required | No |

### Python + R image (both backends)

```bash
# Build (takes longer — installs R + lme4 + pbkrtest)
docker build -t mstrawler-r -f python/Dockerfile.r .

# Run tests
docker run --rm mstrawler-r

# Run analysis with R backend
docker run --rm \
  -v /path/to/your/data:/data \
  -v /path/to/output:/output \
  -w /output \
  mstrawler-r \
  python -c "
from mstrawler.file_converter import convert_pd_file
from mstrawler.pipeline import ms_trawl
import pandas as pd

df = convert_pd_file('/data/PSMs.txt', '/data/Proteins.txt')
sample_file = pd.read_csv('/data/sample_file.csv')
covariate_file = pd.read_csv('/data/covariate_file.csv')

ms_trawl(df, sample_file, covariate_file, seed=42, backend='r')
"

# Check which backends are available
docker run --rm mstrawler-r python -c "
from mstrawler.model_fitting import check_backend
print(check_backend('python'))
print(check_backend('r'))
"
```

| Property | Value |
|----------|-------|
| Base image | `python:3.12-slim` + R |
| Size | ~1.5 GB |
| Backend | Both `"python"` and `"r"` |
| R required | Included in image |

### Running the tutorial example in Docker

```bash
# Using the bundled tutorial data
docker run --rm \
  -v $(pwd)/output:/output \
  -w /output \
  mstrawler-python \
  python -c "
from mstrawler.file_converter import convert_pd_file
from mstrawler.pipeline import ms_trawl
import pandas as pd

df = convert_pd_file(
    '/repo/data/tutorial_example/pd_example_export_PSMs.txt',
    '/repo/data/tutorial_example/pd_example_export_Proteins.txt',
)
sample_file = pd.read_csv('/repo/data/tutorial_example/sample_file.csv')
covariate_file = pd.read_csv('/repo/data/tutorial_example/covariate_file.csv')

ms_trawl(df, sample_file, covariate_file, seed=42)
print('Done! Check /output for CSV files.')
"

ls output/
# Simple.csv  Factor_Age_Old.csv  Factor_Age_Young.csv  Time_Old.csv  Time_Young.csv  ColAdjustmentFactors.csv
```

## Project Structure

```
python/
├── pyproject.toml                     # Package metadata and dependencies
├── README.md                          # This file
├── Dockerfile                         # Pure Python image (no R)
├── Dockerfile.r                       # Python + R image (both backends)
├── src/mstrawler/
│   ├── __init__.py                    # Package version
│   ├── file_converter.py              # Proteome Discoverer → DataFrame
│   ├── preprocessing.py               # LOD imputation, normalization, outliers
│   ├── covariate_setup.py             # Factor/time/circadian configuration
│   ├── restructure.py                 # Wide → long data reshaping
│   ├── model_fitting.py               # Backend dispatcher (python/r)
│   ├── model_fitting_python.py        # Pure Python backend (statsmodels WLS)
│   ├── model_fitting_r.py             # R backend (lme4 via rpy2)
│   ├── results.py                     # Table initialization, FDR, CSV output
│   └── pipeline.py                    # ms_trawl() orchestrator
└── tests/
    ├── e2e/
    │   ├── test_file_converter.py     # 4 tests vs R golden baseline
    │   ├── test_preprocessing.py      # 5 tests vs R golden baseline
    │   └── test_pipeline.py           # 4 tests vs R golden baseline
    └── unit/                          # (future unit tests)
```
