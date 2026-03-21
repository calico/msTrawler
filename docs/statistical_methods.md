# msTrawler: Statistical Methods and Future Directions

## 1. Problem Statement

### What problem does msTrawler solve?

Isobaric-tagged tandem mass spectrometry (e.g., TMT) enables simultaneous quantification of protein abundances across multiple biological samples in a single experiment. Each sample is labeled with a unique chemical tag (e.g., TMT-126, TMT-127N, ..., TMT-133C), and a mass spectrometer measures the relative intensity of each tag per peptide-spectrum match (PSM).

However, analyzing this data is challenging due to three interconnected problems:

1. **Signal-dependent heteroscedasticity** — Low-intensity measurements have disproportionately high variance. A protein with SNR=2 is far noisier than one with SNR=200, yet naive analyses treat them equally.

2. **Limit-of-detection censoring** — Below a certain signal threshold, measurements become unreliable. These values are not truly zero or missing — they are censored observations from a distribution truncated by instrument sensitivity.

3. **Complex experimental designs** — Real proteomics experiments involve:
   - Multiple plexes (batches) requiring cross-plex normalization via bridge channels
   - Categorical factors (e.g., young vs. old)
   - Continuous covariates (e.g., drug dose)
   - Longitudinal time-series with polynomial and circadian components
   - Nested random effects (repeated measures within subjects)

msTrawler addresses all three problems in a unified pipeline: it imputes censored values with appropriate uncertainty, normalizes across plexes, detects outlier scans, and fits weighted mixed-effects models that account for signal-dependent variance — producing protein-level estimates, contrasts, and hypothesis tests with proper FDR control.

---

## 2. Statistical Methods

### 2.1 Limit-of-Detection (LOD) Imputation

**Goal:** Replace below-threshold values with probabilistic imputations that reflect their uncertainty.

**Threshold determination:**
- *Relative mode* (`lod < 1`): A channel is censored if its proportion of total scan intensity falls below `lod`. For example, `lod = 0.01` flags channels contributing < 1% of scan signal.
- *Absolute mode* (`lod >= 1`): A channel is censored if its signal-to-noise ratio (SNR) falls below `lod`.

**Imputation model:**

Censored intensities are drawn from a log-normal distribution:

```
I_imputed = 2^( Normal(log2(LOD/2), σ) )
```

where `σ = 1 / sqrt(scaleSN × LOD / (2 × imputePenalty))`

Key properties:
- The imputed mean is half the detection threshold (conservative)
- The variance scales inversely with signal strength — matching the observed heteroscedasticity
- `imputePenalty` controls how much imputed values are down-weighted in subsequent modeling
- Imputed SNR values are set deterministically (not random), ensuring the downstream weights properly penalize these observations

**Quality gate:** Scans with fewer than `minAbove` channels above LOD are removed entirely — imputing an entire scan would be meaningless.

### 2.2 Geometric Mean Normalization

**Goal:** Remove systematic per-channel, per-plex biases so that relative abundances are comparable across samples.

**Method:**

1. Select stable proteins for normalization (those with low standard deviation, below the `colAdjust` percentile)
2. For each plex × channel combination, compute the log₂ geometric mean across selected proteins
3. Compute a grand mean across all channels
4. Subtract channel-specific deviations from the grand mean

Mathematically, for plex *p* and channel *j*:

```
normFactor(p,j) = mean(log2(I[normRows, j, p])) - grandMean
I_normalized(i,j,p) = 2^(log2(I(i,j,p)) - normFactor(p,j))
```

This ensures all channels within a plex have the same geometric mean across normalization proteins. An optional `ratios` parameter accommodates known dilution factors (e.g., 2× bridge channels).

### 2.3 Outlier Detection (Consensus Residuals)

**Goal:** Identify aberrant peptide-spectrum matches within each protein grouping.

**Method:** For each protein × plex combination, fit two weighted linear models:

| Model | Formula | Purpose |
|-------|---------|---------|
| Channel model | `I ~ channel + 0` (weights ∝ SNR) | Expects consistent ratios across channels |
| Intercept model | `I ~ 1` (weights ∝ SNR) | Expects uniform intensity |

A scan is flagged as an outlier if its standardized residual exceeds the cutoff (default: 3) in **either** model. This consensus approach catches both channel-specific and global outliers.

Flagged scans are either:
- **Removed** (`swapProtein = FALSE`) — the default
- **Reassigned** (`swapProtein = TRUE`) — the protein label is replaced with the peptide sequence, preserving potential post-translational modification signals

### 2.4 Core Statistical Model

**Response variable:** `log2(Intensity)` — the log-transformed reporter ion abundance

**Weighting scheme:** Inverse-variance weights proportional to signal-to-noise:

```
w_ij = scaleSN × SNR_ij
```

This is the key innovation: observations with higher SNR get more weight, properly handling the signal-dependent heteroscedasticity inherent in isobaric tag data.

**Model hierarchy (adaptive selection):**

msTrawler uses a fallback strategy, attempting the most complex model first:

| Level | Model | When used |
|-------|-------|-----------|
| 1 | Two-level `lmer` | Longitudinal data with sufficient replication |
| 2 | One-level `lmer` | Enough subjects for random intercepts |
| 3 | Weighted `lm` | Fallback when mixed-effects estimation fails |

**Full mixed-effects model (Level 1):**

```
log2(I_ijk) = X_ijk β + Z_j u_j + Z_jk v_jk + ε_ijk
```

where:
- *i* = scan, *j* = subject, *k* = sample (timepoint)
- `β` = fixed effects (factors, continuous covariates, time polynomials, circadian terms)
- `u_j ~ N(0, σ²_subject)` = between-subject random intercept
- `v_jk ~ N(0, σ²_sample)` = within-subject, between-sample random intercept
- `ε_ijk ~ N(0, σ²_tech / w_ijk)` = technical noise (inversely weighted by SNR)

**Bridge channel handling:**

When bridge channels are present (recommended for multi-plex experiments), the design matrix is constructed so bridge observations contribute to scan-level estimation but not to sample-level contrasts — their rows in the sample-level design matrix are zeroed out.

### 2.5 Hypothesis Testing

**For mixed-effects models (`lmer`):**
- **Method:** Kenward-Roger approximation via `pbkrtest::KRmodcomp()`
- **Why:** KR adjusts both the F-statistic and denominator degrees of freedom to account for uncertainty in variance component estimation — critical for the small sample sizes typical in proteomics

**For linear models (`lm`):**
- **Method:** F-test via `car::lht()` (linear hypothesis test)

**Time trend tests:**
- Joint hypothesis tests on all time polynomial coefficients (Time, Time², Time³) and optionally circadian terms (sin, cos)
- When a time-category factor exists (e.g., Age), separate trend tests are constructed per category, plus a differential trend test across categories

### 2.6 Multiple Testing Correction

- **Method:** Benjamini-Hochberg FDR (`p.adjust(method = "fdr")`)
- **Scope:** Applied per contrast column across all proteins
- **Numerical handling:** Zero p-values (from numerical precision limits) are replaced with `10^-300` before correction

---

## 3. Future Directions

### 3.1 Statistical Improvements

**Bayesian mixed-effects models**
- Replace frequentist `lmer` + KR approximation with Bayesian hierarchical models (e.g., `brms`, `Stan`)
- Benefits: natural handling of small samples, informative priors for variance components, posterior credible intervals instead of confidence intervals
- Could incorporate protein-level shrinkage (empirical Bayes) similar to `limma`'s approach for microarrays

**Robust regression**
- The current outlier detection + removal approach is two-stage (detect, then remove). A single-stage robust regression (e.g., M-estimation, Huber weights) would simultaneously down-weight outliers during model fitting
- Consider `robustlmm` package for robust mixed-effects models

**Improved LOD modeling**
- Current log-normal imputation is a single draw — consider multiple imputation (Rubin's rules) to properly propagate imputation uncertainty into downstream estimates
- Alternatively, use a censored regression model (Tobit model) that treats below-LOD values as left-censored observations rather than imputed values
- `survival::survreg()` or `brms` with censored likelihoods can handle this natively

**Empirical Bayes shrinkage**
- Borrow strength across proteins (like `limma`'s moderated t-statistic) to improve variance estimates for proteins with few peptides
- This would substantially improve power for low-abundance proteins

**Non-linear time modeling**
- Current polynomial time modeling (degree 1-3) + sinusoidal circadian terms is rigid
- Consider Gaussian process regression or spline-based models for more flexible temporal dynamics
- Functional data analysis (FDA) approaches could model entire time curves

### 3.2 Computational Improvements

**Parallelization**
- The current `protPrep`/`miniTrawl` approach saves .rda files and requires manual orchestration
- Modern alternatives: `future`/`furrr` for transparent parallelization, or `BiocParallel` for Bioconductor-style backends
- Could parallelize at the protein level trivially since proteins are modeled independently

**GPU acceleration**
- For large experiments (thousands of proteins × hundreds of samples), consider GPU-accelerated linear algebra
- `torch` R package or Python with `jax`/`pytorch` for batch model fitting

**Memory efficiency**
- Current approach loads entire dataset into memory
- For very large experiments, consider streaming or chunked processing
- The `arrow`/`parquet` format could replace CSV for intermediate results

### 3.3 Methodological Extensions

**Peptide-level modeling with protein inference**
- Instead of aggregating to protein level, model peptide-level data directly and use protein inference (e.g., shared peptide models) to improve quantification
- Related: isoform-level quantification using peptide uniqueness

**Missing-not-at-random (MNAR) modeling**
- Low-abundance proteins are more likely to be missing — this is MNAR, not MCAR
- Selection models (Heckman-type) or pattern-mixture models could improve estimates for proteins near the detection limit

**Multi-omics integration**
- Extend the model framework to jointly analyze proteomics + transcriptomics + metabolomics
- Shared latent factors across omics layers could improve power

**Differential variability testing**
- Beyond testing mean differences, test for differences in *variance* between conditions
- Relevant for understanding biological heterogeneity (e.g., drug response variability)

**Network-aware testing**
- Incorporate protein-protein interaction networks or pathway membership as prior information
- Bayesian network-regularized models could improve detection of coordinated changes

### 3.4 Software Engineering Improvements

**Python port (partial)**
- Preprocessing (LOD, normalization, outlier detection) translates cleanly to `pandas`/`numpy`/`scipy`
- The main bottleneck: `lme4` + `pbkrtest` have no equivalent-quality Python implementations
- Hybrid approach: Python preprocessing + `rpy2` bridge for mixed-effects models
- Long-term: `statsmodels` mixed-effects is improving; `bambi`/`PyMC` for Bayesian alternatives

**Standardized output format**
- Replace CSV output with structured objects (S4/R6 classes) that support `summary()`, `plot()`, `coef()` methods
- Consider `SummarizedExperiment` or `MSnSet` containers for Bioconductor interoperability

**Comprehensive testing**
- Unit tests for each statistical function (`testthat`)
- Property-based testing: verify that outputs satisfy expected statistical properties (e.g., FDR control under null simulation)
- Integration with continuous integration (GitHub Actions)

---

## 4. Summary

| Component | Current Method | Potential Improvement |
|-----------|---------------|----------------------|
| LOD handling | Single log-normal imputation | Multiple imputation or censored regression (Tobit) |
| Normalization | Geometric mean (stable proteins) | Quantile normalization, or model-based normalization |
| Outlier detection | Consensus standardized residuals | Robust mixed-effects (M-estimation) |
| Core model | Weighted lmer with KR tests | Bayesian hierarchical model with protein-level shrinkage |
| Time modeling | Polynomial + sinusoidal | Gaussian processes or splines |
| Multiple testing | BH FDR | Local FDR (lfdr) or Bayesian FDR |
| Parallelization | Manual .rda splitting | `future`/`furrr` or `BiocParallel` |
| Missing data | Impute + remove | MNAR-aware selection models |
