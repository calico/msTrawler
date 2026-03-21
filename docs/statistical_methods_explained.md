# msTrawler Statistical Methods — A Beginner's Guide

## What does msTrawler do?

Imagine you're a scientist studying how proteins change in mice as they age. You collect tissue samples from young and old mice, and you want to measure thousands of proteins at once.

**TMT (Tandem Mass Tag)** technology lets you do this: you label each sample with a unique chemical tag, mix them together, and run them through a mass spectrometer. The machine reads out the intensity of each tag for each protein fragment (called a "peptide-spectrum match" or PSM), giving you a relative measurement of how much protein is in each sample.

The problem? These measurements are messy:
- Some signals are too faint to trust
- Different batches have different baselines
- Some scans are just wrong (outliers)
- You have complex experimental designs (young vs. old, multiple timepoints, repeated subjects)

**msTrawler is a pipeline that cleans up this data and fits statistical models to answer questions like: "Is this protein more abundant in old mice than young mice?"**

Let's walk through each step with simple examples.

---

## Step 1: Limit-of-Detection (LOD) Imputation

### The problem

Your mass spectrometer can't measure arbitrarily small amounts of protein. Below a certain threshold, the signal is just noise. But you can't simply throw these values away — that would bias your results toward proteins that happen to have high signals.

### An analogy

Think of a bathroom scale that only shows weights above 5 kg. If you weigh a kitten and it reads "0", that doesn't mean the kitten weighs nothing — it means the kitten weighs *somewhere between 0 and 5 kg*. You need to fill in a reasonable guess.

### What msTrawler does

For each "below threshold" value, msTrawler draws a random replacement from a log-normal distribution centered at half the detection limit.

**Example:**

Suppose we have 4 samples and the LOD is SNR = 10:

```
Before imputation:
  Sample A: SNR = 50  (above LOD, keep as-is)
  Sample B: SNR = 120 (above LOD, keep as-is)
  Sample C: SNR = 3   (below LOD! needs imputation)
  Sample D: SNR = 80  (above LOD, keep as-is)
```

For Sample C, msTrawler:
1. Sets the imputed center at LOD/2 = 5
2. Draws a random value: `2^(Normal(log2(5), σ))` where σ is small
3. Also assigns a low SNR (= 5) so this value gets *low weight* in later modeling

```
After imputation:
  Sample A: SNR = 50,  Intensity = 4200  (original)
  Sample B: SNR = 120, Intensity = 9800  (original)
  Sample C: SNR = 5,   Intensity = 310   (imputed — note the low SNR!)
  Sample D: SNR = 80,  Intensity = 6500  (original)
```

The key insight: the imputed value comes with a built-in "uncertainty penalty" — its low SNR means it will count for very little when we fit our models later. This is much better than either dropping it (losing information) or treating it as a real measurement (adding noise).

### Quality gate

If almost *all* samples in a scan are below LOD, the entire scan is removed. There's simply not enough real signal to work with.

---

## Step 2: Normalization

### The problem

Different batches (called "plexes") can have different overall signal levels. Even within a plex, some channels may be systematically brighter or dimmer than others due to labeling efficiency or instrument drift.

### An analogy

Imagine comparing test scores from two different teachers. Teacher A grades on a 100-point scale, Teacher B on a 50-point scale. Before comparing students, you need to put everyone on the same scale.

### What msTrawler does

It uses **geometric mean normalization**: find a set of "stable" proteins (ones that shouldn't change between samples), calculate the geometric mean of their intensities in each channel, and adjust all channels so these means are equal.

**Example:**

Suppose we have 3 channels in a plex, and we've identified 2 stable housekeeping proteins:

```
Before normalization:
                 Channel A   Channel B   Channel C
  Protein 1:      100          200          150
  Protein 2:      80           170          110

  Geometric mean: 89.4         184.4        128.5
  Grand mean:     127.5
```

Channel B is systematically high (184.4 vs. grand mean 127.5). The normalization factor for Channel B is `184.4 / 127.5 = 1.45`. We divide all Channel B values by 1.45.

```
After normalization:
                 Channel A   Channel B   Channel C
  Protein 1:      100          138          150
  Protein 2:      80           117          110
```

Now the stable proteins have similar intensities across channels, and any remaining differences in *other* proteins reflect real biology, not technical artifacts.

---

## Step 3: Outlier Detection

### The problem

Occasionally a single scan produces a wildly wrong measurement — maybe a contaminant co-eluted with the peptide, or there was an electrical glitch. These outliers can distort model estimates for the entire protein.

### An analogy

If you're calculating the average income in a room and Bill Gates walks in, the average becomes meaningless. You need to detect and handle the outlier.

### What msTrawler does

For each protein in each plex, it fits two simple models and checks if any scan has unusually large residuals (the difference between observed and expected).

**Model 1 — Channel model:** "Each channel should have its own consistent level"
```
Intensity ~ Channel    (no intercept, weights = SNR)
```

**Model 2 — Intercept model:** "All channels should be roughly equal"
```
Intensity ~ 1          (intercept only, weights = SNR)
```

A scan is flagged as an outlier if its **standardized residual** exceeds 3 in *either* model. (A standardized residual of 3 means the observation is 3 standard deviations from what the model expects — extremely unlikely by chance.)

**Example:**

```
Protein X, Plex 1 (3 scans, 4 channels):

                Ch-A    Ch-B    Ch-C    Ch-D
  Scan 1:       100     105     98      102     ← looks normal
  Scan 2:       95      110     100     97      ← looks normal
  Scan 3:       100     500     95      101     ← Ch-B is way off!

Channel model residual for Scan 3, Ch-B: standardized = 4.2 → OUTLIER
```

Scan 3 is flagged and removed (or reassigned to peptide-level analysis if the user wants to look for post-translational modifications).

---

## Step 4: The Core Statistical Model

This is the heart of msTrawler. After cleaning the data, we fit a statistical model for *each protein* to estimate abundances and test for differences between groups.

### 4.1 The response variable

We model **log₂(Intensity)** rather than raw intensity. Why?

- Protein abundances span orders of magnitude (some proteins are 1000× more abundant than others)
- Log-transform makes the distribution more symmetric
- Differences on the log scale correspond to *fold changes*, which is what biologists care about
  - A log₂ difference of 1 = 2-fold change
  - A log₂ difference of 2 = 4-fold change

### 4.2 The weighting scheme (the key innovation)

Not all measurements are equally reliable. A measurement with SNR = 100 is much more precise than one with SNR = 5. msTrawler assigns each observation a weight:

```
weight = scaleSN × SNR
```

This is called **inverse-variance weighting**: higher-SNR observations get more influence on the model fit. Think of it as telling the model: "trust the loud, clear signals more than the faint, noisy ones."

**Example:**

```
Observation 1: log2(Intensity) = 8.5,  SNR = 100  →  weight = 100
Observation 2: log2(Intensity) = 7.2,  SNR = 5    →  weight = 5
```

Observation 1 has 20× more influence on the estimated protein abundance. If Observation 2 was an imputed value (from Step 1), this is exactly what we want — the uncertain imputation barely affects the result.

### 4.3 Simple example: No covariates

The simplest case: estimate how much of each protein is in each sample.

```
Model: log2(Intensity) ~ Sample    (one coefficient per sample, weights = SNR)
```

For a protein measured in 3 samples with 10 scans each, this gives us 3 estimates (one per sample) that are weighted averages of the log-intensities across scans.

**Output (Simple.csv):**
```
Protein   est_SampleA   est_SampleB   est_SampleC
ProteinX  12.3          11.8          12.5
ProteinY  8.1           9.4           8.0
```

These are log₂ abundance estimates. ProteinX has roughly equal abundance across samples; ProteinY is higher in Sample B.

### 4.4 Factor covariates: Young vs. Old

Now suppose we have a factor: Age (Young vs. Old). We want to test whether each protein differs between groups.

```
Model: log2(Intensity) ~ Age + ... + (1|Subject)
       weights = SNR
```

This estimates a **contrast**: log₂(Old) − log₂(Young). If this equals 1, it means the protein is 2-fold more abundant in old mice.

**Example output (Factor_Age_Old.csv):**
```
Protein    Est_Young    Pval_Young    Qval_Young
ProteinX   -0.15        0.72          0.89        ← no difference
ProteinY   1.82         0.003         0.01        ← 3.5-fold higher in Young
ProteinZ   -0.95        0.01          0.04        ← 1.9-fold lower in Young
```

The `Est` column is the log₂ fold-change (Young relative to Old, since Old is the reference). The `Pval` is the raw p-value, and `Qval` is the FDR-adjusted p-value.

### 4.5 Mixed-effects models: Why random effects?

In a typical proteomics experiment:
- You have multiple **scans** per protein per sample (technical replicates)
- You may have multiple **timepoints** per subject (biological replicates)

These create a hierarchical structure: scans are nested within samples, samples are nested within subjects. A simple linear model would assume all observations are independent, but they're not — scans from the same sample are more similar to each other than scans from different samples.

**Mixed-effects models** handle this by adding **random effects**:

```
Full model:
  log2(I_ijk) = [Fixed effects] + u_j + v_jk + ε_ijk

  where:
    i = scan (technical replicate)
    j = subject (biological unit)
    k = sample/timepoint

    u_j  ~ Normal(0, σ²_subject)   ← variation between subjects
    v_jk ~ Normal(0, σ²_sample)    ← variation between samples within a subject
    ε_ijk ~ Normal(0, σ²_tech/w)   ← technical noise (weighted by SNR)
```

**Why does this matter?**

Without random effects, you might claim a protein is significantly different between Young and Old mice, when really you just have one unusual mouse driving the entire effect. Random effects properly account for the fact that your "sample size" is the number of mice, not the number of scans.

### 4.6 Adaptive model selection

Not every protein has enough data for the full mixed-effects model. msTrawler adapts automatically:

```
Try Level 1: Two-level mixed model (subject + sample random effects)
  ↓ fails? (not enough subjects or convergence issues)
Try Level 2: One-level mixed model (sample random effects only)
  ↓ fails?
Try Level 3: Simple weighted linear model (no random effects)
```

The output tells you which model was used (`modelFit` column), so you know the level of confidence for each protein.

### 4.7 Time-series modeling

For longitudinal experiments (e.g., samples collected every 2 hours), msTrawler can fit polynomial time trends:

```
Linear:     log2(I) ~ ... + Time
Quadratic:  log2(I) ~ ... + Time + Time²
Cubic:      log2(I) ~ ... + Time + Time² + Time³
```

It can also add **circadian terms** for 24-hour rhythms:

```
log2(I) ~ ... + sin(2π/24 × Time) + cos(2π/24 × Time)
```

And if you have a factor like Age, it can test whether the time trend *differs* between Young and Old mice (interaction terms).

---

## Step 5: Hypothesis Testing

### How p-values are calculated

After fitting the model, we need to test hypotheses like "Is the Age effect significantly different from zero?"

**For mixed-effects models:** msTrawler uses the **Kenward-Roger (KR) approximation**. This is important because the standard methods for calculating p-values (like the Wald test) don't work well with small sample sizes — they give p-values that are too small (making you overconfident). The KR correction adjusts for this.

**For simple linear models:** Standard F-tests are used.

### What the KR approximation does, intuitively

In a mixed-effects model, you're estimating variance components (how much variation comes from subjects vs. samples vs. technical noise). These estimates are themselves uncertain. The KR approximation says: "I'm not sure exactly how much between-subject variance there is, so I should be *less confident* in my p-values." It inflates the p-values to account for this uncertainty — a more honest assessment.

---

## Step 6: Multiple Testing Correction (FDR)

### The problem

If you test 10,000 proteins, even at p < 0.05, you'd expect 500 false positives by chance alone. That's useless.

### The solution: Benjamini-Hochberg FDR

The **False Discovery Rate (FDR)** controls the *proportion* of your significant results that are false positives. At FDR = 0.05, you accept that about 5% of the proteins you call "significant" may be false discoveries — but the other 95% are real.

**Example:**

```
Raw p-values for 5 proteins (sorted):
  ProteinA: p = 0.001
  ProteinB: p = 0.008
  ProteinC: p = 0.03
  ProteinD: p = 0.12
  ProteinE: p = 0.45

BH-adjusted q-values:
  ProteinA: q = 0.005   ← significant at FDR 5%
  ProteinB: q = 0.02    ← significant at FDR 5%
  ProteinC: q = 0.05    ← borderline
  ProteinD: q = 0.15    ← not significant
  ProteinE: q = 0.45    ← not significant
```

The Q-value is always ≥ the raw p-value. It's the minimum FDR at which you would call this protein significant.

---

## Putting It All Together

Here's the complete pipeline for a typical experiment with Young vs. Old mice, sampled at multiple timepoints:

```
Raw PSM data (from mass spectrometer)
  │
  ▼
Step 1: LOD Imputation
  • Replace low signals with conservative random draws
  • Assign low SNR weights to imputed values
  │
  ▼
Step 2: Normalization
  • Equalize channel baselines using stable proteins
  │
  ▼
Step 3: Outlier Detection
  • Flag and remove scans with extreme residuals
  │
  ▼
Step 4: For EACH protein, fit:
  log2(Intensity) ~ Age + Time + Time² + sin + cos + Age:Time + ...
                    + (1|Subject) + (1|Subject:Sample)
                    weights = SNR
  │
  ▼
Step 5: Hypothesis Tests
  • "Is Age effect ≠ 0?" → Kenward-Roger p-value
  • "Is time trend ≠ 0?" → joint F-test on all time terms
  • "Do trends differ by Age?" → interaction test
  │
  ▼
Step 6: FDR Correction
  • Adjust p-values across all proteins
  │
  ▼
Output CSVs:
  • Simple.csv — abundance per sample per protein
  • Factor_Age_Young.csv — Age contrasts (log₂ fold-change, p-value, q-value)
  • Time_Young.csv — time trend estimates and tests for Young mice
  • Time_Old.csv — time trend estimates and tests for Old mice
```

---

## Glossary

| Term | Meaning |
|------|---------|
| **PSM** | Peptide-Spectrum Match — one measurement of one peptide fragment |
| **TMT** | Tandem Mass Tag — chemical labels for multiplexed quantification |
| **Plex** | A batch of samples processed together (typically 15-16 channels) |
| **SNR** | Signal-to-Noise Ratio — how much real signal vs. background noise |
| **LOD** | Limit of Detection — the threshold below which measurements are unreliable |
| **Bridge channel** | A reference sample included in every plex to enable cross-plex normalization |
| **Log₂ fold-change** | The difference in log₂ abundance between two conditions; 1 = 2× change |
| **Random effect** | A model component that accounts for grouped variation (e.g., between subjects) |
| **Fixed effect** | A model component for the variables you're testing (e.g., Age, Time) |
| **KR approximation** | Kenward-Roger method for computing p-values in mixed models with small samples |
| **FDR** | False Discovery Rate — the expected proportion of false positives among significant results |
| **Q-value** | The FDR-adjusted p-value; the minimum FDR at which this result would be significant |
| **Heteroscedasticity** | When measurement variance depends on the signal level (louder signals are more precise) |
