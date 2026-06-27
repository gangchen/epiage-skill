---
name: grimage-calculator
description: >-
  Compute GrimAge (and other epigenetic / DNA-methylation aging clocks like
  Horvath, Hannum, PhenoAge) from a DNA methylation beta-value file.
  Use this whenever the user provides methylation array data (a CSV/TSV of CpG
  sites with beta values, e.g. from an EPIC/450K array) and wants to know their
  biological age, epigenetic age, GrimAge, "甲基化年龄", aging clock, age
  acceleration, or how old their DNA "looks" — even if they just drop a
  methylation file and ask "how old am I biologically". Also triggers on
  "GrimAge", "DNAm age", "epigenetic clock", "表观遗传时钟", "生物年龄",
  "衰老时钟", or comparing several aging clocks on one sample.
---

# GrimAge / Epigenetic Clock Calculator

Computes epigenetic aging clocks from a DNA methylation beta-value file. GrimAge
is the default and primary clock.

**Self-contained.** This skill needs only `pandas` + `numpy` — no biolearn, torch,
seaborn, or network. The clock coefficients and the imputation reference are
vendored under `data/` (extracted from the open-source biolearn library and
trimmed to just the ~1900 CpGs the clocks use, ~150 KB total). The math in
`scripts/compute_grimage.py` is a faithful reimplementation of biolearn's
`GrimageModel` and `LinearMethylationModel`, verified to reproduce biolearn's
outputs to 4 decimals across all five clocks.

## What GrimAge is (so you interpret it correctly)

GrimAge is a **second-generation** clock: instead of predicting chronological age
(like the first-generation Horvath clock), it predicts **mortality/aging risk**,
expressed in year units. It works by first estimating DNAm surrogates of 7
plasma proteins (ADM, B2M, cystatin C, GDF15, leptin, PAI-1, TIMP1) plus
DNAm smoking pack-years, then combining those with age and sex in a survival
model. **V2** additionally adds DNAm logA1C (blood sugar) and logCRP
(inflammation) and is the recommended version.

Because age and sex feed directly into the formula, **GrimAge cannot be computed
without the person's chronological age and biological sex.** Always obtain both
before running. If the user only gives an age range, use the midpoint and run the
`--sensitivity` table so they see how much the answer depends on the exact age.

## Workflow

### 1. Inspect the input file
Look at the first few lines. Two layouts are supported and auto-detected:
- **Long**: two columns — CpG id, then beta value (header names don't matter).
- **Matrix**: first column = CpG id, remaining columns = one or more samples.

Beta values are floats in [0,1]. Note the row count; a typical EPIC array has
~865k probes, 450K ~485k. Coverage of GrimAge's ~1030 clock CpGs is what matters,
and the script reports it.

### 2. Get age and sex
GrimAge requires both. Ask the user if not already known. Sex maps to
`--sex m` or `--sex f`.

### 3. Check dependencies
Only `pandas` and `numpy` are needed, and they're usually already present. Verify:

```bash
python3 -c "import pandas, numpy" 2>/dev/null && echo OK || pip install pandas numpy
```

No venv, no biolearn/torch, no network required. Just run the script with the
system `python3`.

### 4. Run the computation
```bash
python3 <skill-dir>/scripts/compute_grimage.py \
  --input "<betas.csv>" --age <years> --sex <m|f> \
  --sensitivity 40 42 47 49        # optional, when exact age is uncertain
```
Default runs GrimAge V1 + V2. To add other clocks:
`--clocks grimage horvath hannum phenoage`.

The script: loads betas, imputes any missing clock CpGs from the bundled
(trimmed) `data/sesame_450k_median.csv` population reference, runs each model, and
prints a table of value, acceleration (clock − chronological age), and CpG
coverage, plus a JSON line for downstream use.

### 5. Report and interpret
Present the GrimAge value(s) and acceleration. Then convey these caveats, which
genuinely matter for honest interpretation:

- **Open-source reproduction, not the official calculator.** biolearn reimplements
  GrimAge; values track Horvath's official server/Clock Foundation closely but may
  differ slightly (coefficient provenance + preprocessing/imputation differences).
  It's a reliable estimate, not an officially certified number. If the user needs
  a citable/authoritative value, point them to the Horvath DNAm Age calculator or a
  commercial provider (e.g. TruDiagnostic / Clock Foundation).
- **"Acceleration" here = clock − chronological age**, a simple difference. The
  academic *GrimAgeAccel* (residual vs. a same-age cohort) needs a population
  sample and **cannot be computed for a single person** — so don't phrase a +8 as
  "8 years older than peers"; say "epigenetic-predicted age is 8 years above
  chronological age."
- **Don't over-read GrimAge as looking-your-age.** It's a mortality/aging-risk
  score in year units, not "your DNA looks N years old." For the "猜年龄" question,
  Horvath/Hannum are the right clocks.
- **Age sensitivity**: GrimAge moves with the assumed chronological age; if age was
  approximate, show the sensitivity table and ask for the exact age to finalize.
- **Imputation**: a small number of missing clock CpGs filled from a population
  median is normal and low-impact; flag it if coverage drops below ~90%.

## Notes
- Multiple samples (matrix input) are handled — one result row per sample.
- For first-generation clocks (Horvath/Hannum/PhenoAge) age/sex are not required,
  but pass `--age` anyway so acceleration is reported.
- **DunedinPACE is intentionally not included.** Unlike the additive clocks here,
  it requires quantile-normalizing the sample against a gold-standard reference
  distribution, which can't be done faithfully from a single sample without the
  heavier biolearn machinery. If a user specifically needs DunedinPACE (a *pace of
  aging* rate, not an age), tell them it needs the full biolearn install.
- **Provenance.** Coefficients in `data/` are biolearn's, which reimplements the
  published clocks. GrimAge values track Horvath's official calculator closely but
  aren't an officially certified number (see the interpretation caveats above).
