# epiage-skill

Agent skill(s) for computing **epigenetic aging clocks** from a DNA methylation
beta-value file — entirely offline, with only `pandas` + `numpy`.

## Skills in this repo

### `grimage-calculator`

Compute **GrimAge** (V1 & V2) and several other DNA-methylation aging clocks
(**Horvath**, **Hannum**, **PhenoAge**) from a CpG beta-value CSV (e.g. an
Illumina EPIC / 450K array export).

- **Self-contained**: depends only on `pandas` + `numpy`. No `biolearn`, no
  `torch`, no network. Clock coefficients and the imputation reference are
  vendored (trimmed to the ~1900 CpGs the clocks use, ~150 KB).
- **Faithful**: the math reimplements [biolearn](https://bio-learn.github.io/)'s
  `GrimageModel` and `LinearMethylationModel` and reproduces biolearn's outputs
  to 4 decimal places across all five clocks.

## Install

```bash
npx skills add gangchen/epiage-skill
```

This installs the `grimage-calculator` skill into your agent. Once installed,
just give your agent a methylation file and ask for your GrimAge / biological age.

You can also run the script directly:

```bash
python3 grimage-calculator/scripts/compute_grimage.py \
  --input betas.csv --age 45 --sex m \
  --clocks grimage horvath hannum phenoage \
  --sensitivity 40 42 47 49
```

## Input format

A CSV, auto-detected as one of:

- **Long**: two columns — CpG id, then beta value (header names ignored).
  ```
  CpG_site,Beta_value
  cg00000109,0.9238
  cg00000658,0.8628
  ```
- **Matrix**: first column = CpG id, remaining columns = one or more samples.

Beta values are floats in `[0, 1]`.

## Why age & sex are required for GrimAge

GrimAge is a **second-generation, mortality-trained** clock. It estimates DNAm
surrogates of 7 plasma proteins plus smoking pack-years (V2 also adds DNAm A1C &
CRP), then combines them **with chronological age and sex** in a survival model.
Age and sex feed the formula directly, so both are mandatory. (Horvath / Hannum /
PhenoAge don't require them, but passing `--age` lets the tool report
acceleration = clock − chronological age.)

## Interpreting the result — read this

- **Open-source reimplementation, not an official/certified value.** Numbers track
  Horvath's official calculator closely but may differ slightly. For a citable
  number, use the Horvath DNAm Age calculator or a commercial provider.
- **"Acceleration" here = clock − chronological age**, a simple difference. The
  academic *GrimAgeAccel* (residual vs. a same-age cohort) needs a population
  sample and **can't be computed for one person**. So +8 means "epigenetic-
  predicted age is 8 years above chronological age," **not** "8 years older than
  your peers."
- **GrimAge is a risk score in year units**, not "your DNA looks N years old." For
  the "guess my age" question, Horvath/Hannum are the right clocks.
- **Not medical advice.** Research/educational use only.

## Provenance & license

- Skill code: MIT (see [LICENSE](LICENSE)).
- Clock coefficients and the methylation reference are derived from
  **biolearn** (MIT-licensed). See [NOTICE](NOTICE) for full attribution and the
  original clock publications.
- **GrimAge** has commercial-use restrictions (UCLA TDG / the Clock Foundation)
  for cosmetics and life-insurance applications. This repo is a free
  research/educational tool; for commercial licensing contact the Clock Foundation.
