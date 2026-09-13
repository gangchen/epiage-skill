---
name: epigenetic-clocks
description: >-
  Process human whole-blood methylation array IDAT pairs with SeSAMe (sesame),
  including QC and beta export for 450K, EPIC, EPICv2 and MSA when supported by
  the installed package/data. Compute 25 epigenetic aging clocks and 12 exposome
  or health scores from IDAT-derived or existing CSV/TSV beta matrices. Use for
  IDAT preprocessing, biological age, GrimAge, Horvath, PhenoAge, DunedinPACE,
  age acceleration, methylation lifestyle scores, 甲基化芯片、甲基化年龄、
  生物年龄、表观遗传时钟 or 暴露组. Missing clock CpGs use a whole-blood
  imputation panel; interpretation is restricted to human whole blood.
metadata:
  homepage: https://github.com/gangchen/epiage-skill
  openclaw:
    requires:
      bins: [python3]
---

# Epigenetic Clock Calculator

Processes raw IDAT pairs with SeSAMe and computes DNA-methylation aging clocks
from the resulting beta matrix or an existing beta-value file. **25 clocks** and
**12 exposome/health predictors** are available; GrimAge is the headline clock.

**Beta-to-clock data is bundled.** Once `pandas` + `numpy` are installed, calculation
needs no biolearn, torch, seaborn, scipy, or network. Clock coefficients,
normalization references and the blood imputation panel are vendored under
`data/` (~6 MB total). `scripts/compute_clocks.py` faithfully reimplements biolearn's
`GrimageModel`, `LinearMethylationModel`, and the DunedinPACE quantile-
normalization (with a numpy-only `rankdata`), reproducing biolearn's outputs for
all 37 models — 25 aging clocks + 12 exposome/health predictors — verified bit-exact
against biolearn on complete beta input. Masked matrix values are handled per
sample; DunedinPACE uses each sample's observed background plus gold means.

**The complete IDAT-to-clock workflow can run offline.** IDAT preprocessing
requires R plus `sesame` and `sesameData`. Public IDAT references are bundled
(47.3 MiB compressed), matched to `sesame` 1.24.0 / `sesameData` 1.24.0 and tested
with R 4.4.3. R/Python environments are not packaged: install software separately
when needed. At installation and startup, check R, Python, their packages and
required system libraries; list every detected
missing or incompatible item, its purpose, version and size (unknown if not yet
known), then ask whether the user permits downloading and installing. Do not install packages,
initialize an online cache, or rebuild references without that consent. If the
machine has no internet, ask the user to transfer compatible standard
installers/packages or missing reference files from another machine.
Sample data stays local and must not enter Git. For IDAT input or offline setup,
read [references/idat-sesame.md](references/idat-sesame.md) and use
`scripts/preprocess_idat.R`. Existing beta files bypass IDAT preprocessing.

## The models (run `--list-clocks` for the live list)

**Aging clocks (25):**

| family | clocks | unit |
|---|---|---|
| **GrimAge** (2nd-gen, mortality) | `grimagev1`, `grimagev2` | years *(needs age+sex)* |
| **1st-gen chronological** | `horvath`, `horvath2`, `hannum`, `lin`, `vidalbralo`, `weidner`, `garagnani`, `bocklandt` | years |
| **2nd-gen biological age** | `phenoage`, `hrsinchphenoage` | years |
| **Ying 2022 (causality)** | `yingcausage`, `yingdamage`, `yingadaptage` | years |
| **stochastic** | `stoch`, `stocp`, `stocz` | years |
| **tissue-specific** | `pedbe` (pediatric buccal), `cortical` (brain) | years |
| **3rd-gen pace of aging** | `dunedinpace`, `dunedinpoam` | years/year |
| **other markers** | `dnamtl` (telomere kb), `zhang` (mortality), `epitoc1` (mitotic) | varies |

**Exposome / lifestyle & health predictors (12)** — methylation *scores*, not aging
clocks (McCartney 2018 / Reed / disease EWAS):

| group | models | unit |
|---|---|---|
| **exposome** | `smoking`, `alcohol`, `bmi`, `bmi_reed`, `bodyfat`, `hdl`, `ldl`, `totalchol`, `education` | score |
| **health** | `cvd` (coronary heart disease), `alzheimers`, `depression` | risk |

Group aliases for `--clocks`: `all` (everything), `aging` (25 clocks), `grimage`,
`core` (default), `firstgen`, `secondgen`, `thirdgen`, `exposome`, `health`,
`phenotypes` (exposome+health).

**`dunedinpace` needs ~20k background CpGs** for its quantile normalization (not
just its 173 model CpGs). On a sparse input it self-imputes the missing background
from the gold-standard reference; if its reported coverage is well below ~90% the
result is dominated by the reference and unreliable — say so. A full EPIC/450K
export covers it fine.

## Key facts about GrimAge (so you interpret it correctly)

GrimAge is **second-generation, mortality-trained**: it estimates DNAm surrogates
of 7 plasma proteins + smoking pack-years (V2 adds DNAm A1C & CRP), then combines
them **with chronological age and sex** in a survival model. Age and sex feed the
formula directly, so **both are required** for any `grimage*` clock. The other
clocks don't need them, but passing `--age` lets the tool report acceleration.

## Workflow

### 1. Identify and inspect the input
For IDATs (a directory, a prefix, paired `_Grn.idat` / `_Red.idat`, optionally
gzip, or a sample sheet), follow [the SeSAMe workflow](references/idat-sesame.md).
After checking local dependencies, use `--inspect` to validate IDAT input and
identify the human methylation array from its address signature; do not infer a
platform from filenames or probe count, or force an unsupported array manifest.
Then check the matching local annotation and run SeSAMe QCDPB to produce beta
values and QC. Inspect `qc.csv` before imputation or clocks. Preserve masks as
`NA`; do not replace them with zero or remove a locus from every sample because
it failed in one sample. High sample missingness requires review. MSA has no
recommended design mask in the bundled 1.24.0 resources: retain and report
`design_mask_available=FALSE`; never claim design-mask QC was completed.
Preprocessing alone does not need age or sex.

For existing beta files, do not rerun IDAT preprocessing.
Auto-detected layouts: **Long** (two columns: CpG id + beta) or **Matrix** (first
column CpG id, remaining columns = samples); CSV/TSV and gzip are supported.
Betas are floats in [0,1] or explicit `NA`/blank. Duplicate CpG IDs are rejected:
resolve array probe replicates in SeSAMe before calculating clocks. Coverage is
the count of observed, unmasked features for each sample and clock.
For beta-only data, retain the platform as unknown unless provenance establishes
it; CpG count cannot prove a chip model.

### 2. Get age and sex when calculating GrimAge
These are required only for GrimAge. If only an age range is known, use the midpoint and run
`--sensitivity` so the user sees how much the answer depends on the exact age.
The CLI's age and sex apply to every column in a matrix. For multiple people
with different metadata, run separate sample matrices with each person's known
age/sex; the IDAT sample sheet does not supply these values to the clock script.

### 3. Check local clock dependencies and references
```bash
python3 -c "import pandas, numpy; print('Python dependencies available')"
python3 <skill-dir>/scripts/compute_clocks.py --check-resources --clocks all
```
`--check-resources` lists all missing files required by the selected clocks,
including `blood_panel.npz` for methyLImp, and exits without computing scores.
Report missing resources and ask about downloading or transferring them as
above. Restore missing files before calculation. A complete local setup needs
no network.

### 4. Impute and compute clocks
```bash
python3 <skill-dir>/scripts/compute_clocks.py \
  --input "<betas.csv>" --age <years> --sex <m|f> \
  --clocks all                 # or: core (default), grimage, firstgen, secondgen, or specific keys
  # --sensitivity 40 42 47 49  # optional, when exact age is uncertain
python3 <skill-dir>/scripts/compute_clocks.py --list-clocks   # see all keys
```
The script imputes absent and sample-specific masked clock CpGs with the bundled
whole-blood panel (with median fallback for CpGs the panel lacks), runs each model, and prints a
table of value, acceleration (for year-unit clocks), and coverage, plus a JSON
line for downstream use.
`n_missing` counts masked/absent required features, `n_imputed` counts successful
fills, and `n_unresolved` counts missing features without a usable reference.
An unresolved model gets `status=unavailable` and a null JSON value; other models
continue. Do not interpret an unavailable score or call an imputed value measured.

### 5. Report and interpret
Lead with GrimAge, then the comparison. Always convey these caveats:

- **Open-source reimplementation, not the official calculator.** Values track
  Horvath's official server / Clock Foundation closely but aren't a certified
  number. For a citable value, point to the Horvath DNAm Age calculator or a
  commercial provider.
- **"Acceleration" here = clock − chronological age**, a simple difference. The
  academic *AgeAccel* (residual vs. a same-age cohort) needs a population sample
  and **can't be computed for one person** — so +8 means "epigenetic-predicted age
  is 8 years above chronological age," not "8 years older than peers."
- **Generations differ.** 1st-gen (Horvath/Hannum/Lin/…) target chronological age
  and tend to land near the true age. 2nd-gen (GrimAge/PhenoAge) target health
  outcomes and predict mortality/aging better — they can diverge from true age by
  design. Don't read GrimAge as "looks N years old"; it's a risk score in years.
- **Non-year clocks** (`dunedinpoam` pace, `dnamtl` telomere kb, `zhang`/`stocz`
  mortality risk, `epitoc1` mitotic) are not ages — no acceleration is shown.
- **Exposome / health predictors are methylation *scores*, not clinical values.**
  `bmi`/`hdl`/`smoking`/`cvd`/etc. output a DNAm-predicted score/risk (many via a
  sigmoid, so in [0,1]) — NOT your actual BMI, cholesterol, or a diagnosis. Read
  them as relative epigenetic signals, and only where coverage is high — several
  (`education`, `ldl`, `cvd`, `alzheimers`, `depression`) have low coverage on
  sparse MSA data, so distrust those. Never present them as a medical result.
- **Coverage / imputation**: a few missing clock CpGs filled from a population
  median is normal; flag clocks whose coverage drops below ~90%. The output's
  `n_lowconf` column counts imputed CpGs whose blood-reference SD > 0.08 (fills to
  distrust) — a per-CpG signal sharper than coverage alone. A clock with many
  `n_lowconf` fills is unreliable on this sample even if it "ran".
- **Blood-only design.** This skill assumes human WHOLE BLOOD. Clocks are applied
  on blood; the imputation reference is blood. Do not use it on other tissues.
- **Imputation of missing CpGs — methyLImp by default.** Missing clock CpGs (e.g.
  EPIC-trained CpGs absent on an MSA/WeGene export) are imputed with **methyLImp**:
  reduced-rank (PCA) regression that predicts each missing CpG from the sample's
  *observed* CpGs using the inter-CpG correlation structure of a whole-blood panel
  — more accurate than a flat median (~30% lower RMSE when tissue matches). CpGs
  the panel lacks fall back to the blood median; `n_lowconf` (blood SD>0.08) flags
  unreliable fills. The active mode is printed at run start.
  - The panel lives in `data/blood_panel.npz` and **is bundled** (~5 MB; prebuilt
    from GSE40279 = 656 whole-blood 450K samples, reduced to the clock CpGs + 50
    blood PCs + per-CpG median/SD). So methyLImp is active out of the box — it
    cuts imputation RMSE ~10% vs a flat median on a held-out-CpG benchmark. To
    rebuild/customize, run `scripts/build_blood_panel.py` only after permission
    to download its public source dataset. If the panel is missing, report the
    missing resource and ask about downloading/transferring it before proceeding.
    The CLI requires the panel and does not silently switch to global median.
  - Reality check: methyLImp mainly improves the *heavily-imputed* clocks. High-
    coverage clocks (GrimAge/Horvath/PhenoAge) impute few CpGs, so their values
    move <0.2 yr regardless — the confidence flag is the bigger practical win.
- **Not medical advice** — research/educational use only.

## Notes
- Matrix input → one result row per sample per clock.
- **Deliberately excluded** (can't be reproduced faithfully without heavier
  machinery, or aren't aging clocks):
  - PC-clocks (`PCHorvath`…), `AltumAge`, `GPAge` — need PCA rotation / neural nets.
  - Gestational clocks (Knight, Lee, Mayne, Bohlin) — for cord blood / newborns.
  If a user specifically needs one of these, tell them it requires the full
  biolearn install.
- **Provenance.** Coefficients are biolearn's (MIT), which reimplements the
  published clocks. GrimAge has commercial-use restrictions (UCLA TDG / Clock
  Foundation) for cosmetics & life-insurance use.
