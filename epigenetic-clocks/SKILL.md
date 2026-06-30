---
name: epigenetic-clocks
description: >-
  Compute epigenetic / DNA-methylation aging clocks — GrimAge (V1 & V2), Horvath
  (v1 & skin-blood), Hannum, PhenoAge, the Ying causality clocks, DunedinPoAm pace
  of aging, DNAmTL telomere length, and ~24 clocks total — from a DNA methylation
  beta-value file. Use this whenever the user provides methylation array data (a
  CSV/TSV of CpG sites with beta values, e.g. from an Illumina EPIC/450K array) and
  wants their biological age, epigenetic age, GrimAge, DNAm age, age acceleration,
  pace of aging, "甲基化年龄", "生物年龄", "表观遗传时钟", "衰老时钟", or how old
  their DNA "looks" — even if they just drop a methylation file and ask "how old am
  I biologically". Also triggers on "GrimAge", "Horvath clock", "PhenoAge",
  "DunedinPACE/PoAm", "epigenetic clock", or comparing several aging clocks on one
  sample. Designed for human WHOLE BLOOD samples (e.g. WeGene/EPIC/450K/MSA blood
  exports); missing CpGs are imputed by default with methyLImp using a whole-blood
  reference panel.
---

# Epigenetic Clock Calculator

Computes DNA-methylation aging clocks from a beta-value file. **25 clocks** are
available; GrimAge is the headline one.

**Self-contained.** Needs only `pandas` + `numpy` — no biolearn, torch, seaborn,
scipy, or network. Clock coefficients and references are vendored under `data/`
(extracted from the open-source biolearn library, trimmed to the CpGs the clocks
use, ~1 MB total). `scripts/compute_clocks.py` faithfully reimplements biolearn's
`GrimageModel`, `LinearMethylationModel`, and the DunedinPACE quantile-
normalization (with a numpy-only `rankdata`), reproducing biolearn's outputs for
all 25 clocks (verified to <0.005, i.e. rounding only).

## The clocks (run `--list-clocks` for the live list)

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

Group aliases for `--clocks`: `all`, `grimage`, `core` (default), `firstgen`,
`secondgen`, `thirdgen`.

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

### 1. Inspect the input file
Auto-detected layouts: **Long** (two columns: CpG id + beta) or **Matrix** (first
column CpG id, remaining columns = samples). Betas are floats in [0,1]. Note the
row count; coverage of each clock's CpGs is reported per clock.

### 2. Get age and sex
Required for GrimAge. If only an age range is known, use the midpoint and run
`--sensitivity` so the user sees how much the answer depends on the exact age.

### 3. Check dependencies (only pandas + numpy)
```bash
python3 -c "import pandas, numpy" 2>/dev/null && echo OK || pip install pandas numpy
```
No venv, biolearn, torch, or network needed.

### 4. Run
```bash
python3 <skill-dir>/scripts/compute_clocks.py \
  --input "<betas.csv>" --age <years> --sex <m|f> \
  --clocks all                 # or: core (default), grimage, firstgen, secondgen, or specific keys
  # --sensitivity 40 42 47 49  # optional, when exact age is uncertain
python3 <skill-dir>/scripts/compute_clocks.py --list-clocks   # see all keys
```
The script imputes any missing clock CpGs from the bundled (trimmed)
`data/sesame_450k_median.csv` population reference, runs each model, and prints a
table of value, acceleration (for year-unit clocks), and coverage, plus a JSON
line for downstream use.

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
    rebuild/customize, run `scripts/build_blood_panel.py`. If the
    file is ever missing the tool degrades gracefully to global-median and says so.
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
  - Trait & disease predictors (BMI, cholesterol, smoking, alcohol, Alzheimer's,
    CVD, …) — these are biomarker models, not aging clocks.
  If a user specifically needs one of these, tell them it requires the full
  biolearn install.
- **Provenance.** Coefficients are biolearn's (MIT), which reimplements the
  published clocks. GrimAge has commercial-use restrictions (UCLA TDG / Clock
  Foundation) for cosmetics & life-insurance use.
