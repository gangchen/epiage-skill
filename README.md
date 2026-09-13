# epiage-skill

An installable agent skill for computing **epigenetic / DNA-methylation aging
clocks** from a CpG beta-value file or raw methylation-array **IDAT pairs**.
Beta-to-clock calculation needs only `pandas` + `numpy`; IDAT preprocessing uses
R/Bioconductor **SeSAMe**. The necessary public array references are bundled
(47.3 MiB compressed). With compatible dependencies installed, the complete
workflow runs offline:

**Identify input and chip model → check local annotations → SeSAMe beta values
and QC → imputation → epigenetic clocks.**

```bash
npx skills add gangchen/epiage-skill
```

> ⚠️ **Research / educational use only — not a medical device.** This computes
> DNA-methylation *research scores* (biological age, pace of aging, and relative
> disease-risk / lifestyle scores). It does **not** diagnose, treat, or provide
> medical advice, and its outputs are not clinical measurements. Do not make health
> decisions from them; consult a qualified clinician.
>
> 🔒 **Private by design.** Sample processing runs locally; methylation data is
> not uploaded. Installation and startup check local resources first. Missing
> items are listed with their purpose and size (or an explicit unknown size),
> then the agent asks whether downloading is allowed. Nothing is downloaded
> without consent. On a machine without internet, install the missing software
> from standard installation packages transferred from another machine.

## At a glance

- **37 methylation models** in one run — 25 aging clocks (GrimAge V1/V2, Horvath ×2,
  Hannum, PhenoAge, Ying causality clocks, DunedinPACE/PoAm, DNAmTL, …) plus 12
  exposome & health predictors (DNAm smoking, alcohol, BMI, body fat, cholesterol,
  education, and CHD / Alzheimer's / depression risk scores).
- **For human whole blood** — give it a blood methylation export (WeGene / EPIC /
  450K / MSA), get clock results and per-clock coverage. Age + sex are required
  only for GrimAge; age also enables acceleration for clocks measured in years.
- **IDAT support** — paired `.idat`/`.idat.gz`, recursive directories and sample
  sheets; validates IDAT input, identifies HM450, EPIC, EPICv2 or MSA from array
  address signatures, then runs SeSAMe QCDPB, QC and masked beta export with
  replicate-probe collapse. The bundled SeSAMe 1.24.0 references include no
  recommended MSA design mask, so MSA results explicitly flag that QC limitation.
- **Self-contained beta-to-clock calculation** — only `pandas` + `numpy`. Coefficients, the
  DunedinPACE normalization reference, and a whole-blood methyLImp panel are all
  vendored (~6 MB). No biolearn / torch / scipy / network at runtime.
- **Faithful** — reimplements biolearn's clocks, verified to match to <0.005.
- **Smart imputation** — missing CpGs (common on the newer MSA chip) are filled by
  default with **methyLImp** (correlation-based, ~10% lower error than a median on
  a held-out blood benchmark), with a confidence flag on each fill.

Then just hand your agent a methylation file and ask for your biological age.

## Skill: `epigenetic-clocks`

> **For human whole-blood samples.** The clocks and the imputation reference are
> blood-based — designed for blood methylation exports (WeGene / EPIC / 450K /
> MSA). Don't use it on other tissues.

Computes **25 aging clocks** from a methylation beta-value CSV (e.g. an Illumina
EPIC / 450K array export), including:

- **GrimAge** V1 & V2 (2nd-gen, mortality-trained)
- **1st-gen chronological**: Horvath (v1 & skin-blood), Hannum, Lin, Vidal-Bralo,
  Weidner, Garagnani, Bocklandt
- **2nd-gen biological age**: PhenoAge, HRSInCH-PhenoAge
- **Ying 2022 causality clocks**: CausAge, DamAge, AdaptAge
- **Stochastic clocks**: StocH, StocP, StocZ
- **Tissue-specific**: PEDBE (pediatric buccal), Cortical (brain)
- **3rd-gen pace of aging**: DunedinPACE, DunedinPoAm
- **Other markers**: DNAmTL (telomere length), Zhang (mortality), EpiTOC1 (mitotic)

Plus **12 exposome & health predictors** (methylation *scores*, not aging clocks):

- **exposome / lifestyle** (McCartney 2018 / Reed): smoking, alcohol, BMI (×2), body
  fat, HDL / LDL / total cholesterol, education
- **health / disease risk**: coronary heart disease, Alzheimer's, depression

Run `--list-clocks` for the full list. Group aliases for `--clocks`: `all`, `aging`,
`core` (default), `grimage`, `firstgen`, `secondgen`, `thirdgen`, `exposome`,
`health`, `phenotypes`. These predictors are relative DNAm scores (many sigmoid-
squashed to [0,1]) — **not** your actual BMI/cholesterol or a diagnosis.

- **Bundled clock data**: after installing `pandas` + `numpy`, clock calculation
  requires no `biolearn`, `torch`, `scipy`, or network. Coefficients, DunedinPACE's
  20k-probe normalization reference and the blood imputation panel are vendored
  under `epigenetic-clocks/data/` (~6 MB, in addition to the IDAT reference bundle).
- **Faithful**: the math reimplements [biolearn](https://bio-learn.github.io/)'s
  `GrimageModel`, `LinearMethylationModel`, and the DunedinPACE quantile
  normalization (with a numpy-only `rankdata`), verified to reproduce biolearn's
  outputs for all 25 clocks (agreement < 0.005, i.e. rounding only).

## Install

```bash
npx skills add gangchen/epiage-skill
```

Once installed, give your agent a methylation file and ask for your GrimAge /
biological age. You can also run the script directly:

```bash
python3 epigenetic-clocks/scripts/compute_clocks.py \
  --input betas.csv --age 45 --sex m \
  --clocks all                 # or: core (default), grimage, firstgen, secondgen, or specific keys
  # --sensitivity 40 42 47 49  # optional, when exact age is uncertain

python3 epigenetic-clocks/scripts/compute_clocks.py --list-clocks
```

## Input format

### Raw IDAT (SeSAMe)

The bundled references match `sesame` **1.24.0** + `sesameData` **1.24.0**, tested
with R **4.4.3**. R, Python and their dependencies are installed separately;
the skill does not package their environments. Follow the
[IDAT/offline setup guide](epigenetic-clocks/references/idat-sesame.md) if anything
is missing, and ask before downloading and installing it. With dependencies ready:

```bash
Rscript epigenetic-clocks/scripts/preprocess_idat.R --check-deps
python3 epigenetic-clocks/scripts/compute_clocks.py --check-resources --clocks all

# Identify the actual array before processing; no network is used.
Rscript epigenetic-clocks/scripts/preprocess_idat.R \
  --input /path/to/idats --inspect

Rscript epigenetic-clocks/scripts/preprocess_idat.R \
  --input /path/to/idats --output-dir idat-output

# Review idat-output/qc.csv before interpreting clocks.
# For a single person with known age and sex:
python3 epigenetic-clocks/scripts/compute_clocks.py \
  --input idat-output/betas.csv --age 45 --sex m --clocks all
```

The preprocessor also accepts a single IDAT prefix/channel or `--sample-sheet`
with `sample_id,idat_prefix`. Outputs include `input_inspection.csv`, `betas.csv`,
`qc.csv`, sample/file checksums and software versions. Review QC before
imputation and clock interpretation; imputation does not repair a failed sample.
Existing output directories are never replaced. The age/sex arguments apply to
all samples in a matrix: use separate sample runs when people have different
metadata. IDAT preprocessing itself needs neither. Keep IDATs, sample metadata
and results outside Git.

The default run extracts the bundled references to a temporary local cache.
Other SeSAMe versions need a matching cache supplied with `--cache`; preparing
it online requires prior download consent. At startup, check R, Python, their
packages and required system libraries. If the machine cannot connect, the user
can transfer compatible standard installers/packages and the missing reference
files from another machine.

### Beta files

A CSV/TSV (optionally `.gz`), auto-detected as one of:

- **Long**: two columns — CpG id, then beta value (header names ignored).
  ```
  CpG_site,Beta_value
  cg00000109,0.9238
  cg00000658,0.8628
  ```
- **Matrix**: first column = CpG id, remaining columns = one or more samples.

Beta values are floats in `[0, 1]`; masked values remain `NA` or blank. Named
sample columns are preserved (`Beta_value` remains the legacy `Sample` label).
Duplicate CpG rows are rejected; collapse array replicates during preprocessing.
Coverage and imputation counts are per sample, including masked values. Models
whose missing features cannot be filled return `status=unavailable` with a null
JSON value, while other models continue.
An existing beta file can establish that it contains methylation measurements,
but CpG counts alone cannot prove a specific chip model. Record the platform as
unknown unless an array manifest, export metadata or the original IDATs support it.

### Where to get the raw data

You need a per-CpG **methylation beta-value** file. One consumer source is
**[WeGene (微基因)](https://www.wegene.com/)**, whose methylation product lets you
download your raw beta values as a CpG-vs-beta CSV — exactly the **Long** format
above (`CpG_site,Beta_value`). Export it from your WeGene account and pass it
straight to `--input`.

Existing Illumina EPIC/450K beta values work too; raw IDATs can now be processed
with the included SeSAMe script. Note that coverage varies
by source: clocks needing many probes — especially `dunedinpace` (~20k background
CpGs) — are only reliable on a fairly complete export; the tool reports per-clock
coverage so you can tell.

### Imputing missing CpGs (methyLImp, blood)

Newer arrays like the **MSA** chip drop many EPIC-trained clock CpGs. By default the
tool imputes them with **methyLImp** — reduced-rank (PCA) regression that predicts a
missing CpG from your *observed* CpGs using the inter-CpG correlation structure of a
whole-blood reference (~30% lower error than a flat median when tissue matches). The
`n_lowconf` column flags fills whose blood SD > 0.08 (distrust those).

The blood panel (`data/blood_panel.npz`, ~5 MB) **is bundled** — prebuilt from
GSE40279 (656 whole-blood 450K samples), reduced to the clock CpGs + 50 blood PCs +
per-CpG median/SD. So methyLImp is active out of the box; on a held-out-CpG benchmark
it cuts imputation RMSE ~10% vs a flat median. To rebuild/customize:

```bash
python3 epigenetic-clocks/scripts/build_blood_panel.py   # GSE40279, 656 blood samples
```

Rebuilding downloads the public source dataset and is optional; ask for download
consent first. Normal calculation uses the bundled panel without downloading it.

If the panel or another required reference file is missing, the CLI lists all
missing files and stops. Restore them before calculation; the agent asks before
downloading. CpGs absent from an available panel can still use median fallback.
Note: methyLImp mainly helps the heavily-imputed clocks; high-coverage
clocks (GrimAge/Horvath/PhenoAge) move <0.2 yr either way.

## Why age & sex are required for GrimAge

GrimAge is a **2nd-generation, mortality-trained** clock: it estimates DNAm
surrogates of 7 plasma proteins + smoking pack-years (V2 also adds DNAm A1C & CRP),
then combines them **with chronological age and sex** in a survival model. Both are
mandatory for the `grimage*` clocks. The other clocks don't need them, but passing
`--age` lets the tool report acceleration (= clock − chronological age) for the
year-unit clocks.

## Interpreting the result — read this

- **Open-source reimplementation, not an official/certified value.** Numbers track
  Horvath's official calculator closely but may differ slightly. For a citable
  number, use the Horvath DNAm Age calculator or a commercial provider.
- **"Acceleration" = clock − chronological age**, a simple difference. The academic
  *AgeAccel* (residual vs. a same-age cohort) needs a population sample and **can't
  be computed for one person**. So +8 means "epigenetic-predicted age is 8 years
  above chronological age," **not** "8 years older than your peers."
- **Generations differ.** 1st-gen clocks target chronological age and land near
  your true age; 2nd-gen (GrimAge/PhenoAge) target health outcomes and predict
  mortality better — they can diverge from true age by design.
- **Coverage matters.** Each clock reports CpG coverage; clocks heavily imputed on
  a sparse input (coverage < ~90%) are less reliable for that sample.
- **Non-year clocks** (pace, telomere kb, mortality risk, mitotic) are not ages.
- **Not medical advice.** Research/educational use only.

> **DunedinPACE note:** its quantile normalization needs ~20k background CpGs (not
> just its 173 model CpGs). On a sparse input it self-imputes the rest from the
> gold-standard reference; if the reported coverage is well below ~90% the result
> is unreliable. A full EPIC/450K export covers it fine.

## Deliberately not included

PC-clocks / `AltumAge` / `GPAge` (need PCA rotation or neural nets), gestational
clocks (cord blood / newborns). These require additional implementations or the
full biolearn install.

## Provenance & license

- Skill code: MIT (see [LICENSE](LICENSE)).
- Clock coefficients and the methylation reference are derived from **biolearn**
  (MIT). See [NOTICE](NOTICE) for full attribution and the original clock papers.
- Bundled array references come from **sesameData 1.24.0 / ExperimentHub**
  (sesameData: Artistic-2.0). The [resource inventory](epigenetic-clocks/data/sesame-resources.json)
  records accessions, sizes and SHA-256 hashes; see [NOTICE](NOTICE) and the
  [included license](epigenetic-clocks/data/SESAME_DATA_LICENSE.txt).
- **GrimAge** has commercial-use restrictions (UCLA TDG / the Clock Foundation) for
  cosmetics and life-insurance applications. This repo is a free
  research/educational tool; for commercial licensing contact the Clock Foundation.

## Validation

```bash
python3 -m unittest discover -s tests -v
Rscript --vanilla tests/test_preprocess_idat.R
```

The IDAT workflow has been exercised on an actual EPICv2 sample with network
access restricted. Its betas matched those generated using the original
reference cache. All 37 models retained their previous numerical results on
complete synthetic input. Tests also cover sample-specific missingness,
DunedinPACE sample independence, IDAT pairing and missing-resource reporting.
