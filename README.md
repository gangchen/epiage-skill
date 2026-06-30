# epiage-skill

An installable agent skill for computing **epigenetic / DNA-methylation aging
clocks** from a CpG beta-value file — entirely offline, with only `pandas` +
`numpy`.

```bash
npx skills add gangchen/epiage-skill
```

## At a glance

- **25 aging clocks** in one run — GrimAge V1/V2, Horvath ×2, Hannum, PhenoAge,
  Ying causality clocks, DunedinPACE/PoAm, DNAmTL, and more.
- **For human whole blood** — give it a blood methylation export (WeGene / EPIC /
  450K / MSA) plus age + sex, get every clock with acceleration and a per-clock
  reliability flag.
- **Self-contained & offline** — only `pandas` + `numpy`. Coefficients, the
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

Run `--list-clocks` for the full list and selectable keys.

- **Self-contained**: only `pandas` + `numpy`. No `biolearn`, `torch`, `scipy`, or
  network. Coefficients + references are vendored under `epigenetic-clocks/data/`
  (~1 MB; includes DunedinPACE's 20k-probe normalization reference).
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

A CSV, auto-detected as one of:

- **Long**: two columns — CpG id, then beta value (header names ignored).
  ```
  CpG_site,Beta_value
  cg00000109,0.9238
  cg00000658,0.8628
  ```
- **Matrix**: first column = CpG id, remaining columns = one or more samples.

Beta values are floats in `[0, 1]`.

### Where to get the raw data

You need a per-CpG **methylation beta-value** file. One consumer source is
**[WeGene (微基因)](https://www.wegene.com/)**, whose methylation product lets you
download your raw beta values as a CpG-vs-beta CSV — exactly the **Long** format
above (`CpG_site,Beta_value`). Export it from your WeGene account and pass it
straight to `--input`.

Any platform that outputs Illumina EPIC/450K beta values works too (e.g. an
`idat`-derived matrix processed with `minfi`/`sesame`). Note that coverage varies
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

If the panel is ever missing the tool falls back to global-median and prints the
active mode. Note: methyLImp mainly helps the heavily-imputed clocks; high-coverage
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
clocks (cord blood / newborns), and trait/disease predictors (BMI, cholesterol,
smoking, Alzheimer's, …) — the last are biomarker models, not aging clocks. These
require the full biolearn install.

## Provenance & license

- Skill code: MIT (see [LICENSE](LICENSE)).
- Clock coefficients and the methylation reference are derived from **biolearn**
  (MIT). See [NOTICE](NOTICE) for full attribution and the original clock papers.
- **GrimAge** has commercial-use restrictions (UCLA TDG / the Clock Foundation) for
  cosmetics and life-insurance applications. This repo is a free
  research/educational tool; for commercial licensing contact the Clock Foundation.
