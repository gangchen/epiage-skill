# Coefficient and output-unit audit

Checked 2026-09-13. This records a targeted audit prompted by a report about
DNAmTL and output units. It is not a claim that every bundled predictor has
been independently validated against its original training data.

## DNAmTL: corrected intercept sign

The original [Lu et al. 2019 paper](https://www.aging-us.com/article/102173/text)
and [Supplementary Data 1 workbook](https://cdn.aging-us.com/article/102173/supplementary/SD7/0/aging-v11i16-102173-supplementary-material-SD7.xlsx)
define a 140-CpG linear predictor, in kilobases:

`DNAmTL = 7.924780053 + sum(weight[cpg] * beta[cpg])`.

The workbook sheet `DNAmLTL_addgene_addmQTL` has `Intercept` in A7 and
**+7.924780053** in B7. All 140 CpG weights matched the bundled CSV exactly;
the only discrepancy was the CSV's negative intercept. It is now corrected.
Workbook SHA-256: `6079db3cadca4af5f2c5007e7a11a777f40759398ccf5dd8f94d98981ef2ea10`.
An independently transcribed numeric fixture is in
`tests/fixtures/DNAmTL_Lu2019_original.csv` at the repository root.

With identical beta inputs, the correction raises every finite old result by
**15.849560106 kb**. Recompute prior results; changing the unit, taking absolute
values, or clipping negative outputs would not fix the coefficient error.
The correct unit remains **kb**, with no additional inverse transform.
DNAmTL is a methylation-derived telomere surrogate, not a direct telomere measurement.

Synthetic regression anchors (for checking implementation, not biological examples):

| Input at all 140 CpGs | Correct output, kb |
|---|---:|
| beta = 0 | 7.924780053 |
| beta = 0.5 | 8.417606178 |
| beta = 1 | 8.910432303 |

## Garagnani and Bocklandt: beta values, not years

The bundled `Garagnani.csv` and `Bocklandt.csv` each contain one CpG with
coefficient 1 and no intercept; the script applies an identity transform.
Their outputs are therefore exactly the respective CpG beta values:

- `garagnani`: `cg16867657` (ELOVL2).
- `bocklandt`: `cg09809672` (EDARADD).

These two API keys are retained, but now have unit **beta**, category
**single-CpG marker**, and no age acceleration. The legacy `firstgen` selection
alias still includes them; its name does not change their interpretation.

[Garagnani et al. 2012](https://doi.org/10.1111/acel.12005) studied ELOVL2 as an
age-associated methylation marker. The packaged identity score is not an age
calibration. [Bocklandt et al. 2011](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0014821)
published a saliva age regression using two cytosines; that model is not
implemented by this one-CpG identity score. Do not invent a conversion to years.

## StocZ and Zhang are distinct models

[Tong et al. 2024](https://www.nature.com/articles/s43587-024-00600-8) derives
StocZ from the 514 CpGs of Zhang's chronological-age clock and trains stochastic
predictors of age. **StocZ retains years and age acceleration**; its former
`stochastic (mortality)` category and the skill's non-age description were wrong.

The separate `zhang` key uses the **10-CpG mortality score**, not that 514-CpG
age clock. Its weighted sum is dimensionless and can be negative; the output
unit is now **score**, not an age or probability. Its publication year is 2017:
[Zhang et al., Nature Communications](https://doi.org/10.1038/ncomms14617).

## Other checked unit distinctions

- `dunedinpace` and `dunedinpoam` express pace in **years/year**, not age.
  `DunedinPoAm38` uses 46 CpGs; “38” denotes the training age, not the number
  of CpGs. Its negative intercept is present in the author's model and is not
  analogous to the DNAmTL sign error.
- `epitoc1` computes mean beta across 385 CpGs and retains **score**. It does
  not directly return a number of stem-cell divisions or an age in years.
- Exposome and health outputs remain model scores. A sigmoid output on [0,1]
  is not, by that fact alone, a calibrated probability or a physical measurement.
  This audit does not establish clinical calibration for those predictors.
- The Bahado-Singh Alzheimer's paper was published in 2021, not 2022; that
  metadata has been corrected ([original article](https://doi.org/10.1371/journal.pone.0248375)).

Sources: [DunedinPoAm author code](https://github.com/danbelsky/DunedinPoAm38),
[DunedinPACE paper](https://elifesciences.org/articles/73420),
[EpiTOC paper](https://doi.org/10.1186/s13059-016-1064-3).

## McCartney trait scores: removed unsupported sigmoid transforms

The original [McCartney et al. 2018 paper](https://doi.org/10.1186/s13059-018-1514-1)
and [Additional file 1](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13059-018-1514-1/MediaObjects/13059_2018_1514_MOESM1_ESM.xlsx)
provide 4,307 weights across the eight models below. CpG sets match the bundled
CSVs exactly; differences in coefficient precision are under `5e-8`.
The source table contains no intercept rows and does not specify a sigmoid.
Six transforms inherited from biolearn have been removed: `bmi`, `bodyfat`,
`hdl`, `ldl`, `totalchol`, `education`. `smoking` and `alcohol` already used
linear weighted sums. All eight now report **raw weighted DNAm scores**.

Workbook SHA-256: `b6c9721cc4beb8ae7073f61713e104c6467530f529b785dd44857756d563f884`.
The test fixtures `McCartney2018_original.tsv` and `McCartney2018_expected.json`
are transcribed from the original table, not generated from the production CSVs.

| API key | CpGs | Original-table score at beta = 0.5 |
|---|---:|---:|
| bmi | 1109 | -0.7326852944 |
| smoking | 233 | 2.0069389860 |
| alcohol | 450 | -5.5593637390 |
| education | 373 | -2.4688491745 |
| totalchol | 204 | -1.0609031930 |
| hdl | 737 | 2.5504501939 |
| ldl | 233 | -2.7555480414 |
| bodyfat | 968 | -15.3519314083 |

These are synthetic regression anchors, not typical biological values.
Recompute earlier six-model outputs: removing sigmoid changes their numeric scale.
The paper and [supplementary methods](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13059-018-1514-1/MediaObjects/13059_2018_1514_MOESM4_ESM.pdf)
describe log transforms of BMI, smoking pack-years + 1 and alcohol units/week + 1,
an education category code, and adjustment for age, sex and genetic PCs before
training. Therefore, these weighted sums must not be relabeled as kg/m²,
cholesterol concentration, years of education, pack-years or disease probability.
Simply exponentiating a score does not restore the original phenotype. This
change restores the published weighted-score calculation, not absolute clinical
calibration or all original preprocessing stages.

## Temporarily unavailable: CVD and depression

Both keys remain registered for compatibility, but the CLI returns a null value,
`status=unavailable`, and a `status_reason`. Other requested models continue.
Direct `predict_linear` calls also refuse these specifications. Historical files
are retained for provenance and future reconstruction; their presence does not
mean the old formula is supported. No substitute probability or raw score is
silently returned. Thus 35 of 37 registered entries currently support calculation.

### CVD / Westerman

The [author's scoring code and parameters](https://github.com/kwesterman/meth_cvd/tree/master/output/mrs_calculation)
construct raw within-cohort beta-weighted scores, standardize them, and combine
four cohorts for CSL. The bundled 235 weights match only a subset of the WHI
model, multiplied by 0.4764 (maximum rounding difference under `5e-10`). The
author's WHI model contains 320 nonzero CpGs and its CSL weight is 0.47643454.
The bundled representation omits 85 WHI terms, the other three cohort scores,
and their standardization; it does not match the author's complete combined
342-CpG or ComBat 817-CpG models either. Its sigmoid is also unsupported.
Removing that sigmoid alone would leave an incomplete model, so numerical
output is paused until a complete scoring specification and reference test exist.

### Depression / Barbu

The [original paper](https://doi.org/10.1038/s41380-020-0808-3) specifies normalized
**M-value** inputs, including preprocessing and covariate adjustment; the previous
implementation used beta values directly. All 196 weights match
[Supplementary file 2A](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41380-020-0808-3/MediaObjects/41380_2020_808_MOESM4_ESM.xlsx)
to within `5e-10`, but that table has no intercept. The extra 12.2169841 in the
bundled CSV is unsupported by the source and is numerically identical to the
Lin age-model intercept; that equality does not establish how the error arose.
The author's supplementary analysis code starts from precomputed scores and
does not resolve their full construction or calibration. Neither switching to
M values nor setting the intercept to zero alone establishes equivalence, so
this numerical output is also paused pending source-model reconstruction.

Depression source-table SHA-256:
`f66ab5ce6e0f9e8958c1f57ff277a07521f9dc24697bad9bafc8331a25f79b88`.

## Validation and interpretation

Regression tests compare DNAmTL and McCartney coefficient sets to original-source
fixtures, exercise known synthetic outputs, check CLI units and acceleration
behavior, and ensure blocked models return no score while other models continue.
Runtime remains offline and requires only
the existing Python dependencies. Agreement with biolearn is useful for checking
an implementation, but does not establish that upstream coefficients or
metadata match the original paper; StocZ illustrates that distinction.
