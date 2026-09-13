# IDAT → SeSAMe QC/betas → imputation → epigenetic clocks

Use this workflow for human whole-blood methylation arrays. A complete local
installation runs every stage offline. The skill contains the public array and
clock references. R, Python, their packages and system libraries are installed
separately; this skill does not package or export software environments.

## Check resources before downloading

At installation and each startup, inspect what is already available. Classify
the input from its files first, then obtain any missing reader dependencies
before validating IDAT content. For IDAT processing, check R, `sesame`,
`sesameData`, their dependencies, required system libraries and the reference
cache. For clock calculation,
check Python, `numpy`, `pandas`, the selected model coefficient files,
`DunedinPACE_Gold_Means.csv` when needed, and `blood_panel.npz` for methyLImp.
After confirming Python dependencies, list any missing clock/reference files:

```bash
python3 epigenetic-clocks/scripts/compute_clocks.py --check-resources --clocks all
```

When anything is missing or incompatible, present one list with **item, purpose,
required version, size and local status**. Use the sizes below or the manifest;
for a package whose download size is not known, say **unknown**, not zero.
Ask whether the user permits downloading and installing the listed items before making network
requests or running installation/cache initialization commands. For example:

> 缺少 sesame 1.24.0 及其 R 依赖（用于读取 IDAT、计算 beta 和 QC；下载大小待确认）。
> 本地已有 47.3 MiB 的芯片注释包。是否允许下载并安装缺少的依赖？如果这台机器不能联网，
> 可以由你在另一台机器下载兼容的标准安装包，再传入本机安装。

Existing authorization for those downloads remains valid. Without consent,
keep processing offline and report the blocked stage. A machine with no network
cannot download even after consent: the user can download compatible standard
installers/packages or missing reference files on another machine and transfer
them. Do not substitute another chip's manifest, skip QC,
or silently call median fallback methyLImp to work around missing resources.

## Bundled references and compatibility

`data/sesame-reference-cache.tar.gz` is **49,559,929 bytes (47.3 MiB)**. It contains
eight SeSAMe resources plus the ExperimentHub metadata needed for local lookup:

| Resource | Purpose | Stored size before the outer archive |
|---|---|---:|
| `idatSignature` | Identify the array from IDAT probe addresses | 928 B |
| `HM450.address` | 450K probe/address annotation | 7.64 MiB |
| `EPIC.address` | EPIC probe/address annotation | 13.15 MiB |
| `EPICv2.address` | EPICv2 probe/address annotation | 14.84 MiB |
| `MSA.address` | MSA probe/address annotation | 4.84 MiB |
| `KYCG.HM450.Mask.20220123` | Recommended 450K design mask | 1.79 MiB |
| `KYCG.EPIC.Mask.20220123` | Recommended EPIC design mask | 2.80 MiB |
| `KYCG.EPICv2.Mask.20230314` | Recommended EPICv2 design mask | 0.69 MiB |
| ExperimentHub database and index | Resolve cached resource identifiers offline | 10.00 MiB |

Exact bytes, source identifiers and SHA-256 hashes are recorded in
[`data/sesame-resources.json`](../data/sesame-resources.json). Clock coefficients,
normalization references and the blood imputation panel add about 6 MB and are
also bundled; they do not need rebuilding for normal use.

The cache is matched to **`sesame` 1.24.0 + `sesameData` 1.24.0**, tested with
**R 4.4.3**. By default the script extracts it into a temporary local cache and
requires those SeSAMe versions. If using another version, prepare matching
references and pass `--cache /absolute/cache` consistently to inspection,
dependency checks and processing. Do not assume a newer package can use the
bundled cache unchanged.

**MSA limitation:** the bundled 1.24.0 resources have an MSA address table but no
recommended MSA design mask. The script reports a warning,
`design_mask_available=FALSE` and `qc_status=review_design_mask_unavailable`.
Detection, bead-count and signal-related QC still run, but recommended
design-mask QC has not been completed. Preserve this limitation in the report.

## Identify the input and process one platform

Commands below run from the repository root; for an installed skill, use its
actual path. Keep raw IDATs, sample sheets and results outside Git.

```bash
Rscript epigenetic-clocks/scripts/preprocess_idat.R --check-deps

Rscript epigenetic-clocks/scripts/preprocess_idat.R \
  --input /path/to/idats --inspect

Rscript epigenetic-clocks/scripts/preprocess_idat.R \
  --input /path/to/idats --output-dir /path/to/new-idat-output
```

Accepted inputs are a recursive directory, IDAT prefix, one channel filename
ending in `_Grn.idat` or `_Red.idat` (optionally `.gz`), or a CSV sample sheet:

```csv
sample_id,idat_prefix
sample_01,idats/123456_R01C01
sample_02,idats/123456_R02C01
```

```bash
Rscript epigenetic-clocks/scripts/preprocess_idat.R \
  --sample-sheet /path/to/samples.csv --output-dir /path/to/new-idat-output
```

Both red and green files must exist; relative prefixes resolve beside the sample
sheet. Discovery checks channel filenames and pairs. Content inspection uses
SeSAMe's binary IDAT reader and actual probe-address signatures to identify
HM450, EPIC, EPICv2 or MSA. It rejects unsupported or unidentified arrays and
inconsistent paired address sets. It also reads IDAT barcode, chip type and
array position metadata and rejects mismatches between channels. Missing
identity metadata is recorded in `pair_identity_available`. These checks do not
independently prove a person's identity. Use one platform per run; split
mixed-platform inputs.

The identification stage uses the small signature resource, then processing
checks that platform's address annotation and available design mask. A beta-only
CSV/TSV bypasses this step: CpG count or filename alone cannot establish its
specific chip model. Record an unknown platform unless provenance supplies it.

## QC and beta values

The script applies **QCDPB in order**: recommended design masking (`Q`, where
available), Infinium-I channel inference (`C`), dye-bias correction (`D`), pOOBAH
detection masking (`P`), then noob background correction (`B`). Its default
detection cutoff is **0.05** (`p > 0.05` is masked), and the minimum bead count is
**1**. `--detection-p` and `--min-beads` explicitly change these settings. The
upstream functions and preprocessing codes are documented in SeSAMe's
[`R/open.R`](https://github.com/zwdzwd/sesame/blob/devel/R/open.R).

Nonfinite detection values are masked. Masked values remain `NA`; replicates
such as EPICv2 probe suffixes collapse by mean **after masking**, using SeSAMe
`getBetas`. Loci without a usable replicate remain missing. The beta matrix
keeps CpG and relevant non-CpG `ch` loci; samples align by locus name, not row
position. See SeSAMe's
[`R/sesame.R`](https://github.com/zwdzwd/sesame/blob/devel/R/sesame.R).

| Output | Contents |
|---|---|
| `input_inspection.csv` | Sample, inferred platform and available IDAT metadata |
| `betas.csv` | Masked, replicate-collapsed matrix for the clock script |
| `qc.csv` | Raw intensity/channel/dye statistics, detection failures, masks, missingness and design-mask availability |
| `samples.csv` | Sample/pair paths and input MD5 checksums |
| `run_info.txt`, `session_info.txt` | Processing settings and software versions |
| `probe_qc_*.csv.gz` | Optional per-probe mask, p-value and beta, with `--save-probe-qc` |

Review QC before imputation. The **20% missing-locus flag is a review heuristic**,
not a validated pass/fail threshold. Neither `review_metrics` nor producing
numbers certifies that a sample passed QC. High missingness, weak signal,
detection failure and an unavailable design mask need explicit interpretation.
Imputation cannot rescue a poor-quality sample. Upstream QC functions are in
[`R/QC.R`](https://github.com/zwdzwd/sesame/blob/devel/R/QC.R).

Existing output directories are refused. Input data are not modified or uploaded.
Recorded sample paths, identifiers and QC outputs remain private sample data.

## Imputation, then clocks

After QC review, pass `betas.csv` to `compute_clocks.py`. It imputes absent and
sample-specific masked clock CpGs using the bundled whole-blood methyLImp panel,
then computes the selected clocks. Inspect measured coverage, `n_imputed`,
`n_lowconf` and `n_unresolved`; imputed values are estimates, not measurements.
Unresolved models are reported as unavailable while other models continue.
Missing coefficient/reference files, including the methyLImp panel, are listed
together and stop the CLI before calculation. Restore them first. This differs
from a CpG absent from an available panel, which may use median fallback.

```bash
# Single person; age and sex are needed for GrimAge only.
python3 epigenetic-clocks/scripts/compute_clocks.py \
  --input /path/to/new-idat-output/betas.csv --age 45 --sex m --clocks all

# Example without GrimAge or personal metadata.
python3 epigenetic-clocks/scripts/compute_clocks.py \
  --input /path/to/new-idat-output/betas.csv --clocks horvath phenoage
```

`--age` and `--sex` apply to every sample column. For different people, preserve
their sample IDs and run each person's matrix with that person's metadata. The
IDAT sample sheet does not pass age/sex to the clock script. Neither IDAT
preprocessing nor non-GrimAge models require age or sex; age additionally enables
simple age acceleration for year-unit clocks.

## Install missing software or prepare a matching reference cache

When R, Python, packages or system libraries are missing, list them and request
permission to download and install the missing items. If this machine cannot
connect, the user can transfer standard installers/packages matching its
operating system, architecture and R/Python versions. Install those locally,
then repeat the dependency/resource checks. The 47.3 MiB annotation archive
remains part of the skill; no R or Python environment archive is created.

After approval, install R/SeSAMe dependencies following the official
[Bioconductor SeSAMe instructions](https://bioconductor.org/packages/sesame/).
The current Bioconductor release may install a different SeSAMe version; in that
case prepare its matching reference cache explicitly:

```bash
# Downloads public references; run only after download consent.
Rscript epigenetic-clocks/scripts/preprocess_idat.R \
  --init-cache --cache /path/to/sesame-cache

# Subsequent checks and processing remain local.
Rscript epigenetic-clocks/scripts/preprocess_idat.R \
  --check-deps --cache /path/to/sesame-cache

Rscript epigenetic-clocks/scripts/preprocess_idat.R \
  --input /path/to/idats --cache /path/to/sesame-cache \
  --output-dir /path/to/new-idat-output
```

For offline transfer of a custom cache, copy its complete directory including
ExperimentHub metadata and verify `--check-deps` on the target. Keep it paired
with the same package versions. Normal processing disables remote Hub access
and SeSAMe's alternate download fallback; `--init-cache` is the separate explicit
network operation. Rebuilding the blood panel also downloads public data and
requires separate consent if it was not already authorized.
