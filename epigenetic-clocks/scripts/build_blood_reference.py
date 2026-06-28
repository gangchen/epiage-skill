#!/usr/bin/env python3
"""Build the whole-blood imputation reference (per-CpG median + SD) used for
mLiftOver-style imputation + confidence flagging.

Why: imputing a missing clock CpG from a *blood* median is a better guess than a
tissue-agnostic global median for a blood sample, and the per-CpG SD lets us flag
fills that are intrinsically unreliable (mLiftOver uses SD>0.08). The skill works
without this file (it falls back to the bundled sesame global median), but drops
in `data/blood_reference_450k.csv` automatically once you generate it here.

Reference dataset: GSE40279 (Hannum 2013), 656 whole-blood Illumina 450K samples —
the canonical whole-blood reference. ~1.24 GB download; we stream-extract only the
~26k clock/background CpGs, so peak memory stays tiny.

Run on a normal connection (the download is large):
    python build_blood_reference.py            # auto-download + build
    python build_blood_reference.py LOCAL.txt.gz   # use an already-downloaded matrix

Requires: pandas, numpy, and `curl` + `zcat`/`gzip` on PATH.
"""
import os, sys, subprocess, gzip
import pandas as pd, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
GSE_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/"
           "GSE40279_series_matrix.txt.gz")


def clock_cpgs():
    """Union of every CpG/feature used by any vendored clock + DunedinPACE background."""
    feats = set()
    for f in os.listdir(DATA):
        p = os.path.join(DATA, f)
        if f in ("GrimAgeV1.csv", "GrimAgeV2.csv"):
            feats |= {c for c in pd.read_csv(p)["var"] if str(c).startswith(("cg", "ch."))}
        elif f.endswith(".csv") and f not in ("sesame_450k_median.csv", "blood_reference_450k.csv"):
            idx = pd.read_csv(p, index_col=0).index
            feats |= {c for c in idx if str(c).startswith(("cg", "ch."))}
    return feats


def main():
    cpgs = clock_cpgs()
    print(f"Need {len(cpgs)} clock/background CpGs.")
    gz = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA, "GSE40279_series_matrix.txt.gz")
    if not os.path.exists(gz):
        print(f"Downloading GSE40279 (~1.24 GB) -> {gz} ...")
        subprocess.run(["curl", "-L", "-o", gz, GSE_URL], check=True)

    # stream rows, keep only our CpGs + the sample header line
    header, rows = None, {}
    with gzip.open(gz, "rt") as fh:
        for line in fh:
            if header is None and line.startswith('"ID_REF"'):
                header = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")]
                continue
            if line[:4] in ('"cg0', '"cg1', '"cg2', '"ch.') or line.startswith(("cg", "ch.")):
                cpg = line.split("\t", 1)[0].strip().strip('"')
                if cpg in cpgs:
                    vals = [v.strip().strip('"') for v in line.rstrip("\n").split("\t")[1:]]
                    rows[cpg] = pd.to_numeric(pd.Series(vals), errors="coerce").values
    if not rows:
        sys.exit("No CpG rows matched — check the matrix format.")
    mat = pd.DataFrame.from_dict(rows, orient="index")
    out = pd.DataFrame({"median": mat.median(axis=1, skipna=True),
                        "sd": mat.std(axis=1, skipna=True)}).dropna(subset=["median"])
    out.index.name = "CpG"
    dest = os.path.join(DATA, "blood_reference_450k.csv")
    out.to_csv(dest)
    n_samples = (len(header) - 1) if header else mat.shape[1]
    print(f"Wrote {dest}: {len(out)} CpGs from {n_samples} whole-blood samples.")
    print(f"  low-confidence CpGs (SD>0.08): {int((out['sd']>0.08).sum())} "
          f"({(out['sd']>0.08).mean()*100:.1f}%)")


if __name__ == "__main__":
    main()
