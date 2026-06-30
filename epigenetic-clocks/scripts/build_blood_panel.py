#!/usr/bin/env python3
"""Build the whole-blood methyLImp panel used to impute missing clock CpGs.

This skill is designed for human WHOLE BLOOD. Missing clock CpGs are imputed by
default with methyLImp (reduced-rank PCA regression) using the inter-CpG structure
of a blood reference. This script builds that reference once and writes
`data/blood_panel.npz` (the skill auto-detects and uses it).

Reference: GSE40279 (Hannum 2013), 656 whole-blood Illumina 450K samples — the
canonical whole-blood panel. ~1.24 GB download; we stream-extract only the ~26k
clock/background CpGs, so peak memory stays small.

The .npz stores, over the clock CpGs present in blood:
  cpgs   : CpG ids (order of mu/V columns / median / sd)
  mu     : per-CpG mean       (methyLImp baseline)
  V      : top-K PCA loadings (K x nCpG) — the blood correlation basis
  median : per-CpG median     (median fallback)
  sd     : per-CpG SD         (confidence flag)

Run on a normal connection (download is large; NCBI's GEO mirror can be slow):
    python build_blood_panel.py                 # auto-download + build
    python build_blood_panel.py LOCAL.txt.gz    # use an already-downloaded matrix

Requires pandas + numpy and `curl`.
"""
import os, sys, subprocess, gzip
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
GSE_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/"
           "GSE40279_series_matrix.txt.gz")
K = 50  # number of blood principal components retained


def clock_cpgs():
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
    print(f"Need {len(cpgs)} clock/background CpGs for the blood panel.")
    gz = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA, "GSE40279_series_matrix.txt.gz")
    if not os.path.exists(gz):
        print(f"Downloading GSE40279 (~1.24 GB) -> {gz} ...")
        subprocess.run(["curl", "-L", "-o", gz, GSE_URL], check=True)

    header, rows = None, {}
    with gzip.open(gz, "rt") as fh:
        for line in fh:
            if header is None and line.startswith('"ID_REF"'):
                header = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")]
                continue
            if line.startswith(('"cg', '"ch.', "cg", "ch.")):
                cpg = line.split("\t", 1)[0].strip().strip('"')
                if cpg in cpgs:
                    vals = [v.strip().strip('"') for v in line.rstrip("\n").split("\t")[1:]]
                    rows[cpg] = pd.to_numeric(pd.Series(vals), errors="coerce").to_numpy()
    if not rows:
        sys.exit("No CpG rows matched — check the matrix format.")

    mat = pd.DataFrame.from_dict(rows, orient="index").dropna()   # CpG x sample
    X = mat.to_numpy(dtype=np.float64).T                          # sample x CpG
    cpg_ids = list(mat.index)
    mu = X.mean(axis=0)
    Xc = X - mu
    # top-K PCA loadings (blood correlation basis) via SVD of centred matrix
    _, _, Vt = np.linalg.svd(Xc, full_matrices=False)
    V = Vt[:K]                                                    # K x nCpG
    out = os.path.join(DATA, "blood_panel.npz")
    np.savez_compressed(out, cpgs=np.array(cpg_ids), mu=mu.astype(np.float32),
                        V=V.astype(np.float32),
                        median=np.median(X, axis=0).astype(np.float32),
                        sd=X.std(axis=0).astype(np.float32))
    sd = X.std(axis=0)
    print(f"Wrote {out}: {len(cpg_ids)} CpGs, {X.shape[0]} blood samples, K={K} PCs.")
    print(f"  low-confidence CpGs (SD>0.08): {int((sd>0.08).sum())} ({(sd>0.08).mean()*100:.1f}%)")


if __name__ == "__main__":
    main()
