#!/usr/bin/env python3
"""Compute GrimAge (and other epigenetic clocks) from a DNA methylation
beta-value file. SELF-CONTAINED: depends only on pandas + numpy. The clock
coefficients and the imputation reference are vendored under ../data/, extracted
from the open-source biolearn library. No biolearn / torch / seaborn / network
needed at runtime.

The math is a faithful reimplementation of biolearn's GrimageModel and
LinearMethylationModel — verified to reproduce biolearn's outputs exactly.

GrimAge is a 2nd-generation, mortality-trained clock. It estimates DNAm
surrogates of 7 plasma proteins + smoking pack-years (V2 adds DNAm A1C & CRP),
then combines them with age and sex in a survival model. Because age and sex
feed the formula directly, --age and --sex are REQUIRED for GrimAge.

Input (auto-detected):
  - Long  : two columns, CpG id + beta (header names ignored; first col = CpG).
  - Matrix: first column = CpG id, remaining columns = one or more samples.

Usage:
  python compute_grimage.py --input betas.csv --age 45 --sex m
  python compute_grimage.py --input betas.csv --age 45 --sex m \
         --clocks grimage horvath hannum phenoage --sensitivity 40 42 47 49
"""
import argparse, os, sys
import pandas as pd, numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

# clock keyword -> (kind, coefficient file, transform)
#   kind "grim"   : GrimageModel algorithm (needs age+sex)
#   kind "linear" : intercept + sum(coef*beta), then transform
GRIM = {"grimagev1": "GrimAgeV1.csv", "grimagev2": "GrimAgeV2.csv"}
LINEAR = {
    "horvath":  ("Horvath1.csv", "horvath"),   # anti_trafo(sum + 0.696)
    "hannum":   ("Hannum.csv",   "identity"),
    "phenoage": ("PhenoAge.csv", "identity"),   # intercept row in file
}
GROUPS = {  # friendly group -> concrete clock keys
    "grimage": ["grimagev1", "grimagev2"],
    "grimagev1": ["grimagev1"], "grimagev2": ["grimagev2"],
    "horvath": ["horvath"], "hannum": ["hannum"], "phenoage": ["phenoage"],
}
NEEDS_AGE_SEX = {"grimagev1", "grimagev2"}


def anti_trafo(x, adult_age=20):
    return np.where(x < 0, (1 + adult_age) * np.exp(x) - 1, (1 + adult_age) * x + adult_age)


def load_betas(path):
    df = pd.read_csv(path)
    if df.shape[1] < 2:
        sys.exit(f"ERROR: expected >=2 columns in {path}, got {df.shape[1]}")
    df = df.rename(columns={df.columns[0]: "CpG"}).drop_duplicates("CpG").set_index("CpG")
    if df.shape[1] == 1:
        df.columns = ["Sample"]
    return df.apply(pd.to_numeric, errors="coerce")


def impute_missing(dnam, cpgs, ref):
    """Add any missing clock CpGs as rows filled from the population reference."""
    missing = list(dict.fromkeys(c for c in cpgs if c not in dnam.index))  # unique, order-preserving
    if not missing:
        return dnam, missing
    add = pd.DataFrame(index=missing, columns=dnam.columns, dtype=float)
    for s in dnam.columns:
        add[s] = ref.reindex(missing)
    return pd.concat([dnam, add]), missing


def predict_linear(dnam, coef_file, transform, age):
    """biolearn LinearMethylationModel: intercept row=1, sum(coef*beta), transform."""
    coef = pd.read_csv(os.path.join(DATA, coef_file), index_col=0)
    ccol = "CoefficientTraining" if "CoefficientTraining" in coef.columns else coef.columns[0]
    m = dnam.copy()
    m.loc["intercept"] = 1.0
    joined = coef.join(m, how="inner")              # rows = CpGs (+intercept) present in both
    betas = joined[ccol].to_numpy(dtype=float)
    mat = joined.iloc[:, 1:].to_numpy(dtype=float)  # (n_terms, n_samples)
    raw = (mat * betas[:, None]).sum(axis=0)
    vals = anti_trafo(raw + 0.696) if transform == "horvath" else raw
    return pd.Series(np.asarray(vals, dtype=float), index=joined.columns[1:])


def predict_grim(dnam, coef_file, age, sex_code):
    """biolearn GrimageModel. sex_code: 0->Female=1, else Female=0 (male)."""
    coef = pd.read_csv(os.path.join(DATA, coef_file))  # Y.pred, var, beta
    df = dnam.copy()
    df.loc["Age"] = age
    df.loc["Female"] = 1.0 if sex_code == 0 else 0.0
    df.loc["Intercept"] = 1.0

    sub_vals = {}            # sub-clock name -> Series over samples (numpy dot, no pandas align)
    cox = transform = None
    for name, grp in coef.groupby("Y.pred"):
        if name == "COX":
            cox = grp.set_index("var")["beta"]
        elif name == "transform":
            transform = grp.set_index("var")["beta"]
        else:
            cs = grp[grp["var"].isin(df.index)].set_index("var")["beta"]
            mat = df.reindex(cs.index).to_numpy(dtype=float)        # (n_cpg, n_samples)
            sub_vals[name] = pd.Series((mat * cs.to_numpy()[:, None]).sum(axis=0), index=df.columns)

    all_data = pd.DataFrame(sub_vals)
    all_data["Age"] = float(age)
    all_data["Female"] = 1.0 if sex_code == 0 else 0.0
    # COX = dot(sub-clock+Age+Female values, cox coefficients), numpy to avoid align
    cox_mat = all_data.reindex(columns=cox.index).to_numpy(dtype=float)  # (n_samples, n_terms)
    cox_score = pd.Series((cox_mat * cox.to_numpy()[None, :]).sum(axis=1), index=all_data.index)
    Y = (cox_score - transform["m_cox"]) / transform["sd_cox"]
    return Y * transform["sd_age"] + transform["m_age"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--age", type=float)
    ap.add_argument("--sex", choices=["m", "f", "male", "female"])
    ap.add_argument("--clocks", nargs="+", default=["grimage"],
                    help="any of: " + ", ".join(GROUPS))
    ap.add_argument("--sensitivity", type=float, nargs="*")
    args = ap.parse_args()

    # resolve requested clock groups
    keys = []
    for c in args.clocks:
        k = c.lower().replace("-", "").replace("_", "")
        if k not in GROUPS:
            sys.exit(f"ERROR: unknown clock '{c}'. Choose from: {', '.join(GROUPS)}")
        keys.extend(GROUPS[k])
    seen = set(); keys = [k for k in keys if not (k in seen or seen.add(k))]

    sex_code = None
    if any(k in NEEDS_AGE_SEX for k in keys):
        if args.age is None or args.sex is None:
            sys.exit("ERROR: GrimAge requires --age and --sex.")
        sex_code = 0 if args.sex in ("f", "female") else 1

    dnam = load_betas(args.input)
    ref = pd.read_csv(os.path.join(DATA, "sesame_450k_median.csv"), index_col=0).iloc[:, 0]
    samples = list(dnam.columns)
    print(f"Loaded {dnam.shape[0]} CpGs x {dnam.shape[1]} sample(s) from {args.input}")
    if sex_code is not None:
        print(f"Inputs: age={args.age}, sex={args.sex} (sex_code={sex_code})")
    print()

    rows = []
    for k in keys:
        if k in GRIM:
            coef = pd.read_csv(os.path.join(DATA, GRIM[k]))
            cpgs = list(dict.fromkeys(c for c in coef["var"] if str(c).startswith("cg")))
        else:
            cf = LINEAR[k][0]
            cpgs = [c for c in pd.read_csv(os.path.join(DATA, cf), index_col=0).index if str(c).startswith("cg")]
        d2, missing = impute_missing(dnam, cpgs, ref)
        cov = (len(cpgs) - len(missing)) / len(cpgs) * 100

        if k in GRIM:
            vals = predict_grim(d2, GRIM[k], args.age, sex_code)
        else:
            cf, tf = LINEAR[k]
            vals = predict_linear(d2, cf, tf, args.age)

        for s in samples:
            v = float(vals[s])
            accel = (v - args.age) if (args.age is not None) else np.nan
            rows.append(dict(sample=s, clock=k,
                             value=round(v, 2),
                             accel=(round(accel, 2) if not np.isnan(accel) else ""),
                             coverage=f"{cov:.1f}%", n_cpg=len(cpgs), n_imputed=len(missing)))

    res = pd.DataFrame(rows)
    print("=== Results ===")
    print(res.to_string(index=False))

    if args.sensitivity is not None and "grimagev2" in keys:
        cpgs = [c for c in pd.read_csv(os.path.join(DATA, "GrimAgeV2.csv"))["var"] if str(c).startswith("cg")]
        d2, _ = impute_missing(dnam, cpgs, ref)
        ages = sorted(set([args.age] + list(args.sensitivity)))
        print("\n=== GrimAgeV2 sensitivity to chronological age (sample 1) ===")
        for a in ages:
            v = float(predict_grim(d2, "GrimAgeV2.csv", a, sex_code)[samples[0]])
            print(f"  age={a:>4}: GrimAge={v:6.2f}, accel={v-a:+.2f}")

    print("\nJSON:", res.to_json(orient="records"))


if __name__ == "__main__":
    main()
