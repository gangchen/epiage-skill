#!/usr/bin/env python3
"""Compute epigenetic aging clocks from a DNA methylation beta-value file.

SELF-CONTAINED: depends only on pandas + numpy. Clock coefficients and the
imputation reference are vendored under ../data/ (extracted from the open-source
biolearn library, trimmed to the CpGs the clocks use). No biolearn / torch /
seaborn / network needed at runtime. The math faithfully reimplements biolearn's
GrimageModel and LinearMethylationModel and reproduces biolearn's outputs.

24 clocks across several families (see --list-clocks). GrimAge requires --age and
--sex (it is age/sex adjusted). The other clocks don't, but passing --age lets the
tool report acceleration (= clock − chronological age) for the year-unit clocks.

Input (auto-detected): Long (CpG id + beta) or Matrix (CpG id + sample columns).

Usage:
  python compute_clocks.py --input betas.csv --age 45 --sex m
  python compute_clocks.py --input betas.csv --age 45 --sex m --clocks all
  python compute_clocks.py --list-clocks
"""
import argparse, os, sys
import pandas as pd, numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

# key -> dict(file, kind, tf=(type,offset) for linear, cat, unit, year)
#   kind "grim"   : GrimageModel (needs age+sex)
#   kind "linear" : transform(intercept + sum(coef*beta)); tf type "anti"|"lin"
CLOCKS = {
    # --- GrimAge (2nd-gen, mortality-trained; needs age+sex) ---
    "grimagev1": dict(file="GrimAgeV1.csv", kind="grim", cat="GrimAge (mortality)", unit="years", year=2019),
    "grimagev2": dict(file="GrimAgeV2.csv", kind="grim", cat="GrimAge (mortality)", unit="years", year=2022),
    # --- 1st-gen chronological-age clocks ---
    "horvath":    dict(file="Horvath1.csv",  kind="linear", tf=("anti", 0.696),            cat="1st-gen chronological", unit="years", year=2013),
    "horvath2":   dict(file="Horvath2.csv",  kind="linear", tf=("anti", -0.447119319),     cat="1st-gen (skin & blood)", unit="years", year=2018),
    "hannum":     dict(file="Hannum.csv",    kind="linear", tf=("lin", 0.0),               cat="1st-gen chronological", unit="years", year=2013),
    "lin":        dict(file="Lin.csv",       kind="linear", tf=("lin", 0.0),               cat="1st-gen chronological", unit="years", year=2016),
    "vidalbralo": dict(file="VidalBralo.csv",kind="linear", tf=("lin", 84.7),              cat="1st-gen chronological", unit="years", year=2018),
    "weidner":    dict(file="Weidner.csv",   kind="linear", tf=("lin", 38.0),              cat="1st-gen chronological", unit="years", year=2014),
    "garagnani":  dict(file="Garagnani.csv", kind="linear", tf=("lin", 0.0),               cat="1st-gen chronological", unit="years", year=2012),
    "bocklandt":  dict(file="Bocklandt.csv", kind="linear", tf=("lin", 0.0),               cat="1st-gen chronological", unit="years", year=2011),
    # --- tissue-/age-specific ---
    "pedbe":      dict(file="PEDBE.csv",            kind="linear", tf=("anti", -2.1),                 cat="pediatric (buccal)", unit="years", year=2019),
    "cortical":   dict(file="DNAmClockCortical.csv",kind="linear", tf=("anti", 0.577682570446177),   cat="brain cortex", unit="years", year=2020),
    # --- stochastic clocks ---
    "stoch":      dict(file="StocH.csv", kind="linear", tf=("lin", 59.8015666314217), cat="stochastic (Horvath)",  unit="years", year=2024),
    "stocp":      dict(file="StocP.csv", kind="linear", tf=("lin", 92.8310813279039), cat="stochastic (PhenoAge)", unit="years", year=2024),
    "stocz":      dict(file="StocZ.csv", kind="linear", tf=("lin", 64.8077188694894), cat="stochastic (mortality)", unit="years", year=2024),
    # --- 2nd-gen biological-age ---
    "phenoage":        dict(file="PhenoAge.csv",        kind="linear", tf=("lin", 0.0), cat="2nd-gen biological age", unit="years", year=2018),
    "hrsinchphenoage": dict(file="HRSInCHPhenoAge.csv", kind="linear", tf=("lin", 0.0), cat="2nd-gen biological age", unit="years", year=2022),
    # --- Ying 2022 causality-partitioned clocks ---
    "yingcausage":  dict(file="YingCausAge.csv",  kind="linear", tf=("lin", 0.0), cat="Ying causality",  unit="years", year=2022),
    "yingdamage":   dict(file="YingDamAge.csv",   kind="linear", tf=("lin", 0.0), cat="Ying damage",     unit="years", year=2022),
    "yingadaptage": dict(file="YingAdaptAge.csv", kind="linear", tf=("lin", 0.0), cat="Ying adaptation", unit="years", year=2022),
    # --- other aging-related (non-year units; no acceleration) ---
    "zhang":       dict(file="Zhang_10.csv",      kind="linear", tf=("lin", 0.0), cat="mortality risk",   unit="risk",       year=2019),
    "dunedinpace": dict(file="DunedinPACE.csv",   kind="dunedin",                 cat="pace of aging (3rd-gen)", unit="years/year", year=2022),
    "dunedinpoam": dict(file="DunedinPoAm38.csv", kind="linear", tf=("lin", 0.0), cat="pace of aging",    unit="years/year", year=2020),
    "dnamtl":      dict(file="DNAmTL.csv",        kind="linear", tf=("lin", 0.0), cat="telomere length",  unit="kb",         year=2019),
    "epitoc1":     dict(file="EpiTOC1.csv",       kind="linear", tf=("lin", 0.0), cat="mitotic (EpiTOC)", unit="score",      year=2016),
}
NEEDS_AGE_SEX = {k for k, v in CLOCKS.items() if v["kind"] == "grim"}
GROUPS = {
    "all": list(CLOCKS),
    "grimage": ["grimagev1", "grimagev2"],
    "core": ["grimagev1", "grimagev2", "horvath", "hannum", "phenoage"],
    "firstgen": ["horvath", "horvath2", "hannum", "lin", "vidalbralo", "weidner", "garagnani", "bocklandt"],
    "secondgen": ["grimagev1", "grimagev2", "phenoage", "hrsinchphenoage", "yingcausage", "yingdamage", "yingadaptage"],
    "thirdgen": ["dunedinpace", "dunedinpoam"],
}


def anti_trafo(x, adult_age=20):
    return np.where(x < 0, (1 + adult_age) * np.exp(x) - 1, (1 + adult_age) * x + adult_age)


# --- DunedinPACE quantile-normalization (numpy reimpl of biolearn, scipy-free) ---
def _rankdata_avg(a):
    """scipy.stats.rankdata(method='average'), numpy-only."""
    a = np.asarray(a, float)
    sorter = np.argsort(a, kind="quicksort")
    inv = np.empty(len(a), dtype=int); inv[sorter] = np.arange(len(a))
    s = a[sorter]
    obs = np.r_[True, s[1:] != s[:-1]]
    dense = obs.cumsum()[inv]
    count = np.r_[np.flatnonzero(obs), len(a)]
    return 0.5 * (count[dense] + count[dense - 1] + 1)


def _qnorm_to_target(data, target):
    """biolearn quantile_normalize_using_target (per-column, in place on a copy)."""
    st = np.sort(target); data = np.array(data, dtype=float)
    for col in data.T:
        r = _rankdata_avg(col); fl = np.floor(r).astype(int); hi = (r - fl) > 0.4
        col[hi] = 0.5 * (st[fl[hi] - 1] + st[fl[hi]])
        col[~hi] = st[fl[~hi] - 1]
    return data


def _hybrid_impute(dnam, src, required, threshold=0.8):
    """biolearn hybrid_impute: drop sparse rows, fill missing required from src (gold means)."""
    keep = dnam[dnam.notna().mean(axis=1) >= threshold]
    keep = keep.where(keep.notna(), keep.mean(axis=1), axis=0)
    miss = list(set(required) - set(keep.index))
    add = pd.DataFrame.from_dict({c: [src[c]] * dnam.shape[1] for c in miss},
                                 orient="index", columns=dnam.columns)
    return pd.concat([keep, add]).sort_index()


def predict_dunedin(dnam):
    """DunedinPACE: normalize sample to gold-standard distribution, then linear model.
    Returns (values Series, n_background, n_background_missing)."""
    gold = pd.read_csv(os.path.join(DATA, "DunedinPACE_Gold_Means.csv"), index_col=0)["mean"]
    coef = pd.read_csv(os.path.join(DATA, "DunedinPACE.csv"), index_col=0)["CoefficientTraining"]
    present = dnam.index.intersection(gold.index)
    n_bg = len(gold); n_miss = n_bg - len(present)
    filled = _hybrid_impute(dnam.loc[present], gold, list(gold.index))
    target = filled.index.map(gold.to_dict()).tolist()
    norm = pd.DataFrame(_qnorm_to_target(filled.values, target), index=filled.index, columns=filled.columns)
    mp = [c for c in coef.index if str(c).startswith("cg") and c in norm.index]
    vals = norm.loc[mp].multiply(coef.loc[mp], axis=0).sum(axis=0) + coef["intercept"]
    return vals, n_bg, n_miss


def model_cpgs(spec):
    """All CpG features a clock uses (cg + ch. control probes)."""
    if spec["kind"] == "grim":
        v = pd.read_csv(os.path.join(DATA, spec["file"]))["var"]
        feats = [c for c in v if str(c).startswith(("cg", "ch."))]
    else:
        idx = pd.read_csv(os.path.join(DATA, spec["file"]), index_col=0).index
        feats = [c for c in idx if str(c).startswith(("cg", "ch."))]
    return list(dict.fromkeys(feats))


def load_betas(path):
    df = pd.read_csv(path)
    if df.shape[1] < 2:
        sys.exit(f"ERROR: expected >=2 columns in {path}, got {df.shape[1]}")
    df = df.rename(columns={df.columns[0]: "CpG"}).drop_duplicates("CpG").set_index("CpG")
    if df.shape[1] == 1:
        df.columns = ["Sample"]
    return df.apply(pd.to_numeric, errors="coerce")


SD_THRESH = 0.08  # mLiftOver-style: a CpG whose population SD exceeds this is
                  # intrinsically variable, so imputing it from a median is unreliable.


def load_reference():
    """Imputation reference: whole-blood per-CpG median + SD (GSE40279, 656 blood
    samples), falling back to the tissue-agnostic sesame median for any CpG the
    blood panel lacks. Returns (median Series, sd Series). The blood median is a
    better guess for a blood sample than a global median; the SD drives the
    confidence flag on imputed values."""
    sesame = pd.read_csv(os.path.join(DATA, "sesame_450k_median.csv"), index_col=0).iloc[:, 0]
    bpath = os.path.join(DATA, "blood_reference_450k.csv")
    if not os.path.exists(bpath):
        return sesame, None  # graceful fallback to legacy global-median behaviour
    blood = pd.read_csv(bpath, index_col=0)
    ref_med = blood["median"].combine_first(sesame)  # blood preferred, sesame fallback
    return ref_med, blood["sd"]


def impute_missing(dnam, feats, ref_med, ref_sd=None):
    """Add missing clock features as rows filled from the (blood) reference median.
    Returns (filled_df, missing_list, n_lowconf) where n_lowconf counts imputed
    CpGs whose reference SD > SD_THRESH (or is unknown) — i.e. fills to distrust."""
    missing = list(dict.fromkeys(c for c in feats if c not in dnam.index))
    in_ref = [c for c in missing if c in ref_med.index]
    if ref_sd is not None:
        n_lowconf = sum(1 for c in missing
                        if (c not in ref_sd.index) or pd.isna(ref_sd.get(c)) or ref_sd[c] > SD_THRESH)
    else:
        n_lowconf = 0
    if not in_ref:
        return dnam, missing, n_lowconf
    add = pd.DataFrame({s: ref_med.reindex(in_ref) for s in dnam.columns}, index=in_ref)
    return pd.concat([dnam, add]), missing, n_lowconf


def predict_linear(dnam, spec):
    coef = pd.read_csv(os.path.join(DATA, spec["file"]), index_col=0)
    ccol = "CoefficientTraining" if "CoefficientTraining" in coef.columns else coef.columns[0]
    m = dnam.copy()
    m.loc["intercept"] = 1.0
    joined = coef.join(m, how="inner")
    betas = joined[ccol].to_numpy(dtype=float)
    mat = joined.iloc[:, 1:].to_numpy(dtype=float)
    raw = (mat * betas[:, None]).sum(axis=0)
    ttype, off = spec["tf"]
    vals = anti_trafo(raw + off) if ttype == "anti" else (raw + off)
    return pd.Series(np.asarray(vals, dtype=float), index=joined.columns[1:])


def predict_grim(dnam, spec, age, sex_code):
    coef = pd.read_csv(os.path.join(DATA, spec["file"]))  # Y.pred, var, beta
    df = dnam.copy()
    df.loc["Age"] = float(age)
    df.loc["Female"] = 1.0 if sex_code == 0 else 0.0
    df.loc["Intercept"] = 1.0
    sub_vals, cox, transform = {}, None, None
    for name, grp in coef.groupby("Y.pred"):
        if name == "COX":
            cox = grp.set_index("var")["beta"]
        elif name == "transform":
            transform = grp.set_index("var")["beta"]
        else:
            cs = grp[grp["var"].isin(df.index)].set_index("var")["beta"]
            mat = df.reindex(cs.index).to_numpy(dtype=float)
            sub_vals[name] = pd.Series((mat * cs.to_numpy()[:, None]).sum(axis=0), index=df.columns)
    all_data = pd.DataFrame(sub_vals)
    all_data["Age"] = float(age)
    all_data["Female"] = 1.0 if sex_code == 0 else 0.0
    cox_mat = all_data.reindex(columns=cox.index).to_numpy(dtype=float)
    cox_score = pd.Series((cox_mat * cox.to_numpy()[None, :]).sum(axis=1), index=all_data.index)
    Y = (cox_score - transform["m_cox"]) / transform["sd_cox"]
    return Y * transform["sd_age"] + transform["m_age"]


def list_clocks():
    print(f"{'key':16s} {'year':4s} {'unit':11s} category")
    print("-" * 64)
    for k, v in CLOCKS.items():
        star = " *needs age+sex" if v["kind"] == "grim" else ""
        print(f"{k:16s} {v['year']:<4d} {v['unit']:11s} {v['cat']}{star}")
    print("\nGroups for --clocks:", ", ".join(GROUPS))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input")
    ap.add_argument("--age", type=float)
    ap.add_argument("--sex", choices=["m", "f", "male", "female"])
    ap.add_argument("--clocks", nargs="+", default=["core"],
                    help="clock keys or groups (" + ", ".join(GROUPS) + "); default 'core'")
    ap.add_argument("--sensitivity", type=float, nargs="*")
    ap.add_argument("--list-clocks", action="store_true")
    args = ap.parse_args()

    if args.list_clocks:
        list_clocks(); return
    if not args.input:
        sys.exit("ERROR: --input is required (or use --list-clocks).")

    # resolve requested clocks
    keys = []
    for c in args.clocks:
        k = c.lower().replace("-", "").replace("_", "")
        if k in GROUPS:
            keys.extend(GROUPS[k])
        elif k in CLOCKS:
            keys.append(k)
        else:
            sys.exit(f"ERROR: unknown clock/group '{c}'. See --list-clocks.")
    seen = set(); keys = [k for k in keys if not (k in seen or seen.add(k))]

    sex_code = None
    if any(k in NEEDS_AGE_SEX for k in keys):
        if args.age is None or args.sex is None:
            sys.exit("ERROR: GrimAge requires --age and --sex (or drop grimage from --clocks).")
        sex_code = 0 if args.sex in ("f", "female") else 1

    dnam = load_betas(args.input)
    ref_med, ref_sd = load_reference()
    samples = list(dnam.columns)
    print(f"Loaded {dnam.shape[0]} CpGs x {dnam.shape[1]} sample(s) from {args.input}")
    if args.age is not None:
        print(f"Inputs: age={args.age}" + (f", sex={args.sex}" if args.sex else ""))
    print()

    rows = []
    for k in keys:
        spec = CLOCKS[k]
        if spec["kind"] == "dunedin":
            # self-normalizing; coverage measured against the ~20k background probes
            vals, n_feat, missing_n = predict_dunedin(dnam)
            cov = (n_feat - missing_n) / n_feat * 100
            feats = ["bg"] * n_feat; missing = ["m"] * missing_n
            n_lowconf = ""
        else:
            feats = model_cpgs(spec)
            d2, missing, n_lowconf = impute_missing(dnam, feats, ref_med, ref_sd)
            cov = (len(feats) - len(missing)) / len(feats) * 100 if feats else 100.0
            vals = predict_grim(d2, spec, args.age, sex_code) if spec["kind"] == "grim" else predict_linear(d2, spec)
        for s in samples:
            v = float(vals[s])
            accel = (v - args.age) if (args.age is not None and spec["unit"] == "years") else None
            rows.append(dict(sample=s, clock=k, category=spec["cat"], unit=spec["unit"],
                             value=round(v, 2), accel=("" if accel is None else round(accel, 2)),
                             coverage=f"{cov:.0f}%", n_feat=len(feats), n_imputed=len(missing),
                             n_lowconf=n_lowconf))

    res = pd.DataFrame(rows)
    print("=== Results ===")
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(res.to_string(index=False))

    if args.sensitivity is not None and "grimagev2" in keys:
        spec = CLOCKS["grimagev2"]
        d2, _, _ = impute_missing(dnam, model_cpgs(spec), ref_med, ref_sd)
        print("\n=== GrimAgeV2 sensitivity to chronological age (sample 1) ===")
        for a in sorted(set([args.age] + list(args.sensitivity))):
            v = float(predict_grim(d2, spec, a, sex_code)[samples[0]])
            print(f"  age={a:>4}: GrimAge={v:6.2f}, accel={v-a:+.2f}")

    print("\nJSON:", res.to_json(orient="records"))


if __name__ == "__main__":
    main()
