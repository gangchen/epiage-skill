#!/usr/bin/env python3
"""Compute epigenetic aging clocks from a DNA methylation beta-value file.

SELF-CONTAINED: depends only on pandas + numpy. Clock coefficients and the
imputation reference are vendored under ../data/ (extracted from the open-source
biolearn library, trimmed to the CpGs the clocks use). No biolearn / torch /
seaborn / network needed at runtime. The math faithfully reimplements biolearn's
GrimageModel and LinearMethylationModel and reproduces biolearn's outputs.

37 models across several families (see --list-clocks). GrimAge requires --age and
--sex (it is age/sex adjusted). The other clocks don't, but passing --age lets the
tool report acceleration (= clock − chronological age) for the year-unit clocks.

Input (auto-detected): Long (CpG id + beta) or Matrix (CpG id + sample columns).

Usage:
  python compute_clocks.py --input betas.csv --age 45 --sex m
  python compute_clocks.py --input betas.csv --age 45 --sex m --clocks all
  python compute_clocks.py --list-clocks
"""
import argparse, csv, gzip, os, sys
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
    # --- exposome / lifestyle methylation predictors (McCartney 2018, Reed) ---
    #     methylation "scores", not aging clocks; no acceleration. sigmoid outputs
    #     are relative scores in [0,1], identity outputs are raw predictor units.
    "smoking":     dict(file="Smoking.csv",       kind="linear", tf=("lin", 0.0),     cat="exposome: smoking",     unit="score", year=2018),
    "alcohol":     dict(file="Alcohol.csv",       kind="linear", tf=("lin", 0.0),     cat="exposome: alcohol",     unit="score", year=2018),
    "bmi":         dict(file="BMI_McCartney.csv", kind="linear", tf=("sigmoid", 0.0), cat="exposome: BMI",         unit="score", year=2018),
    "bmi_reed":    dict(file="BMI_Reed.csv",      kind="linear", tf=("lin", 0.0),     cat="exposome: BMI (Reed)",  unit="score", year=2020),
    "bodyfat":     dict(file="BodyFatMcCartney.csv",         kind="linear", tf=("sigmoid", 0.0), cat="exposome: body fat",    unit="score", year=2018),
    "hdl":         dict(file="HDLCholesterolMcCartney.csv",  kind="linear", tf=("sigmoid", 0.0), cat="exposome: HDL chol.",   unit="score", year=2018),
    "ldl":         dict(file="LDLCholesterolMcCartney.csv",  kind="linear", tf=("sigmoid", 0.0), cat="exposome: LDL chol.",   unit="score", year=2018),
    "totalchol":   dict(file="TotalCholesterolMcCartney.csv",kind="linear", tf=("sigmoid", 0.0), cat="exposome: total chol.", unit="score", year=2018),
    "education":   dict(file="EducationMcCartney.csv",       kind="linear", tf=("sigmoid", 0.0), cat="exposome: education",   unit="score", year=2018),
    # --- health / disease-risk methylation predictors ---
    "cvd":         dict(file="CVD_Westermann.csv",  kind="linear", tf=("sigmoid", 0.0),   cat="health: coronary heart disease", unit="risk", year=2020),
    "alzheimers":  dict(file="AD_Bahado-Singh.csv", kind="linear", tf=("sigmoid", 0.072), cat="health: Alzheimer's",            unit="risk", year=2022),
    "depression":  dict(file="DepressionBarbu.csv", kind="linear", tf=("lin", 0.0),       cat="health: depression",             unit="risk", year=2020),
}
NEEDS_AGE_SEX = {k for k, v in CLOCKS.items() if v["kind"] == "grim"}
EXPOSOME = ["smoking", "alcohol", "bmi", "bmi_reed", "bodyfat", "hdl", "ldl", "totalchol", "education"]
HEALTH = ["cvd", "alzheimers", "depression"]
AGING = [k for k in CLOCKS if k not in EXPOSOME + HEALTH]
GROUPS = {
    "all": list(CLOCKS),                 # every model (aging + exposome + health)
    "aging": AGING,                      # the 25 aging clocks only
    "grimage": ["grimagev1", "grimagev2"],
    "core": ["grimagev1", "grimagev2", "horvath", "hannum", "phenoage"],
    "firstgen": ["horvath", "horvath2", "hannum", "lin", "vidalbralo", "weidner", "garagnani", "bocklandt"],
    "secondgen": ["grimagev1", "grimagev2", "phenoage", "hrsinchphenoage", "yingcausage", "yingdamage", "yingadaptage"],
    "thirdgen": ["dunedinpace", "dunedinpoam"],
    "exposome": EXPOSOME,
    "health": HEALTH,
    "phenotypes": EXPOSOME + HEALTH,     # all non-aging predictors
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


def predict_dunedin(dnam):
    """DunedinPACE: normalize sample to gold-standard distribution, then linear model.
    Returns (values Series, n_background, per-sample missing counts).

    Fill each sample's absent or masked background probes from the gold means,
    as for a single-sample prediction. Never borrow another sample's betas or
    discard an observed beta because other samples have that probe masked.
    """
    gold = pd.read_csv(os.path.join(DATA, "DunedinPACE_Gold_Means.csv"), index_col=0)["mean"]
    coef = pd.read_csv(os.path.join(DATA, "DunedinPACE.csv"), index_col=0)["CoefficientTraining"]
    observed = dnam.reindex(gold.index).sort_index()
    n_bg = len(gold); n_miss = observed.isna().sum(axis=0)
    filled = observed.where(observed.notna(), gold, axis=0)
    target = gold.reindex(filled.index).to_numpy()
    norm = pd.DataFrame(_qnorm_to_target(filled.values, target), index=filled.index, columns=filled.columns)
    mp = [c for c in coef.index if str(c).startswith("cg")]
    vals = norm.reindex(mp).multiply(coef.loc[mp], axis=0).sum(axis=0, skipna=False) + coef["intercept"]
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
    """Read CpG-row CSV/TSV (optionally gzip); allow explicit missing betas.

    SeSAMe masks are NA values, not zero. Invalid text and out-of-range or
    infinite values are rejected rather than silently treated as masks.
    Duplicate CpG rows must be resolved by the array preprocessing workflow.
    """
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        first_line = handle.readline()
    sep = "\t" if "\t" in first_line else ","
    headers = next(csv.reader([first_line], delimiter=sep), [])
    if len(headers) < 2:
        raise ValueError(f"expected >=2 comma- or tab-delimited columns in {path}")
    sample_headers = headers[1:]
    if any(not name.strip() for name in sample_headers):
        raise ValueError("sample column names must not be empty")
    if len(set(sample_headers)) != len(sample_headers):
        raise ValueError("duplicate sample column names; give every sample a unique name")
    df = pd.read_csv(path, sep=sep, dtype={headers[0]: str})
    if df.shape[1] < 2:
        raise ValueError(f"expected >=2 columns in {path}, got {df.shape[1]}")
    cpgs = df.iloc[:, 0].str.strip()
    if cpgs.isna().any() or cpgs.eq("").any():
        raise ValueError("CpG identifiers must not be empty")
    if cpgs.str.match(r"^cg[0-9]+_").any():
        raise ValueError("array probe suffixes detected; export SeSAMe betas with collapseToPfx=TRUE "
                         "after masking before calculating clocks")
    if not cpgs.str.match(r"^(cg[0-9]+$|ch\.)").any():
        raise ValueError("no recognized Illumina methylation probe IDs (cg.../ch...) in the first column; "
                         "this is not a supported array beta matrix")
    if cpgs.duplicated().any():
        duplicate = cpgs[cpgs.duplicated()].iloc[0]
        raise ValueError(f"duplicate CpG identifier {duplicate!r}; collapse array probe replicates "
                         "during preprocessing (e.g. SeSAMe EPICv2 probe collapse), then retry")
    df = df.iloc[:, 1:].copy()
    df.index = pd.Index(cpgs, name="CpG")
    if df.shape[1] == 1 and df.columns[0] == "Beta_value":
        df.columns = ["Sample"]
    try:
        numeric = df.apply(pd.to_numeric, errors="raise")
    except (ValueError, TypeError) as exc:
        raise ValueError("beta values must be numeric or explicit missing values (NA/blank)") from exc
    values = numeric.to_numpy(dtype=float)
    invalid = (~np.isnan(values)) & ((~np.isfinite(values)) | (values < 0) | (values > 1))
    if invalid.any():
        row, col = np.argwhere(invalid)[0]
        raise ValueError(f"beta for {numeric.index[row]!r}, sample {numeric.columns[col]!r} "
                         "must be finite and within [0, 1], or NA when masked")
    empty_samples = numeric.columns[numeric.notna().sum(axis=0) == 0].tolist()
    if empty_samples:
        raise ValueError(f"no observed beta values for sample(s): {', '.join(empty_samples)}")
    return numeric


SD_THRESH = 0.08  # a CpG whose whole-blood SD exceeds this is intrinsically
                  # variable, so imputing it from a reference is unreliable.


class ImputeRef:
    """Imputation backend for missing CpGs. THIS SKILL IS DESIGNED FOR HUMAN WHOLE
    BLOOD. Missing clock CpGs (e.g. EPIC-trained CpGs absent on an MSA export) are
    imputed by default with **methyLImp** — reduced-rank (PCA) regression that
    predicts a missing CpG from the sample's observed CpGs using the inter-CpG
    correlation structure of a whole-blood reference panel. CpGs the panel lacks
    fall back to the blood median. With no panel installed it degrades gracefully
    to blood-median, then global sesame median.

    Install the blood panel once with scripts/build_blood_panel.py (writes
    data/blood_panel.npz from GSE40279, 656 whole-blood 450K samples)."""

    def __init__(self):
        self.median = pd.read_csv(os.path.join(DATA, "sesame_450k_median.csv"), index_col=0).iloc[:, 0]
        self.sd = None
        self.mu = self.V = None
        self.cidx = {}
        self.mode = "global-median"
        npz = os.path.join(DATA, "blood_panel.npz")
        csv = os.path.join(DATA, "blood_reference_450k.csv")
        if os.path.exists(npz):
            d = np.load(npz, allow_pickle=True)
            cpgs = [str(c) for c in d["cpgs"]]
            self.cidx = {c: i for i, c in enumerate(cpgs)}
            self.mu, self.V = d["mu"], d["V"]
            self.median = pd.Series(d["median"], index=cpgs).combine_first(self.median)
            self.sd = pd.Series(d["sd"], index=cpgs)
            self.mode = "methylimp" if getattr(self.V, "size", 0) else "blood-median"
        elif os.path.exists(csv):
            blood = pd.read_csv(csv, index_col=0)
            self.median = blood["median"].combine_first(self.median)
            self.sd = blood["sd"]
            self.mode = "blood-median"

    def _methylimp(self, sample, targets):
        """methyLImp: predict panel CpGs in `targets` for one sample (a Series of
        its observed betas) via reduced-rank projection onto the blood PCA basis."""
        obs = [c for c in sample.index if c in self.cidx]
        tgt = [c for c in targets if c in self.cidx]
        if not tgt or len(obs) <= self.V.shape[0]:
            return {}
        oidx = [self.cidx[c] for c in obs]
        o = sample.loc[obs].to_numpy(float) - self.mu[oidx]
        z, *_ = np.linalg.lstsq(self.V[:, oidx].T, o, rcond=None)   # blood PC scores
        tidx = [self.cidx[c] for c in tgt]
        pred = self.mu[tidx] + z @ self.V[:, tidx]
        return dict(zip(tgt, np.clip(pred, 0.0, 1.0)))


def impute_missing(dnam, feats, ref):
    """Fill absent rows and sample-specific NA masks without changing observations.

    Return (filled_df, missing_by_sample, n_lowconf_by_sample). Unresolvable
    features remain NA so the caller can flag an unavailable sample/model;
    silently dropping their coefficients would produce an incorrect score.
    """
    feats = list(dict.fromkeys(feats))
    result = dnam.reindex(dnam.index.union(feats, sort=False)).copy()
    missing, n_lowconf = {}, {}
    for s in dnam.columns:
        missing[s] = [c for c in feats if pd.isna(result.at[c, s])]
        filled = ref._methylimp(dnam[s].dropna(), missing[s]) if ref.mode == "methylimp" else {}
        for c in missing[s]:
            result.at[c, s] = filled.get(c, ref.median.get(c, np.nan))
        n_lowconf[s] = sum(
            1 for c in missing[s] if pd.notna(result.at[c, s]) and
            (ref.sd is None or pd.isna(ref.sd.get(c)) or ref.sd[c] > SD_THRESH)
        )
    return result, missing, n_lowconf


def predict_linear(dnam, spec):
    coef = pd.read_csv(os.path.join(DATA, spec["file"]), index_col=0)
    ccol = "CoefficientTraining" if "CoefficientTraining" in coef.columns else coef.columns[0]
    m = dnam.copy()
    m.loc["intercept"] = 1.0
    betas = coef[ccol].to_numpy(dtype=float)
    mat = m.reindex(coef.index).to_numpy(dtype=float)
    raw = (mat * betas[:, None]).sum(axis=0)
    ttype, off = spec["tf"]
    if ttype == "anti":
        vals = anti_trafo(raw + off)
    elif ttype == "sigmoid":
        vals = 1.0 / (1.0 + np.exp(-(raw + off)))
    else:  # "lin"
        vals = raw + off
    return pd.Series(np.asarray(vals, dtype=float), index=dnam.columns)


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
            cs = grp.set_index("var")["beta"]
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


def check_clock_resources(keys):
    required = {CLOCKS[k]["file"] for k in keys}
    required.add("sesame_450k_median.csv")
    if any(CLOCKS[k]["kind"] != "dunedin" for k in keys):
        required.add("blood_panel.npz")
    if "dunedinpace" in keys:
        required.add("DunedinPACE_Gold_Means.csv")
    missing = [name for name in sorted(required)
               if not os.path.isfile(os.path.join(DATA, name)) or os.path.getsize(os.path.join(DATA, name)) == 0]
    if missing:
        raise ValueError("missing local clock/imputation resources:\n  - " + "\n  - ".join(missing) +
                         "\nAsk the user whether to download the missing resources, or import a complete "
                         "offline skill bundle. No download was attempted.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input")
    ap.add_argument("--age", type=float)
    ap.add_argument("--sex", choices=["m", "f", "male", "female"])
    ap.add_argument("--clocks", nargs="+", default=["core"],
                    help="clock keys or groups (" + ", ".join(GROUPS) + "); default 'core'")
    ap.add_argument("--sensitivity", type=float, nargs="*")
    ap.add_argument("--list-clocks", action="store_true")
    ap.add_argument("--check-resources", action="store_true",
                    help="check local coefficients and imputation references without processing a sample")
    args = ap.parse_args()

    if args.list_clocks:
        list_clocks(); return
    if not args.input and not args.check_resources:
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

    try:
        check_clock_resources(keys)
    except ValueError as exc:
        sys.exit(f"ERROR: {exc}")
    if args.check_resources:
        print(f"Local coefficient and imputation resources ready for {len(keys)} model(s).")
        return

    sex_code = None
    if any(k in NEEDS_AGE_SEX for k in keys):
        if args.age is None or args.sex is None:
            sys.exit("ERROR: GrimAge requires --age and --sex (or drop grimage from --clocks).")
        sex_code = 0 if args.sex in ("f", "female") else 1

    try:
        dnam = load_betas(args.input)
    except (ValueError, OSError) as exc:
        sys.exit(f"ERROR: {exc}")
    if args.age is not None and (not np.isfinite(args.age) or args.age < 0):
        sys.exit("ERROR: --age must be finite and nonnegative.")
    if args.sensitivity and any(not np.isfinite(a) or a < 0 for a in args.sensitivity):
        sys.exit("ERROR: --sensitivity ages must be finite and nonnegative.")
    ref = ImputeRef()
    samples = list(dnam.columns)
    print(f"Loaded {dnam.shape[0]} CpGs x {dnam.shape[1]} sample(s) from {args.input}")
    if args.age is not None:
        print(f"Inputs: age={args.age}" + (f", sex={args.sex}" if args.sex else ""))
    _mode = {"methylimp": "methyLImp (blood panel)", "blood-median": "blood median",
             "global-median": "global median (no blood panel installed — see build_blood_panel.py)"}[ref.mode]
    print(f"Tissue: whole blood | imputation of missing CpGs: {_mode}")
    print()

    rows = []
    for k in keys:
        spec = CLOCKS[k]
        if spec["kind"] == "dunedin":
            # self-normalizing; coverage measured against the ~20k background probes
            vals, n_feat, missing_counts = predict_dunedin(dnam)
            unresolved = {s: 0 for s in samples}
            n_lowconf = {s: "" for s in samples}
        else:
            feats = model_cpgs(spec)
            d2, missing, n_lowconf = impute_missing(dnam, feats, ref)
            n_feat = len(feats)
            missing_counts = {s: len(missing[s]) for s in samples}
            unresolved = d2.reindex(feats).isna().sum(axis=0)
            vals = predict_grim(d2, spec, args.age, sex_code) if spec["kind"] == "grim" else predict_linear(d2, spec)
        for s in samples:
            v = float(vals[s])
            missing_n = int(missing_counts[s])
            unresolved_n = int(unresolved[s])
            cov = (n_feat - missing_n) / n_feat * 100 if n_feat else 100.0
            available = np.isfinite(v) and not unresolved_n
            if not available:
                v = np.nan
                print(f"WARNING: {s}/{k} unavailable: {unresolved_n} required CpGs "
                      "have no usable value/reference, or prediction is nonfinite.", file=sys.stderr)
            accel = (v - args.age) if (available and args.age is not None and spec["unit"] == "years") else None
            rows.append(dict(sample=s, clock=k, category=spec["cat"], unit=spec["unit"],
                             value=round(v, 2), accel=(np.nan if not available else
                                 "" if accel is None else round(accel, 2)),
                             coverage=f"{cov:.0f}%", n_feat=n_feat, n_missing=missing_n,
                             n_imputed=missing_n - unresolved_n, n_unresolved=unresolved_n,
                             n_lowconf=n_lowconf[s], status="ok" if available else "unavailable"))

    res = pd.DataFrame(rows)
    print("=== Results ===")
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(res.to_string(index=False))

    if args.sensitivity is not None and "grimagev2" in keys:
        spec = CLOCKS["grimagev2"]
        d2, _, _ = impute_missing(dnam, model_cpgs(spec), ref)
        print("\n=== GrimAgeV2 sensitivity to chronological age (sample 1) ===")
        for a in sorted(set([args.age] + list(args.sensitivity))):
            v = float(predict_grim(d2, spec, a, sex_code)[samples[0]])
            print(f"  age={a:>4}: GrimAge={v:6.2f}, accel={v-a:+.2f}")

    print("\nJSON:", res.to_json(orient="records"))


if __name__ == "__main__":
    main()
