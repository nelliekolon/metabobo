#!/usr/bin/env python3
"""Step 5: does a model trained on one study work on another, and does harmonization help?

usage: python3 ml_crossstudy.py <crossstudy_dir> <train_study>=<table_dir> <test_study>=<table_dir> ...
Targets: sex (logistic regression, AUC) and age (ridge regression, MAE), the only labels every
study shares. Train study must be the reference used by match_studies.py (feature_map.tsv).

For each test study, only reference features matched in that study are used (feature_map.tsv).
Conditions:
  within-study CV  train/test inside the training study (5-fold, repeated), upper bound
  none             train on reference log2 values, predict test study's log2 values as they are
  per-study z      each feature standardized within each study before training / predicting
Test predictions are averaged per person (repeat visits are not independent), and a study's
plasma and serum are scored separately. 95% CIs: bootstrap over people (2000 resamples).
"""
import os, sys
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegressionCV, RidgeCV
from sklearn.model_selection import RepeatedStratifiedKFold, RepeatedKFold
from sklearn.metrics import roc_auc_score, mean_absolute_error

RNG = np.random.default_rng(0)
ALPHAS = np.logspace(-1, 4, 20)


def load(study, tdir):
    L = pd.read_csv(os.path.join(tdir, "feature_table_filtered_log2.tsv"), sep="\t", index_col=0)
    sh = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(tdir)), "samples.tsv"),
                     sep="\t", dtype=str, keep_default_na=False).set_index("sample")
    sh = sh.loc[[c for c in L.columns if sh.loc[c, "role"] == "subject"]]
    sh["age_num"] = pd.to_numeric(sh.age, errors="coerce")
    return L[sh.index].T, sh            # samples x features


def zscore(X):
    return (X - X.mean()) / (X.std() + 1e-9)


def per_person(pred, sh):
    df = pd.DataFrame({"pred": pred, "subject": sh.subject.values, "sex": sh.sex.values, "age": sh.age_num.values})
    return df.groupby("subject").agg(pred=("pred", "mean"), sex=("sex", "first"), age=("age", "first"))


def boot(fn, y, p, n=2000):
    idx = np.arange(len(y)); vals = []
    for _ in range(n):
        b = RNG.choice(idx, len(idx))
        try:
            vals.append(fn(y[b], p[b]))
        except ValueError:
            pass
    return np.percentile(vals, [2.5, 97.5]) if vals else (np.nan, np.nan)


def sex_model():
    return LogisticRegressionCV(Cs=np.logspace(-4, 1, 12), cv=5, scoring="roc_auc", max_iter=5000)


def score_sex(P):
    P = P[P.sex.isin(["male", "female"])]
    y = (P.sex == "male").astype(int).values; p = P.pred.values
    if y.min() == y.max():
        return f"n={len(y)} people, one sex only: not testable"
    auc = roc_auc_score(y, p); lo, hi = boot(roc_auc_score, y, p)
    return f"AUC {auc:.2f} [{lo:.2f}-{hi:.2f}]  (n={len(y)}: {y.sum()} M / {len(y)-y.sum()} F)"


def score_age(P, train_mean):
    P = P[P.age.notna()]
    y, p = P.age.values, P.pred.values
    mae = mean_absolute_error(y, p); base = mean_absolute_error(y, np.full_like(y, train_mean))
    lo, hi = boot(mean_absolute_error, y, p)
    r = np.corrcoef(y, p)[0, 1] if len(y) > 2 else np.nan
    return f"MAE {mae:.1f} y [{lo:.1f}-{hi:.1f}] vs {base:.1f} y predicting the training mean; r={r:.2f} (n={len(y)}, age {y.min():.0f}-{y.max():.0f})"


def within_cv(X, sh):
    out = {}
    s = sh.sex.isin(["male", "female"]).values
    Xs, ys = X[s], (sh.sex[s] == "male").astype(int).values
    pr = np.zeros(len(ys)); cnt = np.zeros(len(ys))
    for tr, te in RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0).split(Xs, ys):
        m = sex_model().fit(Xs.iloc[tr], ys[tr]); pr[te] += m.predict_proba(Xs.iloc[te])[:, 1]; cnt[te] += 1
    out["sex"] = score_sex(per_person(pr / cnt, sh[s]))
    a = sh.age_num.notna().values
    Xa, ya = X[a], sh.age_num[a].values
    pr = np.zeros(len(ya)); cnt = np.zeros(len(ya))
    for tr, te in RepeatedKFold(n_splits=5, n_repeats=5, random_state=0).split(Xa):
        m = RidgeCV(alphas=ALPHAS).fit(Xa.iloc[tr], ya[tr]); pr[te] += m.predict(Xa.iloc[te]); cnt[te] += 1
    out["age"] = score_age(per_person(pr / cnt, sh[a]), ya.mean())
    return out


def main():
    cdir = sys.argv[1]
    specs = [a.split("=", 1) for a in sys.argv[2:]]
    fmap = pd.read_csv(os.path.join(cdir, "feature_map.tsv"), sep="\t", index_col=0)
    train_name, train_dir = specs[0]
    Xr, shr = load(train_name, train_dir)
    rep = [f"train: {train_name} ({len(shr)} samples, {shr.subject.nunique()} people)"]
    cv_all = within_cv(zscore(Xr), shr)
    rep.append(f"\nwithin-study CV on {train_name} (all {Xr.shape[1]} features):")
    rep += [f"  {k:4s} {v}" for k, v in cv_all.items()]

    for name, tdir in specs[1:]:
        Xt_all, sht = load(name, tdir)
        pairs = fmap[name].dropna()
        ref_f = [f for f in pairs.index if f in Xr.columns]
        tst_f = [pairs[f] for f in ref_f if pairs[f] in Xt_all.columns]
        ref_f = [f for f in ref_f if pairs[f] in Xt_all.columns]
        A = Xr[ref_f]
        rep.append(f"\n=== test: {name}  ({len(ref_f)} shared features)")
        cv = within_cv(zscore(A), shr)
        rep.append(f"  within-study CV on {train_name}, shared features only (upper bound):")
        rep += [f"    {k:4s} {v}" for k, v in cv.items()]
        for matrix in sorted(sht.matrix.unique()):
            sel = (sht.matrix == matrix).values
            B = Xt_all.loc[sel, tst_f].copy(); B.columns = ref_f; shm = sht[sel]
            for cond in ["none", "per-study z"]:
                Atr = zscore(A) if cond == "per-study z" else A
                Bte = zscore(B) if cond == "per-study z" else B
                s = shr.sex.isin(["male", "female"]).values
                if cond == "none":   # standardize with TRAINING statistics only
                    mu, sd = A.mean(), A.std() + 1e-9
                    Atr, Bte = (A - mu) / sd, (B - mu) / sd
                m = sex_model().fit(Atr[s], (shr.sex[s] == "male").astype(int))
                sexr = score_sex(per_person(m.predict_proba(Bte)[:, 1], shm))
                a = shr.age_num.notna().values
                r = RidgeCV(alphas=ALPHAS).fit(Atr[a], shr.age_num[a])
                ager = score_age(per_person(r.predict(Bte), shm), shr.age_num[a].mean())
                rep.append(f"  {matrix:13s} {cond:12s} sex {sexr}")
                rep.append(f"  {'':13s} {'':12s} age {ager}")
    with open(os.path.join(cdir, "ml_report.txt"), "w") as fh:
        fh.write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
