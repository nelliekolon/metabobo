#!/usr/bin/env python3
"""Step 4: match features across studies (m/z + mapped retention time).

usage: python3 match_studies.py <out_dir> <ref_study>=<table_dir> <study>=<table_dir> ...
Each table_dir must hold feature_info.tsv and feature_table_filtered_log2.tsv (check_table.py).
Only features that passed each study's 80% filter are matched.

Per study vs the reference:
 1. Anchors: the 1500 most intense features on each side; m/z matches within 5 ppm that are
    mutually unique (exactly one candidate in each direction).
 2. RT map: isotonic (order-preserving) regression of ref RT on study RT through the anchors,
    refit after dropping anchors with residual > 4 MAD. Chromatography differs between studies,
    so only elution ORDER is assumed to be shared, not the time scale.
 3. Matching: all features, m/z within 10 ppm and |mapped RT - ref RT| <= RT tolerance
    (3 x robust SD of anchor residuals, floor 6 s); one-to-one, best combined score first.
Validation: the known compounds below are NOT used for fitting. For each, the mapped RT of the
study's hit is compared with the reference's hit.

Outputs in out_dir: rt_map_<study>.tsv (anchors + fit), matches_<study>.tsv,
feature_map.tsv (one row per reference feature, matched feature ID per study), match_report.txt
"""
import os, sys
import numpy as np, pandas as pd
from sklearn.isotonic import IsotonicRegression

KNOWN = {"creatinine": 114.0662, "proline": 116.0706, "betaine": 118.0863, "leucine/isoleucine": 132.1019,
         "carnitine": 162.1125, "phenylalanine": 166.0863, "hippuric acid": 180.0655, "tyrosine": 182.0812,
         "caffeine": 195.0877, "tryptophan": 205.0972, "cortisol": 363.2166, "LPC 16:0": 496.3398,
         "LPC 18:2": 520.3398}
N_ANCHOR, ANCHOR_PPM, MATCH_PPM, MIN_RT_TOL = 1500, 5.0, 10.0, 6.0


def load(table_dir):
    info = pd.read_csv(os.path.join(table_dir, "feature_info.tsv"), sep="\t", index_col=0)
    L = pd.read_csv(os.path.join(table_dir, "feature_table_filtered_log2.tsv"), sep="\t", index_col=0)
    f = info.loc[L.index, ["mz", "rt_s"]].copy()
    f["inten"] = L.median(axis=1)
    return f


def ppm_pairs(a, b, tol):
    """all (i, j) with |mz_a - mz_b| / mz_b <= tol ppm; a, b are DataFrames."""
    bs = b.sort_values("mz"); bm = bs.mz.values
    out = []
    for i, m in zip(a.index, a.mz.values):
        lo, hi = np.searchsorted(bm, [m * (1 - tol * 1e-6), m * (1 + tol * 1e-6)])
        out += [(i, j) for j in bs.index[lo:hi]]
    return out


def fit_rt_map(q, ref):
    qa = q.nlargest(N_ANCHOR, "inten"); ra = ref.nlargest(N_ANCHOR, "inten")
    pairs = ppm_pairs(qa, ra, ANCHOR_PPM)
    cq = pd.Series([i for i, _ in pairs]).value_counts(); cr = pd.Series([j for _, j in pairs]).value_counts()
    uniq = [(i, j) for i, j in pairs if cq[i] == 1 and cr[j] == 1]
    A = pd.DataFrame({"q": [i for i, _ in uniq], "ref": [j for _, j in uniq]})
    A["rt_q"] = q.loc[A.q, "rt_s"].values; A["rt_ref"] = ref.loc[A.ref, "rt_s"].values
    keep = np.ones(len(A), bool)
    for _ in range(3):
        iso = IsotonicRegression(out_of_bounds="clip").fit(A.rt_q[keep], A.rt_ref[keep])
        res = A.rt_ref - iso.predict(A.rt_q)
        mad = np.median(np.abs(res[keep] - np.median(res[keep])))
        keep = np.abs(res) <= 4 * 1.4826 * mad + 1e-9
    A["used"] = keep; A["mapped"] = iso.predict(A.rt_q); A["resid"] = A.rt_ref - A.mapped
    # In-sample residuals are ~0 (isotonic fits its anchors exactly), so the tolerance comes from
    # 5-fold cross-validated residuals: each anchor predicted by a map fitted without it.
    used = A.index[keep]; rng = np.random.default_rng(0); fold = rng.integers(0, 5, len(used))
    cv = pd.Series(np.nan, index=A.index)
    for k in range(5):
        tr, te = used[fold != k], used[fold == k]
        m = IsotonicRegression(out_of_bounds="clip").fit(A.rt_q[tr], A.rt_ref[tr])
        cv[te] = A.rt_ref[te] - m.predict(A.rt_q[te])
    A["cv_resid"] = cv
    rt_tol = max(MIN_RT_TOL, 3 * 1.4826 * np.median(np.abs(cv[used] - np.median(cv[used]))))
    return iso, A, rt_tol


def match(q, ref, iso, rt_tol):
    q = q.assign(rt_mapped=iso.predict(q.rt_s))
    cand = []
    for i, j in ppm_pairs(q, ref, MATCH_PPM):
        dppm = abs(q.mz[i] - ref.mz[j]) / ref.mz[j] * 1e6
        drt = abs(q.rt_mapped[i] - ref.rt_s[j])
        if drt <= rt_tol:
            cand.append((dppm / MATCH_PPM + drt / rt_tol, i, j, dppm, drt))
    cand.sort()
    used_q, used_r, out = set(), set(), []
    for s, i, j, dppm, drt in cand:
        if i not in used_q and j not in used_r:
            used_q.add(i); used_r.add(j)
            out.append(dict(feature=i, ref_feature=j, ppm=round(dppm, 2), rt_diff=round(drt, 1), score=round(s, 3)))
    return pd.DataFrame(out), q


def best_hit(f, mz):
    h = f[np.abs(f.mz - mz) / mz * 1e6 < 5]
    return h.inten.idxmax() if len(h) else None


def main():
    out_dir = sys.argv[1]
    specs = [a.split("=", 1) for a in sys.argv[2:]]
    ref_name, ref_dir = specs[0]
    ref = load(ref_dir)
    os.makedirs(out_dir, exist_ok=True)
    fmap = pd.DataFrame(index=ref.index); fmap["mz"] = ref.mz; fmap["rt_s"] = ref.rt_s
    rep = [f"reference: {ref_name} ({len(ref)} filtered features)"]
    for name, tdir in specs[1:]:
        q = load(tdir)
        iso, A, rt_tol = fit_rt_map(q, ref)
        A.to_csv(os.path.join(out_dir, f"rt_map_{name}.tsv"), sep="\t", index=False)
        M, qm = match(q, ref, iso, rt_tol)
        M.to_csv(os.path.join(out_dir, f"matches_{name}.tsv"), sep="\t", index=False)
        fmap[name] = M.set_index("ref_feature").feature.reindex(fmap.index)
        cvr = A.cv_resid[A.used]
        rep.append(f"\n{name} ({len(q)} filtered features)")
        rep.append(f"  anchors: {A.used.sum()} used of {len(A)} unique m/z pairs; cross-validated RT error: "
                   f"median |err| {cvr.abs().median():.1f} s, 90th pct {cvr.abs().quantile(0.9):.1f} s "
                   f"-> RT tolerance {rt_tol:.1f} s")
        rep.append(f"  RT range {q.rt_s.min():.0f}-{q.rt_s.max():.0f} s mapped onto ref "
                   f"{A.mapped.min():.0f}-{A.mapped.max():.0f} s")
        rep.append(f"  matched: {len(M)} features ({len(M)/len(q):.0%} of this study, {len(M)/len(ref):.0%} of reference)")
        rep.append("  known compounds (held out; not used to fit the RT map):")
        # strongest hit = what check_table.py reports; it can be an isomer, so also list every
        # candidate at that m/z (mapped RT) to tell a mapping error from a wrong-peak pick.
        for cname, mz in KNOWN.items():
            hq, hr = best_hit(qm, mz), best_hit(ref, mz)
            if hq is None or hr is None:
                rep.append(f"    {cname:20s} not in both filtered tables"); continue
            got = M[(M.feature == hq)]
            same = len(got) and got.ref_feature.iloc[0] == hr
            cands = qm[np.abs(qm.mz - mz) / mz * 1e6 < 5].sort_values("rt_mapped")
            near = (np.abs(cands.rt_mapped - ref.rt_s[hr]) <= rt_tol).any()
            verdict = ("matched" if same else
                       "wrong peak picked (another candidate maps within tolerance)" if near else "MAPPING MISS")
            rep.append(f"    {cname:20s} ref {ref.rt_s[hr]:4.0f} s | strongest hit {qm.rt_s[hq]:4.0f} -> {qm.rt_mapped[hq]:4.0f} s "
                       f"| {verdict} | candidates mapped: {', '.join(f'{x:.0f}' for x in cands.rt_mapped)}")
    others = [n for n, _ in specs[1:]]
    n_in = fmap[others].notna().sum(axis=1)
    rep.append(f"\nreference features matched in k other studies: " +
               ", ".join(f"k={k}: {(n_in == k).sum()}" for k in range(len(others) + 1)))
    fmap.to_csv(os.path.join(out_dir, "feature_map.tsv"), sep="\t")
    with open(os.path.join(out_dir, "match_report.txt"), "w") as fh:
        fh.write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
