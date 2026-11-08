#!/usr/bin/env python3
"""Step 3: sanity checks on a feature table.
usage: python3 check_table.py <table_dir> [min_fraction]"""
import sys, os, numpy as np, pandas as pd
d = sys.argv[1]; frac = float(sys.argv[2]) if len(sys.argv) > 2 else 0.8
tab = pd.read_csv(os.path.join(d, "feature_table.tsv"), sep="\t", index_col=0)
info = pd.read_csv(os.path.join(d, "feature_info.tsv"), sep="\t", index_col=0)
keep = tab.notna().mean(axis=1) >= frac
X = tab[keep]
print(f"features kept (detected in >={frac:.0%} of samples): {keep.sum()} of {len(tab)}")
L = np.log2(X)
L = L.fillna(L.min(axis=1) - 1, axis=0) if False else L.apply(lambda r: r.fillna(r.min() - 1), axis=1)
L = L - L.median(axis=0) + L.median(axis=0).median()          # median normalization per sample
print("per-sample features detected: min %d, median %d, max %d" % tuple(np.percentile(tab.notna().sum(), [0, 50, 100])))

# Sample sheet from make_samples.py (<study>/samples.tsv, next to the table folder);
# falls back to MSV000084556-style names (<person>_<visit>) if there is none.
sheet_path = os.path.join(os.path.dirname(os.path.abspath(d)), "samples.tsv")
if os.path.exists(sheet_path):
    sheet = pd.read_csv(sheet_path, sep="\t", dtype=str, keep_default_na=False).set_index("sample")
    sheet = sheet.reindex(L.columns)
else:
    sheet = pd.DataFrame({"subject": [c.split("_")[0] for c in L.columns],
                          "timepoint": [c.split("_")[1] for c in L.columns],
                          "matrix": "", "role": "subject"}, index=L.columns)
subj, role = sheet.subject.fillna(""), sheet.role.fillna("subject")
Z = L.sub(L.mean(axis=1), axis=0).div(L.std(axis=1) + 1e-9, axis=0)
C = np.corrcoef(Z.T.values)

# 0) Technical precision: pooled QC injections should barely vary
if (role == "pool").sum() >= 3:
    raw = X.loc[:, (role == "pool").values]
    cv = raw.std(axis=1) / raw.mean(axis=1)
    print(f"pooled QC ({(role == 'pool').sum()} runs): median feature CV {cv.median():.0%}, "
          f"{(cv < 0.3).mean():.0%} of features under 30%")

# 1) Do a person's own samples look more alike than other people's? (people with >=2 samples only)
people = (role == "subject").values
s = np.asarray(subj.values, dtype=object)
same = (s[:, None] == s[None, :]) & people[:, None] & people[None, :]
if same.sum() > people.sum():
    iu = np.triu_indices_from(C, 1)
    both = (people[:, None] & people[None, :])[iu]
    w, b = C[iu][same[iu]], C[iu][both & ~same[iu]]
    print(f"\nsample-sample correlation: same person {w.mean():.2f} (n={len(w)}), different people {b.mean():.2f} (n={len(b)})")
    idx = np.where(people)[0]
    repeat = [i for i in idx if same[i].sum() > 1]
    Cp = np.where(np.eye(len(C)) == 1, -9, C)[:, idx]
    nn_same = np.mean([s[idx[np.argmax(Cp[i])]] == s[i] for i in repeat])
    chance = np.mean([(same[i].sum() - 1) / (len(idx) - 1) for i in repeat])
    print(f"nearest neighbour is the same person for {nn_same:.0%} of {len(repeat)} samples (chance ~{chance:.0%})")
else:
    print("\nwithin-person check skipped: no person has more than one sample")

# 2) PCA
U, S, Vt = np.linalg.svd(Z.T.values - Z.T.values.mean(0), full_matrices=False)
ev = S**2 / (S**2).sum()
print(f"PCA variance explained PC1-5: {np.round(ev[:5]*100,1)}")
pcs = pd.DataFrame(U[:, :2] * S[:2], index=L.columns, columns=["PC1", "PC2"])
tp = sheet.timepoint.fillna("")
for label, g in [("timepoint", tp), ("matrix", sheet.matrix.fillna("")), ("role", role)]:
    if g.nunique() > 1:
        print(f"mean PC1/PC2 by {label}:\n", pcs.groupby(g.values).mean().round(1).to_string())

# 3) Known plasma compounds ([M+H]+), within 5 ppm
known = {"creatinine": 114.0662, "proline": 116.0706, "betaine": 118.0863, "leucine/isoleucine": 132.1019,
         "carnitine": 162.1125, "phenylalanine": 166.0863, "hippuric acid": 180.0655, "tyrosine": 182.0812,
         "caffeine": 195.0877, "tryptophan": 205.0972, "cortisol": 363.2166, "LPC 16:0": 496.3398,
         "LPC 18:2": 520.3398, "PC 34:2": 758.5694}
print("\nknown compounds (feature detected in n samples, RT s):")
for name, mz in known.items():
    hit = info[(abs(info.mz - mz) / mz * 1e6) < 5].sort_values("n_detected", ascending=False)
    if len(hit):
        h = hit.iloc[0]; print(f"  {name:20s} m/z {mz:.4f}: n={int(h.n_detected):3d}, RT {h.rt_s:.0f}")
    else:
        print(f"  {name:20s} m/z {mz:.4f}: not found")
L.to_csv(os.path.join(d, f"feature_table_filtered_log2.tsv"), sep="\t")
pcs.assign(subject=subj, timepoint=tp).to_csv(os.path.join(d, "pca_scores.tsv"), sep="\t")
