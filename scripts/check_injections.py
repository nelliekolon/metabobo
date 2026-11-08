#!/usr/bin/env python3
"""Compare repeat injections of the same sample (MTBLS2469 runs each sample as _P_2, _P_3, _P_4).

usage: python3 check_injections.py <mzML_dir> <feature_dir> <out.tsv>
For every sample with all three injections: number of spectra per injection, number of
features, and how well injections agree (features matched at 10 ppm / 15 s, Pearson r of
log2 intensities). If _P_2 is a different acquisition, r(3,2) will sit below r(3,4) and its
spectrum count will differ, across most samples, not just one.
"""
import glob, os, re, sys
import numpy as np, pandas as pd
import pyopenms as oms


def n_spectra(path):
    # count <spectrum ...> tags without loading the file
    n = 0
    with open(path, "rb") as fh:
        for line in fh:
            if b"<spectrum " in line:
                n += 1
    return n


def load(path):
    fm = oms.FeatureMap(); oms.FeatureXMLFile().load(path, fm)
    return np.array([(f.getMZ(), f.getRT(), f.getIntensity()) for f in fm])


def agree(a, b):
    b = b[np.argsort(b[:, 0])]
    pairs = []
    for m, r, i in a:
        lo, hi = np.searchsorted(b[:, 0], [m * (1 - 1e-5), m * (1 + 1e-5)])
        c = b[lo:hi]; c = c[np.abs(c[:, 1] - r) < 15]
        if len(c):
            pairs.append((i, c[np.argmax(c[:, 2]), 2]))
    p = np.log2(np.array(pairs))
    return len(p), np.corrcoef(p.T)[0, 1]


def main():
    mz_dir, feat_dir, out = sys.argv[1:4]
    groups = {}
    for f in glob.glob(os.path.join(feat_dir, "*.featureXML")):
        m = re.match(r"(.+)_P_(\d)$", os.path.basename(f)[:-11])
        if m:
            groups.setdefault(m[1], {})[m[2]] = f
    rows = []
    for base, inj in sorted(groups.items()):
        if set(inj) != {"2", "3", "4"}:
            continue
        F = {k: load(v) for k, v in inj.items()}
        n34, r34 = agree(F["3"], F["4"]); n32, r32 = agree(F["3"], F["2"])
        spec = {k: n_spectra(os.path.join(mz_dir, os.path.basename(v)[:-11] + ".mzML")) for k, v in inj.items()}
        rows.append(dict(sample=base, kind="pool" if "QCPool" in base else "subject",
                         spectra_2=spec["2"], spectra_3=spec["3"], spectra_4=spec["4"],
                         features_2=len(F["2"]), features_3=len(F["3"]), features_4=len(F["4"]),
                         matched_34=n34, r_34=round(r34, 3), matched_32=n32, r_32=round(r32, 3)))
    d = pd.DataFrame(rows); d.to_csv(out, sep="\t", index=False)
    print(f"{len(d)} samples with injections 2, 3 and 4")
    for kind, g in d.groupby("kind"):
        print(f"\n{kind} (n={len(g)}):")
        print(f"  spectra    inj2 median {g.spectra_2.median():.0f}, inj3 {g.spectra_3.median():.0f}, inj4 {g.spectra_4.median():.0f}")
        print(f"  features   inj2 median {g.features_2.median():.0f}, inj3 {g.features_3.median():.0f}, inj4 {g.features_4.median():.0f}")
        print(f"  agreement  r(3,4) median {g.r_34.median():.2f}, r(3,2) median {g.r_32.median():.2f}; "
              f"r(3,2) < r(3,4) in {(g.r_32 < g.r_34).mean():.0%} of samples")
        print(f"  inj2 has fewer spectra than inj3 in {(g.spectra_2 < g.spectra_3).mean():.0%} of samples")


if __name__ == "__main__":
    main()
