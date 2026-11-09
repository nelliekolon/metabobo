#!/usr/bin/env python3
"""Figure data: pooled-QC drift by run day, before and after batch correction (held out).

usage: python3 batch_drift.py <table_dir> <out.tsv>
Same split as batch_correct.py: offsets are fitted on alternate pools in each batch; this writes,
for every HELD-OUT pool, its median log2 deviation from the per-feature median of the OTHER
held-out pools (leave-one-out, so a pool is never compared with itself; features detected in
>= 80% of held-out pools), before and after applying those offsets.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from batch_correct import offsets, fill_nearest, apply


def main():
    tdir, out = sys.argv[1:3]
    tab = pd.read_csv(os.path.join(tdir, "feature_table.tsv"), sep="\t", index_col=0)
    sheet = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(tdir)), "samples.tsv"), sep="\t",
                        dtype=str, keep_default_na=False).set_index("sample").reindex(tab.columns)
    batch, pools = sheet.batch, list(tab.columns[sheet.role == "pool"])
    L = np.log2(tab)
    fit, test = [], []
    for b in sorted(batch[pools].unique()):
        cols = sorted(c for c in pools if batch[c] == b)
        fit += cols[0::2]; test += cols[1::2]
    O = fill_nearest(offsets(L, fit, batch), sorted(batch.unique()))
    det = L[test].notna().mean(axis=1) >= 0.8
    rows = []
    for stage, M in (("before", L), ("after", apply(L, O, batch))):
        T = M.loc[det, test]
        for c in test:
            dev = (T[c] - T.drop(columns=c).median(axis=1)).median()
            rows.append(dict(stage=stage, pool=c, batch=batch[c], deviation_log2=round(float(dev), 3)))
    pd.DataFrame(rows).sort_values(["stage", "batch", "pool"], ascending=[False, True, True]).to_csv(out, sep="\t", index=False)
    print(pd.DataFrame(rows).groupby(["stage", "batch"]).deviation_log2.median())


if __name__ == "__main__":
    main()
