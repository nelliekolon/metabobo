#!/usr/bin/env python3
"""Step 2b: pooled-QC batch correction, for studies whose pools span several batches.

usage: python3 batch_correct.py <table_dir> <out_dir>
Reads <table_dir>/feature_table.tsv and ../samples.tsv (columns role, batch). For every feature
and batch, on log2 intensities:
    offset[batch] = median(pools in that batch) - median(all pools)
and every run in the batch (pools and people) is shifted by -offset. A batch with no pools, or
a feature detected in fewer than 2 of a batch's pools, borrows the offset from the nearest
batch in time that has one; if none has, that feature is left uncorrected in that batch.
Writes out_dir/feature_table.tsv (corrected, linear scale, NaN kept as NaN) and copies
feature_info.tsv, so check_table.py runs on out_dir unchanged. Also writes batch_report.txt.

Validation is held out: offsets are fitted on half of each batch's pools and the pool CV is
measured on the other half, so the reported improvement is not graded on the pools it was fitted to.
"""
import os, shutil, sys
import numpy as np, pandas as pd

MIN_POOLS = 2


def offsets(L, pools, batch):
    """feature x batch offsets (log2), fitted on the given pool columns."""
    ref = L[pools].median(axis=1)
    out = {}
    for b in sorted(batch.unique()):
        cols = [c for c in pools if batch[c] == b]
        if len(cols) == 0:
            continue
        sub = L[cols]
        off = sub.median(axis=1) - ref
        off[sub.notna().sum(axis=1) < MIN_POOLS] = np.nan
        out[b] = off
    return pd.DataFrame(out)


def fill_nearest(O, batches):
    """For each batch (including ones with no pools), use its own offset, else the nearest batch's."""
    days = {b: pd.Timestamp(b) for b in batches}
    filled = {}
    for b in batches:
        order = sorted(O.columns, key=lambda x: abs((days[b] - pd.Timestamp(x)).days))
        col = pd.Series(np.nan, index=O.index)
        for x in order:
            col = col.fillna(O[x])
        filled[b] = col.fillna(0.0)
    return pd.DataFrame(filled)


def apply(L, O, batch):
    C = L.copy()
    for b in batch.unique():
        cols = batch.index[batch == b]
        C[cols] = L[cols].sub(O[b], axis=0)
    return C


def cv(M):
    lin = 2 ** M
    return (lin.std(axis=1) / lin.mean(axis=1)).median()


def main():
    tdir, out_dir = sys.argv[1], sys.argv[2]
    tab = pd.read_csv(os.path.join(tdir, "feature_table.tsv"), sep="\t", index_col=0)
    sheet = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(tdir)), "samples.tsv"),
                        sep="\t", dtype=str, keep_default_na=False).set_index("sample").reindex(tab.columns)
    batch, role = sheet.batch, sheet.role
    if (batch == "").any():
        sys.exit("every run needs a batch in samples.tsv")
    pools = list(tab.columns[role == "pool"])
    L = np.log2(tab)
    rep = []

    # held-out validation: fit on alternate pools within each batch, score the others
    fit, test = [], []
    for b in sorted(batch[pools].unique()):
        cols = sorted(c for c in pools if batch[c] == b)
        fit += cols[0::2]; test += cols[1::2]
    O_fit = fill_nearest(offsets(L, fit, batch), sorted(batch.unique()))
    detected = L[test].notna().mean(axis=1) >= 0.8
    before = cv(L.loc[detected, test])
    after = cv(apply(L, O_fit, batch).loc[detected, test])
    rep.append(f"held-out pools (n={len(test)}, offsets fitted on the other {len(fit)}): "
               f"median CV {before:.0%} -> {after:.0%}")

    # final correction uses all pools
    O_all = offsets(L, pools, batch)
    no_pool = sorted(set(batch.unique()) - set(O_all.columns))
    O = fill_nearest(O_all, sorted(batch.unique()))
    C = apply(L, O, batch)
    keep = L[pools].notna().mean(axis=1) >= 0.8
    rep.append(f"all pools (n={len(pools)}), in-sample: median CV {cv(L.loc[keep, pools]):.0%} -> "
               f"{cv(C.loc[keep, pools]):.0%}")
    rep.append(f"batches: {batch.nunique()}, with pools: {O_all.shape[1]}; "
               f"without pools (nearest-batch offset borrowed): {', '.join(no_pool) or 'none'}")
    flagged = sheet.index[(batch.isin(no_pool)) & (role == "subject")]
    rep.append(f"runs corrected with a borrowed offset: {len(flagged)} " + (f"({', '.join(flagged)})" if len(flagged) else ""))
    big = O_all.abs().median().round(2)
    rep.append("median |offset| per batch (log2): " + ", ".join(f"{k} {v}" for k, v in big.items()))

    os.makedirs(out_dir, exist_ok=True)
    (2 ** C).to_csv(os.path.join(out_dir, "feature_table.tsv"), sep="\t")
    shutil.copy(os.path.join(tdir, "feature_info.tsv"), os.path.join(out_dir, "feature_info.tsv"))
    O.to_csv(os.path.join(out_dir, "batch_offsets_log2.tsv"), sep="\t")
    with open(os.path.join(out_dir, "batch_report.txt"), "w") as fh:
        fh.write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
