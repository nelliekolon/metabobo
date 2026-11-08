#!/usr/bin/env python3
"""Eval step 1: measure method metadata from one raw-derived mzML per study (the answer key's
raw-file fields; see SPEC.md).

usage: python3 extract_truth.py <files.tsv> <out.tsv>
files.tsv: columns study, mzml (one representative run per study, converted with ALL MS levels).
Writes one row per study: ion_mode, scan_range_low, scan_range_high, run_length_min,
instrument_model, ms2_acquired, plus how each was measured (scan_range_source) and counts.

scan range: the MS1 scan window recorded in the file; if the file has none (some MassIVE files
were converted through mzXML), the lowest/highest observed m/z over all MS1 spectra, marked
'observed (approx)': observed peaks sit inside the window, so this can only understate it.
"""
import csv, os, sys
import pyopenms as oms

POL = {1: "positive", 2: "negative"}


def measure(path):
    e = oms.MSExperiment(); opt = oms.PeakFileOptions(); opt.setFillData(False)
    fh = oms.MzMLFile(); fh.setOptions(opt); fh.load(path, e)
    sp = e.getSpectra()
    ms1 = [s for s in sp if s.getMSLevel() == 1]
    pols = {POL.get(s.getInstrumentSettings().getPolarity(), "unknown") for s in sp} - {"unknown"}
    wins = [(w.begin, w.end) for s in ms1 for w in s.getInstrumentSettings().getScanWindows()]
    if wins:
        low, high, src = min(w[0] for w in wins), max(w[1] for w in wins), "scan window"
        n_windows = len({(round(a), round(b)) for a, b in wins})
    else:                                  # no window recorded: load MS1 peaks, take observed extremes
        e2 = oms.MSExperiment(); opt2 = oms.PeakFileOptions(); opt2.setMSLevels([1])
        fh2 = oms.MzMLFile(); fh2.setOptions(opt2); fh2.load(path, e2)
        mz = [p for s in e2.getSpectra() for p in (s.get_peaks()[0][[0, -1]] if s.size() else [])]
        low, high, src, n_windows = min(mz), max(mz), "observed (approx)", 0
    return dict(
        ion_mode="both" if len(pols) > 1 else (pols.pop() if pols else "unknown"),
        scan_range_low=round(low, 1), scan_range_high=round(high, 1), scan_range_source=src,
        distinct_scan_windows=n_windows,
        run_length_min=round(max(s.getRT() for s in sp) / 60, 2),
        instrument_model=e.getInstrument().getName() or "unknown",
        ms2_acquired="yes" if any(s.getMSLevel() >= 2 for s in sp) else "no",
        n_spectra=len(sp), n_ms1=len(ms1))


def main():
    files, out = sys.argv[1], sys.argv[2]
    rows = []
    with open(files) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = measure(r["mzml"])
            rows.append(dict(study=r["study"], mzml=os.path.basename(r["mzml"]), **m))
            print(rows[-1], flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
