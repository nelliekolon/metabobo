#!/usr/bin/env python3
"""Eval step 4: score extraction runs (see SPEC.md).

usage: python3 score.py <runs_dir> <inputs_dir> <truth.tsv[,truth2.tsv]> <targets.tsv> <out_dir>
Truth files are keyed by target (or study); a blank truth value is not scored.
Scores the raw-file fields automatically against the measured truth, checks every quote
against the input text, and writes per-answer rows (scored.tsv) plus a summary (score_report.txt).

Per answer:
  stated         value is not "not stated"
  quote_found    the quote occurs verbatim in the input (after Unicode/whitespace normalization)
  unsupported    stated, but the quote is empty or not found: the main hallucination measure
  correct        (raw-file fields only) value matches the measured truth within tolerance
  name_leak      the raw file name itself gives the answer away (e.g. _POS for ion mode)
"""
import csv, glob, json, os, re, sys, unicodedata
from collections import Counter, defaultdict

RAW_FIELDS = ["ion_mode", "scan_range_low", "scan_range_high", "run_length_min", "instrument_model", "ms2_acquired"]
# ms2_acquired is reported but not scored: one raw file shows whether THAT run recorded MS2, while
# papers describe the study (MTBLS2469: MS2 only in injection 2), so the two answer different questions
SCORED = [f for f in RAW_FIELDS if f != "ms2_acquired"]
NOT_STATED = {"not stated", "", "n/a", "unknown", "not reported"}


def norm(s):
    s = unicodedata.normalize("NFKC", s or "")
    s = s.replace(" ", " ").replace(" ", " ").replace(" ", " ")
    s = re.sub(r"[‐-―−]", "-", s)
    return re.sub(r"\s+", " ", s).strip().lower()


def number(v):
    m = re.search(r"\d[\d ,]*\.?\d*", (v or "").replace(" ", "").replace(" ", ""))
    return float(re.sub(r"[ ,]", "", m.group())) if m else None


def instrument(v):
    v = norm(v)
    for w in ["thermo fisher scientific", "thermo scientific", "thermo", "orbitrap", "mass spectrometer", "(", ")", "-", "™", "®"]:
        v = v.replace(w, " ")
    return re.sub(r"\s+", " ", v).strip()


def correct(field, value, truth):
    if field == "ion_mode":
        return norm(value) == truth
    if field == "ms2_acquired":
        return norm(value)[:3] in ("yes", "no") and norm(value).startswith(truth)
    if field == "instrument_model":
        # the mzML records the model family (e.g. "Q Exactive"); a more specific answer
        # in that family ("Q Exactive HF") counts as compatible, not provably right
        a, t = instrument(value), instrument(truth)
        return a == t or a.startswith(t + " ")
    x, t = number(value), float(truth)
    if x is None:
        return False
    if field == "scan_range_low":
        return abs(x - t) <= 5
    if field == "scan_range_high":
        return abs(x - t) <= 0.05 * t
    if field == "run_length_min":
        return abs(x - t) <= max(1.0, 0.10 * t)


def leak(field, raw_file):
    f = raw_file.upper()
    return ((field == "ion_mode" and re.search(r"(^|_)(POS|NEG)(_|\.|$)", f) is not None)
            or (field in ("chromatography", "column") and re.search(r"C18|HILIC|T3", f) is not None)
            or (field == "sample_type" and re.search(r"PLASMA|SERUM", f) is not None))


def main():
    runs, inputs, truth_f, targets_f, out_dir = sys.argv[1:6]
    truth = {}
    for tf in truth_f.split(","):
        for r in csv.DictReader(open(tf), delimiter="\t"):
            truth[r["study"]] = r
    rows = []
    for path in sorted(glob.glob(os.path.join(runs, "*", "*.json"))):
        rec = json.load(open(path))
        if rec.get("answer") is None:
            rows.append(dict(model=rec["model_requested"], target=rec["target"], condition=rec["condition"],
                             field="(all)", value="", stated=False, quote_found=False, unsupported=False,
                             correct="", name_leak=False, excluded="refusal or unparsable"))
            continue
        src = norm(open(os.path.join(inputs, rec["study"], f"{rec['condition']}.txt")).read())
        t = truth.get(rec["target"]) or truth.get(rec["study"], {})
        for field, a in rec["answer"].items():
            value, quote = a.get("value", ""), a.get("quote", "")
            stated = norm(value) not in NOT_STATED
            found = bool(norm(quote)) and norm(quote) in src
            row = dict(model=rec["model_requested"], target=rec["target"], condition=rec["condition"], field=field,
                       value=value, confidence=a.get("confidence", ""), stated=stated, quote_found=found,
                       unsupported=stated and not found, correct="", truth="",
                       name_leak=leak(field, rec["raw_file"]),
                       excluded="served by fallback model" if rec.get("fallback_used") else "")
            if field in SCORED and t and str(t.get(field, "")).strip() not in ("", "unknown"):
                row["truth"] = t[field]
                row["correct"] = correct(field, value, t[field]) if stated else ""
            elif field == "ms2_acquired" and t:
                row["truth"] = t.get(field, "")
            rows.append(row)
    os.makedirs(out_dir, exist_ok=True)
    keys = ["model", "target", "condition", "field", "value", "confidence", "truth", "stated", "quote_found",
            "unsupported", "correct", "name_leak", "excluded"]
    with open(os.path.join(out_dir, "scored.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, keys, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)

    use = [r for r in rows if not r["excluded"] and r["field"] != "(all)"]
    rep = [f"{len(rows)} answers scored ({sum(1 for r in rows if r['excluded'])} excluded)\n"]
    # majority baseline on raw-file fields: always answer the most common true value
    rep.append("majority baseline (most common true value across the pilot acquisitions):")
    for f in SCORED:
        vals = [truth[k][f] for k in truth if str(truth[k].get(f, "")).strip() not in ("", "unknown")]
        if not vals:
            continue
        top, n = Counter(vals).most_common(1)[0]
        rep.append(f"  {f:17s} always '{top}': {n}/{len(vals)} correct")
    for model in sorted({r["model"] for r in use}):
        rep.append(f"\n=== {model}")
        rep.append(f"  {'condition':10s} {'stated':>8s} {'unsupported':>12s} {'raw-field acc':>14s} {'(excl. name leaks)':>19s}")
        for cond in ["methods", "fulltext", "repository"]:
            g = [r for r in use if r["model"] == model and r["condition"] == cond]
            if not g:
                continue
            st = [r for r in g if r["stated"]]
            raw = [r for r in st if r["correct"] != ""]
            raw_nl = [r for r in raw if not r["name_leak"]]
            acc = lambda xs: f"{sum(r['correct'] for r in xs)}/{len(xs)}" if xs else "-"
            rep.append(f"  {cond:10s} {len(st):3d}/{len(g):<4d} {sum(r['unsupported'] for r in g):6d}/{len(st):<5d} "
                       f"{acc(raw):>14s} {acc(raw_nl):>19s}")
        wrong = [r for r in use if r["model"] == model and r["correct"] is False]
        if wrong:
            rep.append("  wrong raw-field answers:")
            for r in wrong:
                rep.append(f"    {r['target']:14s} {r['condition']:10s} {r['field']:17s} said {r['value']!r:12s} "
                           f"truth {r['truth']!r}  (quote found: {r['quote_found']}, confidence {r['confidence']})")
        uns = [r for r in use if r["model"] == model and r["unsupported"]]
        if uns:
            rep.append("  unsupported answers (quote not in the text):")
            for r in uns:
                rep.append(f"    {r['target']:14s} {r['condition']:10s} {r['field']:17s} said {r['value'][:50]!r}")
    open(os.path.join(out_dir, "score_report.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
