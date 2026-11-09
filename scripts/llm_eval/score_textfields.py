#!/usr/bin/env python3
"""Eval step 6b: score the text-only fields (chromatography, column, sample_type) against the
hand-checked answer key, and sort every wrong answer by cause with hand labels.

usage: python3 score_textfields.py <runs_dir> <truth_textfields.tsv> <hand_labels_text.tsv> <out_dir>

Key (truth_textfields.tsv): one row per target and field, drafted from the repositories'
structured submitter fields (fetch_structured.py), then corrected by hand where they were wrong
(note column says how). A value of "uncertain" is not scored.
Answers are normalized before comparison: chromatography to RP / HILIC; column to its chemistry
(BEH C18, HSS T3, BEH Amide, ...; dimensions are not scored); sample type to plasma / serum /
whole blood / mixed / urine.
Wrong-answer categories (hand labels): text_wrong (the text given states the model's value and
it is wrong), text_conflict (the paper and the repository's structured fields disagree and the
model followed its text), attribution (the model reported another acquisition, or the study as a
whole, where the text cannot say which applies to this file), misread (the text does not support
the answer), prompt_ambiguous (the field's options did not fit, e.g. "mixed" for plasma + stool).
"""
import csv, glob, json, os, re, sys
from collections import Counter, defaultdict

CHEM = [("BEH Amide", r"beh\s*amide|amide.*beh"), ("Accucore Amide", r"accucore.*amide"), ("CSH C18", r"csh"),
        ("HSS T3", r"hss\s*t3|\bt3\b"), ("BEH C18", r"beh\s*c18|beh\b"), ("Hypersil GOLD", r"hypersil\s*gold"),
        ("Zorbax Eclipse Plus C18", r"zorbax|eclipse"), ("Kinetex C18", r"kinetex"), ("CORTECS C18", r"cortecs"),
        ("Dr. Maisch C18", r"maisch|reprosil"), ("HILIC", r"hilic")]


def stated(s):
    return bool(s) and "not stated" not in s.lower()


def column(s):
    s = (s or "").lower()
    return next((k for k, p in CHEM if re.search(p, s)), "other") if stated(s) else ""


def chromatography(s):
    s = (s or "").lower()
    if not stated(s):
        return ""
    return "HILIC" if "hilic" in s else "RP" if re.search(r"revers|\brp\b|c18", s) else "other"


def sample_type(s):
    s = (s or "").lower()
    if not stated(s):
        return ""
    for k, p in [("urine", r"urine"), ("whole blood", r"whole"), ("mixed", r"mixed|plasma.*serum|serum.*plasma"),
                 ("plasma", r"plasma"), ("serum", r"serum")]:
        if re.search(p, s):
            return k
    return "other"


NORM = {"chromatography": chromatography, "column": column, "sample_type": sample_type}


def main():
    runs, key_f, labels_f, out_dir = sys.argv[1:5]
    key = {(r["target"], r["field"]): r for r in csv.DictReader(open(key_f), delimiter="\t")}
    labels = {(r["target"], r["field"], r["answer"]): r for r in csv.DictReader(open(labels_f), delimiter="\t")}
    rows = []
    for path in sorted(glob.glob(os.path.join(runs, "*", "*.json"))):
        rec = json.load(open(path))
        if not rec.get("answer") or (rec["target"], "column") not in key or rec.get("fallback_used"):
            continue
        for f, norm in NORM.items():
            k = key[(rec["target"], f)]
            if k["value"] == "uncertain":
                continue
            a = rec["answer"][f]
            v = norm(a["value"])
            row = dict(model=rec["model_requested"], target=rec["target"], condition=rec["condition"], field=f,
                       value=a["value"], answer=v, key=k["value"], stated=bool(v), correct="", category="", note="")
            if v:
                row["correct"] = v == k["value"]
                if not row["correct"]:
                    lab = labels.get((rec["target"], f, v), {})
                    row["category"], row["note"] = lab.get("category", "UNLABELED"), lab.get("note", "")
            rows.append(row)
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "scored_textfields.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

    rep = []
    for f in NORM:
        vals = [k["value"] for (t, ff), k in key.items() if ff == f and k["value"] != "uncertain"]
        top, n = Counter(vals).most_common(1)[0]
        rep.append(f"majority baseline {f:15s} always '{top}': {n}/{len(vals)}")
    for m in sorted({r["model"] for r in rows}):
        rep.append(f"\n=== {m}")
        for f in NORM:
            g = [r for r in rows if r["model"] == m and r["field"] == f]
            st = [r for r in g if r["stated"]]
            ok = sum(r["correct"] is True for r in st)
            rep.append(f"  {f:15s} stated {len(st):3d}/{len(g):<3d} correct {ok}/{len(st)} ({ok / len(st):.0%})  "
                       f"wrong by cause: {dict(Counter(r['category'] for r in st if r['correct'] is False))}")
        cats = Counter(r["category"] for r in rows if r["model"] == m and r["correct"] is False)
        rep.append(f"  all three fields, wrong by cause: {dict(cats)}")
    open(os.path.join(out_dir, "score_textfields_report.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))
    un = [r for r in rows if r["category"] == "UNLABELED"]
    if un:
        print("\nUNLABELED:", sorted({(r["target"], r["field"], r["answer"]) for r in un}))


if __name__ == "__main__":
    main()
