#!/usr/bin/env python3
"""Eval step 5: sort every wrong raw-field answer by where the error came from.

usage: python3 classify_errors.py <scored.tsv> <runs_dir> <inputs_dir> <out.tsv> [hand_labels.tsv]

A wrong answer is one of:
  text_wrong       the quote states the model's value and the true value appears nowhere in the
                   text: the model copied the text faithfully; the text disagrees with the raw file
  attribution      the text states the true value too (or the study ran both ion modes): the model
                   picked the setting of another acquisition / reported the study, not the file
  misread          the quote does not state the value the model gave
  review           none of the rules decide it; needs a hand check
Each row also keeps the quote, so the auto category can be checked and overridden by hand
(hand_category / note columns). hand_labels.tsv (target, field, value_regex, category, note)
fills them; the first row whose regex matches the answer wins. Hand categories add:
  text_vague       the text gives a less specific or several-methods value (Q Exactive for Q Exactive
                   Plus; "covered 70-1000 m/z" across methods whose windows differ)
  definition       run length: the text gives the gradient program, the file recorded a different span
  task_ambiguous   no raw file name was given (Workbench targets) and the study ran both ion modes
"""
import csv, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score import norm, number, instrument   # same normalization as the scorer

FIELDS = {"ion_mode", "scan_range_low", "scan_range_high", "run_length_min", "instrument_model"}


def nums(s):
    return {float(x.replace(",", "")) for x in re.findall(r"(?<![\w.])\d[\d,]*(?:\.\d+)?", s)}


def states(field, value, text):
    """does this text state this value?"""
    t = norm(text)
    if field == "ion_mode":
        v = norm(value)
        if v == "both":
            return bool(re.search(r"both|positive and negative|negative and positive|switching", t))
        return v[:3] in t                                  # "pos"/"neg" (positive, negative, ESI+ is missed)
    if field == "instrument_model":
        a = instrument(value)
        return bool(a) and (a in instrument(text) or a.replace(" ", "") in instrument(text).replace(" ", ""))
    x = number(value)
    return x is not None and any(abs(x - n) < 1e-6 for n in nums(t))


def main():
    scored, runs, inputs, out = sys.argv[1:5]
    rows = []
    for r in csv.DictReader(open(scored), delimiter="\t"):
        if r["field"] not in FIELDS or r["correct"] != "False":
            continue
        rec = json.load(open(os.path.join(runs, r["model"], f"{r['target']}__{r['condition']}.json")))
        quote = rec["answer"][r["field"]]["quote"]
        text = open(os.path.join(inputs, rec["study"], f"{r['condition']}.txt")).read()
        f, v, tr = r["field"], r["value"], r["truth"]
        truth_val = tr if f in ("ion_mode", "instrument_model") else str(int(round(float(tr))))
        q_states = states(f, v, quote)
        t_has_truth = states(f, truth_val, text)
        if f == "ion_mode" and norm(v) == "both":
            cat = "attribution"
        elif not q_states:
            cat = "misread"
        elif t_has_truth:
            cat = "attribution"
        else:
            cat = "text_wrong"
        rows.append(dict(model=r["model"], target=r["target"], condition=r["condition"], field=f, value=v,
                         truth=tr, auto_category=cat, quote_states_value=q_states, text_has_truth=t_has_truth,
                         quote=re.sub(r"\s+", " ", quote)[:300], hand_category="", note=""))
    if len(sys.argv) > 5:
        labels = list(csv.DictReader(open(sys.argv[5]), delimiter="\t"))
        for r in rows:
            hit = next((l for l in labels if l["target"] == r["target"] and l["field"] == r["field"]
                        and re.search(l["value_regex"], r["value"])), None)
            r["hand_category"], r["note"] = (hit["category"], hit["note"]) if hit else ("UNLABELED", "")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    from collections import Counter
    for m in sorted({r["model"] for r in rows}):
        print(m, "auto:", dict(Counter(r["auto_category"] for r in rows if r["model"] == m)))
        print(m, "hand:", dict(Counter(r["hand_category"] for r in rows if r["model"] == m)))


if __name__ == "__main__":
    main()
