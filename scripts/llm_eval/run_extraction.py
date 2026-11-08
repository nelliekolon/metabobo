#!/usr/bin/env python3
"""Eval step 3: ask a model to extract method metadata from one text condition (see SPEC.md).

usage: python3 run_extraction.py <inputs_dir> <targets.tsv> <runs_dir> <model> [<model> ...]
targets.tsv: study, target, raw_file. One row per acquisition to describe (a study with two
acquisitions, e.g. two columns, has two rows). Conditions run are whichever of methods.txt,
fulltext.txt, repository.txt exist in <inputs_dir>/<study>/.

Each call is saved to <runs_dir>/<model>/<target>__<condition>.json (request settings, raw
response text, parsed answer, usage, which model actually served it). Existing results are
skipped, so the script is safe to rerun. Reads ANTHROPIC_API_KEY from .env if not set.
"""
import csv, json, os, sys, time
import anthropic

FIELDS = {
    "ion_mode": "positive, negative or both",
    "scan_range_low": "lowest m/z of the full-scan (MS1) range, as a number",
    "scan_range_high": "highest m/z of the full-scan (MS1) range, as a number",
    "run_length_min": "length of one LC-MS run (acquisition or gradient) in minutes, as a number",
    "instrument_model": "mass spectrometer model, e.g. 'Q Exactive HF'",
    "ms2_acquired": "yes or no: were fragmentation (MS2/MS/MS) spectra recorded in this run",
    "chromatography": "reverse phase, HILIC or other",
    "column": "LC column name and dimensions",
    "sample_type": "plasma, serum, whole blood, or mixed",
    "group_definitions": "what the study's sample groups or timepoints are, in one sentence",
}
SCHEMA = {
    "type": "object",
    "properties": {f: {"type": "object",
                       "properties": {"value": {"type": "string"},
                                      "quote": {"type": "string"},
                                      "confidence": {"type": "string", "enum": ["low", "medium", "high"]}},
                       "required": ["value", "quote", "confidence"], "additionalProperties": False}
                   for f in FIELDS},
    "required": list(FIELDS), "additionalProperties": False,
}
SYSTEM = """You extract analytical-method metadata for one LC-MS metabolomics acquisition from the text you are given.

Rules:
- Use only the text provided. Do not use outside knowledge about the study, lab, or instrument defaults.
- The text may describe several experiments (other omics layers, other columns or ion modes). Report values only for the acquisition that produced the named raw file. If the text does not make clear which settings belong to that acquisition, answer "not stated".
- If a field is not stated for that acquisition, set value to "not stated", quote to "", and confidence to "high" if you are sure it is absent.
- quote must be copied verbatim from the text (an exact substring, at most about 300 characters) and must contain the evidence for value.
- Numbers: give just the number in value (e.g. "70", "1050", "15")."""

FALLBACK_MODELS = {"claude-opus-5-5"}   # models that accept server-side refusal fallback


def load_key():
    if os.environ.get("ANTHROPIC_API_KEY"):
        return
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".env")
    for line in open(here):
        k, _, v = line.strip().partition("=")
        if k == "ANTHROPIC_API_KEY":
            os.environ[k] = v


def user_prompt(raw_file, condition, text):
    fields = "\n".join(f"- {f}: {d}" for f, d in FIELDS.items())
    label = {"methods": "the paper's methods section", "fulltext": "the full paper",
             "repository": "the data repository's description of the study"}[condition]
    return (f"Raw file: {raw_file}\n\nFields:\n{fields}\n\n"
            f"Text ({label}):\n<text>\n{text}\n</text>")


def call(client, model, prompt):
    kwargs = dict(model=model, max_tokens=16000, system=SYSTEM,
                  messages=[{"role": "user", "content": prompt}],
                  output_config={"effort": "medium", "format": {"type": "json_schema", "schema": SCHEMA}})
    if model in FALLBACK_MODELS:
        # Refusal fallback on by default; any answer served by another model is flagged below
        # and excluded from this model's scores.
        return client.beta.messages.create(betas=["server-side-fallback-2026-07-01"], fallbacks="default", **kwargs)
    return client.messages.create(**kwargs)


def main():
    inputs, targets, runs = sys.argv[1:4]
    models = sys.argv[4:]
    load_key()
    client = anthropic.Anthropic()
    for model in models:
        os.makedirs(os.path.join(runs, model), exist_ok=True)
        for t in csv.DictReader(open(targets), delimiter="\t"):
            for cond in ["methods", "fulltext", "repository"]:
                src = os.path.join(inputs, t["study"], f"{cond}.txt")
                out = os.path.join(runs, model, f"{t['target']}__{cond}.json")
                if not os.path.exists(src) or os.path.getsize(src) == 0 or os.path.exists(out):
                    continue
                text = open(src).read()
                rec = dict(study=t["study"], target=t["target"], raw_file=t["raw_file"], condition=cond,
                           model_requested=model, effort="medium", input_chars=len(text))
                t0 = time.time()
                try:
                    r = call(client, model, user_prompt(t["raw_file"], cond, text))
                except anthropic.BadRequestError as e:
                    print(f"{model} {t['target']} {cond}: BAD REQUEST {e.message}", flush=True); continue
                except anthropic.APIStatusError as e:
                    print(f"{model} {t['target']} {cond}: API ERROR {e.status_code} {e.message}", flush=True); continue
                rec["seconds"] = round(time.time() - t0, 1)
                rec["stop_reason"] = r.stop_reason
                rec["model_served"] = r.model
                iters = getattr(r.usage, "iterations", None) or []
                rec["fallback_used"] = any(getattr(i, "type", "") == "fallback_message" for i in iters)
                rec["usage"] = {"input": r.usage.input_tokens, "output": r.usage.output_tokens}
                rec["request_id"] = r._request_id
                text_out = next((b.text for b in r.content if b.type == "text"), "")
                rec["raw_text"] = text_out
                try:
                    rec["answer"] = json.loads(text_out) if r.stop_reason != "refusal" else None
                except json.JSONDecodeError:
                    rec["answer"] = None
                json.dump(rec, open(out, "w"), indent=1)
                print(f"{model} {t['target']:15s} {cond:10s} {r.stop_reason:9s} in {r.usage.input_tokens:6d} "
                      f"out {r.usage.output_tokens:5d} {rec['seconds']:5.1f}s"
                      + ("  [FALLBACK]" if rec["fallback_used"] else ""), flush=True)


if __name__ == "__main__":
    main()
