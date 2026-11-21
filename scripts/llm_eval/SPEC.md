# LLM metadata-extraction eval: spec (draft, 2026-11-08)

## Question
How accurately can an LLM extract analytical-method metadata for a public metabolomics study
from its paper, and how often does it state values the paper does not support?

Why it matters: missing or wrong method metadata is the bottleneck for combining public studies
(see the metabopool write-up: RT scales, scan ranges and sample types differ between studies and are
rarely recorded per sample). If LLMs can fill it in reliably, harmonization scales to hundreds of
studies; if they hallucinate plausible values, they make it worse in a way that is hard to see.

What makes this eval unusual: several answers can be checked against the data itself, not
against another annotator's reading of the same paper.

## Fields

Each field has one source of truth, chosen to be independent of the paper wherever possible.

| Field | Values | Truth from | Scoring |
| --- | --- | --- | --- |
| ion_mode | positive / negative / both | raw file (spectrum polarity) | exact |
| scan_range_low | m/z | raw file (scan window) | within 5 m/z |
| scan_range_high | m/z | raw file (scan window) | within 5% |
| run_length_min | minutes | raw file (last scan time) | within 1 min or 10% |
| instrument_model | e.g. Q Exactive HF-X | raw file (instrument metadata) | normalized exact |
| ms2_acquired | yes / no | raw file (MS levels present) | exact |
| chromatography | reverse phase / HILIC / other | repository field, checked by hand | exact |
| column | name and dimensions | repository field, checked by hand | manual: correct / partial / wrong |
| sample_type | plasma / serum / whole blood / mixed | repository field + sample sheet | exact |
| group_definitions | free text | paper, checked by hand | manual: correct / partial / wrong |

Every field may also be answered `not stated`.

A study with several acquisitions (e.g. positive and negative mode, or two columns) is scored
against the acquisition whose raw file was measured; the prompt says which.

## What the model gets
Conditions (same prompt and fields in each):
1. **methods**: the paper's methods section only.
2. **full text**: the whole open-access paper.
3. **repository**: the repository's own description and protocol text only (no paper).

Condition 3 separates "the LLM can read" from "the information is in the paper at all".

Output: JSON, one object per field: `value`, `quote` (verbatim supporting text from the input),
`confidence` (low / medium / high).

## Metrics
- **Accuracy when stated**: share of fields where the model gave a value and it is correct.
- **Coverage**: share of fields where the model gave a value rather than `not stated`.
- **Unsupported answers**: the model gave a value but its `quote` does not occur in the input
  (checked by string match after whitespace normalization), or occurs but does not contain the
  value. This is the main hallucination measure, and it needs no human judgement.
- **Wrong but quoted**: value is supported by the quote yet contradicts the raw file. These are
  cases where the paper itself disagrees with the data (or describes a different acquisition);
  reported separately, not counted as model errors.
- **Calibration**: accuracy by stated confidence.
- **Majority baseline**: for each field, the score of always answering the most common true value
  in the set (e.g. "positive", "Q Exactive"). A model only gets credit for beating it. The 4 pilot
  studies measured so far are all positive-mode Q Exactive, so on those fields the pilot cannot
  tell reading from guessing; the full set must be chosen to include negative mode, HILIC and
  other instruments (the survey's `instrument` column helps).

All metrics per field and overall, with bootstrap intervals over studies.

## Studies
- **Pilot (5)**: MSV000084556, MSV000085256, MSV000094097, MTBLS2469, MTBLS2542. Raw files
  already downloaded for 4; MTBLS2542 needs one file.
- **Full set (30-50)**: from `processed/paper_survey.tsv`: open-access paper + structured
  repository fields (76 candidates), after checking that the paper is the study's original paper
  (not a later reuse) and downloading one raw file per study (about 200-300 MB each).

## Answer key
Lives in `eval_private/` (gitignored) until results are final, so it cannot leak into the eval.
Raw-file fields are filled by a script; repository and paper fields by hand, each with the
source and a `verified` flag. `eval_private/gold/` is an unverified starting draft.

## Build order
1. `extract_truth.py`: raw-file fields from one mzML per study (pyOpenMS). Use a conversion
   with all MS levels: the pipeline's MTBLS2469 mzML files are MS1-only (`-L=1`), so they
   cannot answer `ms2_acquired`; re-convert one file per MetaboLights study without `-L`.
2. `fetch_text.py`: methods section and full text from Europe PMC; repository text from the
   MetaboLights / Workbench / MassIVE APIs.
3. `run_extraction.py`: prompt, call the Claude API, save raw responses (cached, re-runnable).
4. `score.py`: metrics above, per condition, field and model.
Pilot on 5 studies end to end before scaling.

## Open decisions
- Models: one strong and one small model (e.g. Opus and Haiku) is enough to show a pattern.
- Whether to repeat each call (temperature 0 vs several samples) to measure consistency.
