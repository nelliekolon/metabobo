# metabobo

**Can a model trained on one public metabolomics study work on another?**

Labs publish raw data from hundreds of human blood metabolomics studies, but almost every model built on them is trained and tested inside a single study. Pooling studies would give far bigger cohorts than any one lab can collect. I wanted to find out what it actually takes.

I ran four public studies (127 people, about 40 GB of raw LC-MS data) through one pipeline with the same settings, matched features across them, and tested whether a simple model transfers from one study to the next. Then I tested whether an LLM can fill in the method details needed to scale this up, checking its answers against the raw files.

The full write-up is in [WRITEUP.md](WRITEUP.md).

## What I found

**The shortcuts were there every time.** Most of the work went into catching things a model could learn instead of biology, each found by a check on the data rather than by a final number looking odd:

- In one study, run day explained the main signal as well as COVID severity did.
- In the training study, sex and diagnosis were tangled (14 of 19 women had MCI).
- One of three "replicate" injections per sample was a different experiment (a fragmentation run).
- Saliva and stool pools sat inside a blood study and dominated its PCA.

**Few features survive the change of study.** Only 88 of the reference study's 7,602 well-detected features are found in all four studies.

**Transfer is promising but underpowered.** A sex classifier with AUC 0.72 in its own study scores 0.68 on another study's plasma, 0.81 after per-study standardization, and 0.58 on serum from the same people. With 20 test people, the intervals overlap. Age never beats predicting the mean.

**LLMs extract method details faithfully; the papers are the problem.** Across 45 studies, I compared what Claude (Opus 5.5 and Haiku 5.5) extracted from each paper with values read directly from the raw files. On the fields a raw file could check, Opus never misread or misattributed a value, and quotes that weren't in the text were rare (9 of 981). But about 1 in 4 studies' own text contradicts its raw data, so a faithful extractor inherits those errors. Checking against the raw data is what catches them.

## Repo layout

| Path | What's there |
| --- | --- |
| `scripts/` | The pipeline, steps 0-5 (see [docs/PIPELINE.md](docs/PIPELINE.md)) |
| `scripts/llm_eval/` | The LLM metadata-extraction eval: study selection, text fetching, extraction, scoring ([SPEC.md](scripts/llm_eval/SPEC.md)) |
| `processed/<study>/` | Sample sheets and QC reports for each study |
| `processed/crossstudy/` | Feature matching and transfer results (`match_report.txt`, `ml_report.txt`) |
| `processed/paper_survey.tsv` | 136 candidate public studies and their papers |

Large files (raw data, feature tables) are not in the repo. `ftp_download.py` re-downloads the raw data, and the scripts regenerate everything else.

## Running it

Python 3 with the versions in `requirements.txt` (pyOpenMS, scikit-learn, pandas < 3). Step-by-step commands and study-specific notes are in [docs/PIPELINE.md](docs/PIPELINE.md).

## Limitations

This is a pilot, not a finished resource. The main limit is size: 127 people in all and 20 in the key test set, so every cross-study number has a wide interval. The four studies share no disease, so only sex and age could be tested. Everything is positive ion mode with one fixed noise cutoff across instruments.

## Data

MassIVE: MSV000084556, MSV000085256, MSV000094097. MetaboLights: MTBLS2469. Original papers are listed in `processed/paper_survey.tsv`.

Built with Claude Code.
