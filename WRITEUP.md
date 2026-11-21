# What breaks when you combine public metabolomics data?

Labs around the world measure small molecules in people's blood, called metabolites, in response to disease, diet, exercise, etc, and many publish their raw data online. If you could combine those studies, you'd have many more data points than one study alone has. I wanted to see if that could work: if you train a model on one study, does it still work on another study's data?

I ran four public human blood metabolomics studies (127 people) through one pipeline with the same settings, to find out what it takes for a model trained on one study to work on another.

Most of what I found is what goes wrong first. In one study, run day explained as much of the signal as disease severity did. In the training cohort, sex and diagnosis were tangled. One injection per sample turned out to be a fragmentation run, and saliva and stool data sat inside a blood study. Only 88 of the reference study's 7,602 well-detected features appear in all four studies. A small transfer test is promising but underpowered: a sex classifier that reaches AUC 0.72 in its own study scores 0.68 on another study's plasma, 0.81 after per-study standardization and 0.58 on serum from the same people, with 20 people and overlapping intervals.

The metadata needed to harmonize more studies is often wrong too. In a follow-up on 45 studies, an LLM extracted method details almost exactly as written, but about 1 in 4 studies' text contradicts its own raw files, so checks against the raw data, not better extraction, are what catch those errors.

## The problem

Public repositories hold hundreds of human blood metabolomics studies, but almost every model built on them is trained and tested inside one study. Studies differ in instrument, chromatography, scan range, sample preparation and lab, so the same molecule shows up at a different time and intensity in each. Pooling studies would give far larger cohorts than any one lab can collect, yet nobody knows how much of a model's accuracy is biology and how much is the study it came from.

This is a distribution-shift problem, the same one ML faces whenever a model meets data from a new source. I asked a narrow version of it: **does a model trained on one public study still work on another, and does harmonization help?**

## What I built

One pipeline takes each study from raw public files to a QC'd feature table, with the same parameters for every study.

It downloads from the repositories' own servers (resumable, about 40 GB in all) and converts Thermo .raw files to mzML for the MetaboLights study. Each study gets one sample sheet recording person, visit, sample type, pool, batch and whether a run is included, and the rules that read those facts out of file names live in one place. Features are detected with OpenMS at 5 ppm and a fixed noise cutoff, linked across runs, kept if found in 80% of them, and summarized in a QC report. Where pooled QC samples span several run days, a pooled-QC correction removes each day's offset.

| Study | Repository | People | Runs used | Sample type | Clinical label |
| --- | --- | --- | --- | --- | --- |
| MSV000084556 | MassIVE | 20 | 60 | plasma | rheumatoid arthritis, 3 visits |
| MSV000085256 | MassIVE | 78 | 78 | plasma | MCI vs normal |
| MSV000094097 | MassIVE | 21 | 83 + 17 pools | plasma and serum | osteoarthritis, 2 visits |
| MTBLS2469 | MetaboLights | 8 | 22 + 28 pools | plasma | COVID severity |

The four studies share no disease, so the only labels they all carry are sex and age. MSV000085256, the largest, is the reference that every other study is matched to and the one models are trained on.

## What broke, and how I caught it

Most of the work went here. Each of these problems would have silently corrupted a result, and each was caught by a check on the data rather than by a final number looking odd.

**Retention-time alignment made things worse.** Runs in MSV000084556 drift by under 2 seconds, yet the aligned feature tables came out worse than the unaligned ones: alignment distorted 8 runs. It is now off by default.

**The wrong sample types.** MSV000094097 includes saliva and stool pools. They dominated both the PCA (the first two components were sample type, not people) and the 80% detection filter, so the sample sheets are now blood only. Even within blood, plasma against serum outweighs the differences between people: the first component split plasma at -16 and serum at +16, so I score the two separately.

**One injection was a different experiment.** Each MTBLS2469 sample was injected three times. Injection 2 had 25% fewer MS1 scans and agreed less with its neighbours (r 0.68 against 0.82, in all 50 samples). Converting with every MS level showed why: it held 3,424 MS2 scans, so it was the fragmentation run. I use injection 3 only.

**Run day was confounded with disease severity.** In MTBLS2469 the pooled QC samples varied by a median CV of 64% across run days, and run day explained 99% of the first principal component while severity explained 94%: a model could learn either one and look the same. A pooled-QC correction shifts every run by its day's offset.

![Held-out pooled QC samples in MTBLS2469 by run day, before and after batch correction](docs/figures/batch_drift.png)

*13 held-out pools; offsets fitted on the other 15. Data: `processed/figures/batch_drift_MTBLS2469.tsv`.*

The correction is fitted on half of each day's pools and judged on the other half, so it is never graded on the pools it was fitted to. On those held-out pools the CV falls from 64% to 47%. Scatter within a day remains: the correction removes each day's offset, not run-to-run noise.

**A fixed noise cutoff means different things on different instruments.** MTBLS2469 yields about 50,000 features per run, against 10,000 to 20,000 elsewhere. I kept the cutoff fixed, so every study is treated alike, and count it as a limitation.

**Sex and diagnosis are tangled in the training cohort.** In MSV000085256, 14 of the 19 women have MCI, and diagnosis alone predicts sex at AUC 0.72. A sex classifier could score well by learning diagnosis. Within each diagnosis group the sex signal holds (0.69 within MCI, 0.68 within normal), so it is not only that.

**The metadata is sometimes wrong.** One repository labels healthy controls "Alzheimer's disease". This is what led to the follow-up below.

The pattern across all of these: the shortcut a model could have learned (run day, sample type, diagnosis) was there every time.

## Matching features across studies

A model trained on one study can only use features that can be found again in another. Only 88 of the reference study's 7,602 well-detected features are found in all four studies; 564 are found in at least two.

All four studies use reverse-phase chromatography, so compounds elute in the same order but on different time scales, from about 5 to 20 minutes. I map each study's retention times onto the reference with an order-preserving fit through unambiguous, high-intensity m/z matches, and set the match tolerance from the fit's cross-validated error. Fourteen known compounds are held out of the fit as an independent check.

![Retention-time maps from each study onto the reference](docs/figures/rt_map.png)

*Anchors used in the fit only. Data: `processed/crossstudy/rt_map_*.tsv`.*

| Study vs reference | Features matched | Known compounds matched |
| --- | --- | --- |
| MSV000084556 | 1,309 (50%) | 5 of 5 |
| MSV000094097 | 541 (19%) | 4 of 6 |
| MTBLS2469 | 1,058 (17%) | 5 of 10 (4 on the wrong peak, 1 missed) |

MSV000084556 and the reference likely share a lab and a method, which is why they overlap so much more. Different labs on similar instruments share a small fraction of their features. MTBLS2469 runs a longer gradient, and its first 7 minutes map onto the reference's first 2, where its anchors scatter most and where its known compounds land on the wrong peak.

## A pilot: does a model survive the change of study?

Sex and age are the only labels all four studies share. Models are trained on MSV000085256 (76 people with known sex) using only the features matched in each test study, and scored once per person. MSV000084556 (2 men) and MTBLS2469 (all male) cannot test sex, which leaves MSV000094097, where the same people gave both plasma and serum (20 in the test set).

| Sex classifier, trained on MSV000085256 | AUC [95% CI] |
| --- | --- |
| Within its own study, all features | 0.72 [0.56-0.86] |
| Within its own study, the 541 features shared with MSV000094097 | 0.71 [0.56-0.85] |
| MSV000094097 plasma, no harmonization | 0.68 [0.42-0.91] |
| MSV000094097 plasma, per-study z-scoring | 0.81 [0.60-0.98] |
| MSV000094097 serum, no harmonization | 0.68 [0.43-0.93] |
| MSV000094097 serum, per-study z-scoring | 0.58 [0.31-0.86] |

*Per-person scoring, 2,000 bootstrap resamples. Full output: `processed/crossstudy/ml_report.txt`.*

A sex classifier transfers to a new study's plasma, and simple per-study standardization helps, but a switch to serum in the same people breaks it. Within its own study the model reaches AUC 0.72 on all features and 0.71 on the 541 features it shares with MSV000094097, so restricting to shared features costs little. On MSV000094097 plasma it scores 0.68 untouched and 0.81 after z-scoring each study separately; on serum it falls to 0.58.

Two cautions. The plasma gain is promising, not proven: with 20 people the intervals overlap. And the confound from the training cohort is checked rather than assumed away: the model still separates sex within MCI (0.69) and within normal (0.68). Age is a null result: no model beats guessing the training mean, even within one study.

## Follow-up: can an LLM fill in the missing metadata?

Harmonizing more studies needs each one's method details (ion mode, scan range, column, sample type, groups), and those live in free text that is often vague or wrong. So I asked whether an LLM can extract them, and checked its answers against something other than the text: the raw files themselves.

The short answer is yes, as faithfully as the text allows. On the fields a raw file can check, Opus 5.5 never misattributed or misread a value. But about 1 in 4 studies' own text contradicts its raw data, so a faithful extractor inherits those errors, and checking against raw files is what catches them.

**Setup.** 45 public LC-MS blood studies with an open-access paper: 39 MetaboLights studies where I measured the truth from one raw file (ion mode, scan window, run length, instrument), and 6 Metabolomics Workbench studies scored against the repository's structured fields. Two models (Opus 5.5 and Haiku 5.5) read three texts per study (the methods section, the full paper, the repository description) and returned 10 fields as structured JSON, each with a verbatim quote. Every quote is checked against the input text automatically, and every wrong answer (157) was labelled by hand by its cause. The design is in [scripts/llm_eval/SPEC.md](scripts/llm_eval/SPEC.md).

**Most wrong answers are the text's fault.** 10 of 39 studies state at least one value their own raw file contradicts: papers naming a SCIEX QTRAP, an Agilent Q-TOF or a Waters QTOF for data recorded on a Thermo Q Exactive, or a scan range of 50-1500 for a 70-1050 window. In one study the repository is wrong where the paper is right. Metabolon's stock sentence "the scan range varied slightly between methods but covered 70-1000 m/z" is wrong for its negative-mode method, which records 80-1000, and both models copy it in every such study.

**The models' own errors are rare and of one kind.** On the raw-file fields, Haiku's 7 errors are all attribution: reporting "both" ion modes for a single-mode file, or the positive-mode gradient for a negative-mode file. Opus made none. Where a text describes two methods that could both fit the file, both models pick one instead of answering "not stated" (7 such errors for Opus and 8 for Haiku on the text-only fields). Quotes that are not in the text are rare too: 9 of 981 for Opus and 12 of 1,138 for Haiku, mostly ion mode read from "NEG" in a file name.

**Fields the raw file cannot check.** For these I built an answer key from the repositories' structured submitter fields (assay and sample tables the models never saw), corrected by hand. Accuracy is 98% and 94% (Opus, Haiku) for chromatography, 94% and 91% for column chemistry, and 93% and 91% for sample type. Here too the key and the paper disagree: in 5 studies the paper says plasma while the sample table says serum. For the groups each study compares, graded against references I wrote from the sample groups and the paper, 133 of 135 Opus answers and 127 of 135 Haiku answers are complete. The one wrong answer is Opus guessing at groups the text leaves undefined, which I count as worse than "not stated". A second rater, blind to my grades, agreed on 17 of 20 answers.

Accuracy is similar whichever text the model reads (75-84%), and both models are far above always giving the most common value (13-55% per field). The whole run cost about $13, and Haiku costs about 1/40 as much as Opus for nearly the same accuracy.

## Limitations

The main limitation is size: 127 people in all and 20 in the key test set, so every cross-study number has a wide interval. The four studies share no disease, so only sex and age could be tested, and sex in one test study only. Everything is positive ion mode with one fixed noise cutoff across instruments, and the MTBLS2469 correction borrows offsets for 4 samples run on days without pools.

For the LLM evaluation, only wrong answers were hand-checked, so a right answer could be right by luck. The Workbench prompts named no raw file, so their ion-mode answers are excluded. And the group references are mine, with the blind second rating covering only 20 answers from 15 studies.

## Next steps

- Scale to more of the 136 candidate studies; MTBLS2542 is next (320 runs, 51 GB).
- Harmonize better than z-scoring, for example by matching on reference compounds or pooled QCs.
- Use the LLM extractor to fill in method metadata across the candidate studies, with raw-file checks flagging where the text is wrong.

## Code and data

- Code, sample sheets and QC reports: this repo
- Data: MassIVE MSV000084556, MSV000085256, MSV000094097; MetaboLights MTBLS2469
- Original papers for each study: DOIs in `processed/paper_survey.tsv`
