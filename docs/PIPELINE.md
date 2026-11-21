# Pipeline details

Working notes for re-running the pipeline. For what the project found, see the [README](../README.md) and [WRITEUP.md](../WRITEUP.md).

## Folders
- pilot_data/<study>/        raw downloads (mzML, or Thermo .raw for MetaboLights)
- processed/<study>/features_v1/        one .featureXML per run (step 1)
- processed/<study>/table_v2_noalign/   CURRENT feature table for the study (step 2-3)
- scripts/

## Environment
Run everything with .venv/bin/python (versions in requirements.txt). Keep pandas < 3:
pandas 3 breaks check_table.py.

## Pipeline (fixed parameters, same for every study)
0. scripts/make_samples.py pilot_sample_metadata.tsv processed
   Writes processed/<study>/samples.tsv: person, visit, matrix, pool vs subject, clinical
   group, include flag (blood only; saliva/stool pools and MTBLS2469 _test runs excluded).
   Per-study filename rules live in RULES at the top of the script.
0b. (MetaboLights .raw only) scripts/convert_raw.py pilot_data/<study> pilot_data/<study>_mzML --samples processed/<study>/samples.tsv
   ThermoRawFileParser (tools/, v2.0.0-dev osx-arm64): MS1 only, centroided, plain mzML.
1. scripts/detect_features.py <mzML dir> <out dir>
   OpenMS FeatureFinderMetabo on MS1: 5 ppm, noise 3000, SNR 3, single-trace features kept.
2. scripts/build_table.py <feature dir> <out dir> --no-align --samples processed/<study>/samples.tsv
   QT feature linking (10 ppm, 15 s). RT alignment is OFF by default for now: on
   MSV000084556 the runs drift <2 s and PoseClustering alignment broke 8 runs.
2b. (studies with pools spread over several batches) scripts/batch_correct.py <table dir> <out dir>
   Per feature, shifts each batch (samples.tsv "batch") so its pools match the overall pool median
   (log2). Batches without pools borrow the nearest batch's offset. Reports held-out pool CV.
3. scripts/check_table.py <table dir> 0.8
   Keeps features found in >=80% of samples, log2, min-imputation, median normalization;
   reports pooled-QC CV, within-person similarity, PCA (by visit/matrix/role) and
   known-compound hits. Reads ../samples.tsv; writes check_report.txt via tee.
4. scripts/match_studies.py processed/crossstudy MSV000085256=<table dir> <study>=<table dir> ...
   Matches each study's filtered features to a reference study: m/z within 10 ppm plus an
   order-preserving (isotonic) RT map fitted on unambiguous high-intensity m/z anchors.
   RT tolerance from cross-validated anchor error. Known compounds are held out as a check.
   Output: feature_map.tsv (reference feature -> matched feature per study), match_report.txt.
5. scripts/ml_crossstudy.py processed/crossstudy MSV000085256=<table dir> <study>=<table dir> ...
   Sex (logistic) and age (ridge), trained on MSV000085256, tested on each other study on
   shared features; no harmonization vs per-study z-scoring; per-person scoring; bootstrap CIs.

## Study notes
- MSV000085256: one plasma sample per person (MCI vs Normal); scan range m/z 100+.
- MSV000094097: plasma AND serum from the same people; PC1 is plasma vs serum.
  Scan range starts at m/z 150, so creatinine/proline/betaine/leucine are out of range.
- MTBLS2469: each sample injected 3x; _P_2 is the MS2 (DDA) run, _P_3/_P_4 are MS1-only. _P_2
  has fewer MS1 scans and agrees less in all 50 samples (scripts/check_injections.py), so only
  _P_3 is used. Signal drifts between run days (pool CV
  64%); batch_correct.py brings held-out pool CV to 47%. CURRENT table: table_v3_batchcorr.
  Only 8 people over 9 run days, and COVID severity is confounded with run day: not usable
  as a training set for severity.
- RT differs between studies (e.g. tryptophan 56 s vs 140 s): match across studies by
  m/z plus an RT mapping built from known compounds, not raw RT.

## Downloading
python3 ftp_download.py pilot_download_usis.tsv pilot_data [--only STUDY]
