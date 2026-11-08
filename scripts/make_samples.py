#!/usr/bin/env python3
"""Step 0: one sample sheet per study, so later steps never parse filenames themselves.

usage: python3 make_samples.py <metadata.tsv> <processed_dir>
Writes <processed_dir>/<study>/samples.tsv with columns:
  sample     run name (file basename without extension; matches feature table columns)
  subject    person ID ("pool" for pooled QC samples)
  timepoint  visit within a person, if the study has repeat visits
  matrix     blood plasma / blood serum / saliva / stool
  role       subject or pool
  group      study-specific clinical group (e.g. MCI/Normal, Mild/Severe)
  health, disease, sex, age   from the Pan-ReDU metadata
  injection  repeat-injection number, where a study injects each sample more than once
  batch      acquisition batch as an ISO date (run day), where the filenames record it
  include    1 if used in the blood feature table, 0 otherwise

Subject/timepoint rules live in RULES below, one per study; everything else comes from metadata.
"""
import csv, os, re, sys

MISSING = {"", "missing value", "not applicable"}


def clean(v):
    return "" if v is None or v.strip() in MISSING else v.strip()


def rule_msv084556(path, name):
    # Raw_data/PD001_T1.mzML -> person PD001, visit T1
    m = re.match(r"(PD\d+)_(T\d+)$", name)
    return dict(subject=m[1], timepoint=m[2])


def rule_msv085256(path, name):
    # MCI/324_Plasma[_<acquisition timestamp>].mzML -> one sample per person
    m = re.match(r"(\d+)_Plasma(_\d{14})?$", name)
    return dict(subject=m[1], timepoint="", group=path.split("/")[-2])


def rule_msv094097(path, name):
    # Blood/Plasma/VAOAD_015_Day15_S.mzML, OAD_005_Day15_P_2.mzML (reinjection), Pool_blood_serum_3
    m = re.match(r"((?:VA)?OAD)_(\d+)_(Basal|Day\d+)_[PS](?:_(\d+))?$", name)
    if m:
        return dict(subject=f"{m[1]}_{m[2]}", timepoint=m[3], group=m[1])
    low = name.lower()
    matrix = ("blood plasma" if "plasma" in low else "blood serum" if "serum" in low
              else "saliva" if "saliva" in low else "stool" if "stool" in low else "")
    return dict(subject="pool", timepoint="", matrix=matrix)


def rule_mtbls2469(path, name):
    # Patient samples/Reg_x/04102020_Covid_TP_Mild_9030_1_P325_P_2 -> person 9030, visit 1, group Mild
    # QCPool/17102020_QCPool_Set4_P_3[_test]
    # Every sample was injected 3 times (_P_2/_3/_4). Injection 3 is the one used: on a test sample
    # _3 and _4 agreed (r=0.83) but _2 had fewer spectra, more features and a shifted LPC RT (r=0.71).
    # Filenames start with the run day (DDMMYYYY), used as the batch for batch_correct.py.
    day = f"{name[4:8]}-{name[2:4]}-{name[0:2]}"
    m = re.match(r"\d{8}_Covid(?:Pos|_TP_(.+?))_(\d{4})_(\d+)_P\d+_P_(\d+)$", name)
    if m:
        return dict(subject=m[2], timepoint=m[3], group=m[1] or "CovidPos", injection=m[4],
                    use=m[4] == "3", batch=day)
    m = re.match(r"\d{8}_QCPool_Set\d+_P_(\d+)(_test)?$", name)
    if m:
        return dict(subject="pool", timepoint="", injection=m[1], use=m[1] == "3" and not m[2],
                    batch=day)
    raise ValueError(name)


RULES = {"MSV000084556": rule_msv084556, "MSV000085256": rule_msv085256,
         "MSV000094097": rule_msv094097, "MTBLS2469": rule_mtbls2469}
COLS = ["sample", "subject", "timepoint", "matrix", "role", "group", "health", "disease",
        "sex", "age", "injection", "batch", "include"]


def main():
    meta_file, out_root = sys.argv[1], sys.argv[2]
    by_study = {}
    with open(meta_file) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by_study.setdefault(r["ATTRIBUTE_DatasetAccession"], []).append(r)

    for ds, rows in by_study.items():
        if ds not in RULES:
            print(f"{ds}: no naming rule yet, skipped")
            continue
        out = []
        for r in rows:
            path = r["filename"].split("/", 1)[1]               # drop leading "f."
            name = os.path.splitext(os.path.basename(path))[0]
            p = RULES[ds](path, name)
            role = "pool" if p["subject"] == "pool" else "subject"
            matrix = p.get("matrix") or clean(r["UBERONBodyPartName"])
            if not matrix and ds.startswith("MTBLS"):
                matrix = "blood plasma"
            out.append(dict(
                sample=name, subject=p["subject"], timepoint=p.get("timepoint", ""),
                matrix=matrix, role=role, group=p.get("group", ""),
                health=clean(r["HealthStatus"]), disease=clean(r["DOIDCommonName"]),
                sex=clean(r["BiologicalSex"]), age=clean(r["AgeInYears"]),
                injection=p.get("injection", ""), batch=p.get("batch", ""),
                # blood only, and one injection per sample where a study has several
                include=int(matrix.startswith("blood") and p.get("use", True))))
        d = os.path.join(out_root, ds)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "samples.tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, COLS, delimiter="\t")
            w.writeheader()
            w.writerows(sorted(out, key=lambda x: x["sample"]))
        inc = [x for x in out if x["include"]]
        n_subj = len({x["subject"] for x in inc if x["role"] == "subject"})
        n_pool = sum(x["role"] == "pool" for x in inc)
        print(f"{ds}: {len(out)} runs, {len(inc)} included "
              f"({len(inc) - n_pool} from {n_subj} people, {n_pool} pools)")


if __name__ == "__main__":
    main()
