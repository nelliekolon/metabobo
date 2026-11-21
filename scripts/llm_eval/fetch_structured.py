#!/usr/bin/env python3
"""Eval step 6a: draft the answer key for the text-only fields from the repositories' structured
submitter fields (see SPEC.md): chromatography, column, sample_type, group_definitions.

usage: python3 fetch_structured.py <targets.tsv> <inputs_dir> <out.tsv>

MetaboLights: the acquisition's assay table (a_*.txt; chromatography and column columns, on the
row of the target's raw file) and the study's sample table (s_*.txt; organism part, and the
factor-value columns, which define the groups). The models never saw these tables: their
repository text is the study description and protocols. Workbench: the mwTab's CH: and
SUBJECT_SAMPLE_FACTORS lines (part of the repository text the models saw).
This is a draft: every value is checked by hand against the paper before it is used to score.
"""
import csv, ftplib, re, sys, time, urllib.request
from collections import Counter

UA = {"User-Agent": "metapool-eval/0.1 (academic research script)"}
EBI = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public"


def get(url):
    for t in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            err = e; time.sleep(3 * (t + 1))
    raise err


def table(text):
    rows = list(csv.reader(text.splitlines(), delimiter="\t"))
    return [h.strip() for h in rows[0]], rows[1:]


def cols(head, pat):
    return [i for i, h in enumerate(head) if re.search(pat, h, re.I)]


def mtbls(ftp, t):
    study = t["study"]
    head, rows = table(get(f"{EBI}/{study}/{t['acquisition']}"))
    base = t["raw_file"].split("/")[-1]
    row = next((r for r in rows if any(base in c for c in r)), rows[0] if rows else [])
    pick = lambda pat: "; ".join(dict.fromkeys(row[i].strip() for i in cols(head, pat)
                                               if i < len(row) and row[i].strip()))
    out = dict(chromatography=pick(r"column type|chromatography type"),
               column=pick(r"column model"), chrom_instrument=pick(r"chromatography instrument"))
    names = [n.split("/")[-1] for n in ftp.nlst(f"/pub/databases/metabolights/studies/public/{study}")]
    s_files = [n for n in names if n.startswith("s_") and n.endswith(".txt")]
    if s_files:
        sh, srows = table(get(f"{EBI}/{study}/{s_files[0]}"))
        part = Counter(r[i].strip() for r in srows for i in cols(sh, r"organism part") if i < len(r) and r[i].strip())
        fac = cols(sh, r"^factor value")
        groups = Counter(" | ".join(f"{sh[i][13:-1]}={r[i].strip()}" for i in fac if i < len(r)) for r in srows)
        out["sample_type"] = "; ".join(f"{k} ({v})" for k, v in part.most_common())
        # the organism part of the sample behind THIS raw file (studies can mix urine, plasma, serum)
        sn = next((row[i].strip() for i in cols(head, r"^sample name$") if i < len(row)), "")
        sni, opi = cols(sh, r"^sample name$"), cols(sh, r"organism part")
        hit = next((r for r in srows if sni and sni[0] < len(r) and r[sni[0]].strip() == sn), None)
        out["file_sample"] = sn
        out["file_sample_type"] = hit[opi[0]].strip() if hit and opi and opi[0] < len(hit) else ""
        out["groups"] = "; ".join(f"{k} ({v})" for k, v in groups.most_common(12))
        out["n_samples"] = len(srows)
    return out


def workbench(t, inputs):
    txt = open(f"{inputs}/{t['study']}/repository.txt").read()
    line = lambda key: "; ".join(dict.fromkeys(m.strip() for m in re.findall(rf"^{key}\s+(.+)$", txt, re.M)))
    fac = re.findall(r"^SUBJECT_SAMPLE_FACTORS\s+\S*\s+\S+\s+(.+?)(?:\t|$)", txt, re.M)
    return dict(chromatography=line("CH:CHROMATOGRAPHY_TYPE"), column=line("CH:COLUMN_NAME"),
                chrom_instrument=line("CH:INSTRUMENT_NAME"),
                sample_type=line("SU:SAMPLE_TYPE") or line("CO:SAMPLE_TYPE"),
                groups="; ".join(f"{k} ({v})" for k, v in Counter(f.strip() for f in fac).most_common(12)),
                n_samples=len(fac))


def main():
    targets, inputs, out = sys.argv[1:4]
    ftp = ftplib.FTP("ftp.ebi.ac.uk", timeout=120); ftp.login()
    rows = []
    for t in csv.DictReader(open(targets), delimiter="\t"):
        try:
            d = mtbls(ftp, t) if t["source"] == "MetaboLights" else workbench(t, inputs)
        except Exception as e:
            d = dict(error=str(e)[:200])
        rows.append(dict(target=t["target"], study=t["study"], **d))
        print(t["target"], {k: str(v)[:70] for k, v in d.items()}, flush=True)
    ftp.quit()
    keys = ["target", "study", "chromatography", "column", "chrom_instrument", "sample_type", "file_sample",
            "file_sample_type", "groups",
            "n_samples", "error"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, keys, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
