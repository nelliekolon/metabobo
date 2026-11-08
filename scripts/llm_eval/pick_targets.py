#!/usr/bin/env python3
"""Pick one acquisition (and, for MetaboLights, one raw file) per study for the full eval set.

usage: python3 pick_targets.py <eval_candidates.tsv> <out_targets.tsv> [--exclude ID ...]
Studies with paper_status 'review' and any --exclude IDs are left out.
Per study, among its acquisitions, the choice balances the set as it goes: HILIC first when the
study has it and the set has fewer HILIC than reverse-phase runs, then the ion mode the set has
fewer of. MetaboLights: the acquisition's assay table (a_*.txt) gives its raw files; the first
raw file listed is used and its size recorded. Workbench: no single-file raw download, so the
answer key comes from the repository's structured fields (truth_tier = repository).
"""
import csv, json, re, sys, time, urllib.request, ftplib

UA = {"User-Agent": "metabobo-eval/0.1 (academic research script)"}
EBI = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public"


def get(url):
    for t in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            err = e; time.sleep(3 * (t + 1))
    raise err


SKIP = re.compile(r"blank|solvent|wash|urine|feces|fecal|stool|tissue|cell|std|standard|test|condition", re.I)


def ftp_assays(ftp, study):
    """assay tables as they are actually named on the server (the API's names can differ)."""
    names = [n.split("/")[-1] for n in ftp.nlst(f"/pub/databases/metabolights/studies/public/{study}")]
    out = []
    for fn in names:
        if not (fn.startswith("a_") and fn.endswith(".txt")):
            continue
        low = fn.lower()
        ion = "positive" if re.search(r"pos", low) else "negative" if re.search(r"neg", low) else ""
        chrom = "HILIC" if "hilic" in low else "reverse phase" if re.search(r"reverse|rp|c18", low) else ""
        out.append([fn, ion, chrom, "lipid" if "lipid" in low else ""])
    return out


def raw_files(study, assay_file):
    """raw files of an assay (or its derived mzML/mzXML when no raw file is listed); blank,
    solvent, wash and non-blood runs moved to the end."""
    rows = list(csv.reader(get(f"{EBI}/{study}/{assay_file}").splitlines(), delimiter="\t"))
    head = [h.strip().lower() for h in rows[0]]
    raw = [i for i, h in enumerate(head) if h == "raw spectral data file"]
    der = [i for i, h in enumerate(head) if h == "derived spectral data file"]
    def pick(cols, exts):
        return [r[i].strip() for r in rows[1:] for i in cols
                if i < len(r) and r[i].strip().lower().endswith(exts)]
    files = pick(raw, (".raw", ".mzml", ".mzxml")) or pick(der, (".mzml", ".mzxml"))
    return sorted(dict.fromkeys(files), key=lambda f: (bool(SKIP.search(f)), "qc" not in f.lower()))


def main():
    cand_f, out_f = sys.argv[1], sys.argv[2]
    excl = set(sys.argv[sys.argv.index("--exclude") + 1:]) if "--exclude" in sys.argv else set()
    cands = [r for r in csv.DictReader(open(cand_f), delimiter="\t")
             if r["paper_status"] != "review" and r["study"] not in excl]
    n = {"positive": 0, "negative": 0, "HILIC": 0, "RP": 0}
    ftp = ftplib.FTP("ftp.ebi.ac.uk", timeout=120); ftp.login(); ftp.voidcmd("TYPE I")
    rows = []
    for c in cands:
        if c["source"] == "MetaboLights":
            acq = ftp_assays(ftp, c["study"])
            acq = [a for a in acq if a[3] != "lipid"] or acq          # metabolomics before lipidomics
        else:
            acq = json.loads(c["acquisitions"]) if c["acquisitions"].endswith("]") else []
        acq = [a for a in acq if a[1] in ("positive", "negative")] or acq
        if not acq:
            print(f"{c['study']}: no acquisitions listed, skipped"); continue
        def key(a):
            hil = "hilic" in (a[2] or "").lower()
            want_hilic = n["HILIC"] < n["RP"] / 2
            non_blood = bool(re.search(r"urine|feces|fecal|stool|tissue|cell", a[0], re.I))
            return (non_blood, 0 if hil == want_hilic else 1, n.get(a[1], 0))
        a = sorted(acq, key=key)[0]
        files = []
        if c["source"] == "MetaboLights":
            for cand in sorted(acq, key=key):            # first acquisition that lists usable files
                try:
                    files = raw_files(c["study"], cand[0])
                except Exception as e:
                    print(f"{c['study']}: {cand[0]} unreadable ({e})")
                    continue
                if files and not SKIP.search(files[0]):     # an acquisition with a usable blood run
                    a = cand; break
                files = files if files else []
        row = dict(study=c["study"], source=c["source"], target=f"{c['study']}_{(a[1] or 'x')[:3].upper()}",
                   acquisition=a[0], ion_mode_repo=a[1], chromatography_repo=a[2],
                   instrument_repo=a[3] if c["source"] != "MetaboLights" else c["instrument"],
                   doi=c["doi"], pmcid=c["pmcid"], raw_file="", raw_bytes="",
                   truth_tier="raw file" if c["source"] == "MetaboLights" else "repository")
        if c["source"] == "MetaboLights":
            if not files:
                print(f"{c['study']}: no raw files in {a[0]}, skipped"); continue
            row["raw_file"] = files[0]
            try:
                row["raw_bytes"] = ftp.size(f"/pub/databases/metabolights/studies/public/{c['study']}/{files[0]}")
            except ftplib.all_errors:
                try:
                    row["raw_bytes"] = ftp.size(f"/pub/databases/metabolights/studies/public/{c['study']}/FILES/{files[0]}")
                    row["raw_file"] = "FILES/" + files[0]
                except ftplib.all_errors:
                    row["raw_bytes"] = ""
        n[a[1]] = n.get(a[1], 0) + 1
        n["HILIC" if "hilic" in (a[2] or "").lower() else "RP"] += 1
        rows.append(row)
        print(f"{row['target']:18s} {row['chromatography_repo'][:16]:16s} {row['raw_file'][-45:]:45s} "
              f"{(int(row['raw_bytes'])/1e6 if row['raw_bytes'] else 0):7.0f} MB", flush=True)
    ftp.quit()
    with open(out_f, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    print(f"\n{len(rows)} targets; balance {n}")


if __name__ == "__main__":
    main()
