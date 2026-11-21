#!/usr/bin/env python3
"""Choose the full eval set: for each candidate study, find its ORIGINAL open-access paper and
list its acquisitions (ion mode, chromatography, instrument) so the set can be chosen for variety.

usage: python3 select_studies.py <paper_survey.tsv> <out.tsv>
Candidates: open-access paper found by Europe PMC + >= 3 structured method fields in the repository.

Original-paper rule (a paper that cites the accession is often a later reuse, see pilot_papers.tsv):
  1. the repository's own linked paper, if it is open access, else
  2. among Europe PMC hits: published within 2 years of the dataset, sharing >= 2 content words
     with the dataset title; best overlap wins.
  Anything else is 'review' (checked by hand), never silently accepted.
Acquisitions: Workbench analyses (chromatography_type, ion_mode, instrument); MetaboLights assay
files (ion mode and column from the assay name) and raw file names (to pick one file per acquisition).
"""
import csv, json, re, sys, time, urllib.parse, urllib.request

UA = {"User-Agent": "metabopool-eval/0.1 (academic research script)"}
STOP = set("the a an of and in for with to from by on at is are as using study analysis human plasma serum "
           "blood metabolomics metabolomic metabolite metabolites lc ms lcms untargeted patients".split())


def get(url, as_json=True):
    for t in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
                b = r.read().decode("utf-8", "replace")
                return json.loads(b) if as_json else b
        except Exception as e:
            err = e; time.sleep(3 * (t + 1))
    return {"_error": str(err)}


def words(s):
    return {w for w in re.findall(r"[a-z0-9]{3,}", (s or "").lower()) if w not in STOP}


def epmc(query):
    u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
        {"query": query, "format": "json", "resultType": "core", "pageSize": 25})
    d = get(u)
    return d.get("resultList", {}).get("result", []) if isinstance(d, dict) else []


def repo(acc):
    """title, year, linked DOIs, acquisitions [(label, ion_mode, chromatography, instrument)]."""
    if acc.startswith("MTBLS"):
        c = get(f"https://www.ebi.ac.uk/metabolights/ws/studies/public/study/{acc}").get("content", {})
        dois = [p.get("doi", "").lower() for p in c.get("publications") or [] if p.get("doi")]
        acq = []
        for a in c.get("assays") or []:
            fn = (a.get("fileName") or "").lower()
            ion = "positive" if re.search(r"pos", fn) else "negative" if re.search(r"neg", fn) else ""
            chrom = "HILIC" if "hilic" in fn else "reverse phase" if re.search(r"rp|c18|reverse", fn) else ""
            acq.append((a.get("fileName", ""), ion, chrom, a.get("platform", "")))
        rel = c.get("studyPublicReleaseDate") or c.get("studySubmissionDate") or ""
        if isinstance(rel, (int, float)):        # epoch milliseconds
            rel = time.strftime("%Y", time.gmtime(rel / 1000))
        return c.get("title", ""), str(rel)[:4], dois, acq
    if acc.startswith("ST"):
        s = get(f"https://www.metabolomicsworkbench.org/rest/study/study_id/{acc}/summary")
        a = get(f"https://www.metabolomicsworkbench.org/rest/study/study_id/{acc}/analysis")
        a = {"1": a} if isinstance(a, dict) and "analysis_id" in a else (a if isinstance(a, dict) else {})
        acq = [(v.get("analysis_id", ""), (v.get("ion_mode") or "").lower(), v.get("chromatography_type", ""),
                v.get("ms_instrument_name", "")) for v in a.values() if isinstance(v, dict) and v.get("analysis_type") == "MS"]
        return s.get("study_title", ""), (s.get("submission_date") or "")[:4], [], acq
    return "", "", [], []


def original_paper(acc, title, year, repo_dois):
    hits = epmc(f'"{acc}"')
    for d in repo_dois:
        hits += epmc(f'DOI:"{d}"')
    seen, cands = set(), []
    for h in hits:
        key = (h.get("doi") or h.get("pmcid") or h.get("id", "")).lower()
        if key in seen:
            continue
        seen.add(key)
        oa = h.get("isOpenAccess") == "Y" and bool(h.get("pmcid"))
        ov = len(words(title) & words(h.get("title", "")))
        dy = abs(int(h.get("pubYear") or 0) - int(year)) if year and h.get("pubYear") else 99
        cands.append(dict(doi=(h.get("doi") or "").lower(), pmcid=h.get("pmcid", ""), oa=oa, overlap=ov,
                          dyear=dy, title=h.get("title", ""), year=h.get("pubYear", ""),
                          linked=(h.get("doi") or "").lower() in repo_dois))
    linked = [c for c in cands if c["linked"] and c["oa"]]
    if linked:
        return linked[0], "repository-linked", cands
    ok = sorted([c for c in cands if c["oa"] and c["dyear"] <= 2 and c["overlap"] >= 2], key=lambda c: -c["overlap"])
    if ok:
        return ok[0], "title+year match", cands
    return None, "review", cands


def main():
    survey, out = sys.argv[1], sys.argv[2]
    rows = []
    for r in csv.DictReader(open(survey), delimiter="\t"):
        if int(r["epmc_open_access"] or 0) == 0 or int(r["repo_method_fields"] or 0) < 3:
            continue
        acc = r["accession"]
        title, year, dois, acq = repo(acc)
        year = year or r["year"]
        paper, how, cands = original_paper(acc, title, year, dois)
        ions = sorted({a[1] for a in acq if a[1]})
        chroms = sorted({a[2] for a in acq if a[2]})
        rows.append(dict(study=acc, source=r["source"], year=year, fluid=r["fluid"], instrument=r["instrument"],
                         title=title[:150], paper_status=how,
                         doi=paper["doi"] if paper else "", pmcid=paper["pmcid"] if paper else "",
                         paper_title=(paper["title"][:150] if paper else ""), paper_year=paper["year"] if paper else "",
                         n_candidates=len(cands), n_acquisitions=len(acq),
                         ion_modes=";".join(ions), chromatography=";".join(chroms),
                         acquisitions=json.dumps(acq)[:2000]))
        print(f"{acc:12s} {how:18s} ions={rows[-1]['ion_modes']:18s} chrom={rows[-1]['chromatography']:25s} "
              f"{rows[-1]['paper_title'][:60]}", flush=True)
        time.sleep(0.5)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
