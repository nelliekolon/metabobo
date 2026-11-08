#!/usr/bin/env python3
"""Survey: which of the 136 candidate studies have (a) a findable paper with open full text,
(b) repository-declared method metadata, for an LLM metadata-extraction eval.
usage: survey_papers.py <pilot_slice_136_datasets.csv> <out.tsv>"""
import csv, json, re, sys, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor

UA = {"User-Agent": "metabobo-survey/0.1 (academic research script)"}


def get(url, as_json=True, tries=3):
    for t in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                b = r.read().decode("utf-8", "replace")
                return json.loads(b) if as_json else b
        except Exception as e:
            err = e
            time.sleep(3 * (t + 1))
    return {"_error": str(err)}


def europepmc(acc):
    # papers whose full text (or abstract) mentions the accession, e.g. in a data availability statement
    q = urllib.parse.quote(f'"{acc}"')
    d = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json&pageSize=25&resultType=lite")
    res = d.get("resultList", {}).get("result", []) if isinstance(d, dict) else []
    oa = [r for r in res if r.get("isOpenAccess") == "Y" or r.get("inPMC") == "Y"]
    return len(res), len(oa), ";".join(r.get("doi", r.get("id", "")) for r in res[:3])


def repo_info(acc, source):
    out = {"repo_paper": "", "repo_method_fields": 0}
    if acc.startswith("MTBLS"):
        d = get(f"https://www.ebi.ac.uk/metabolights/ws/studies/public/study/{acc}")
        c = d.get("content", {}) if isinstance(d, dict) else {}
        pubs = [p.get("doi") or p.get("pubmedId") or p.get("title", "")[:40] for p in (c.get("publications") or [])]
        out["repo_paper"] = ";".join(x for x in pubs if x)
        protos = c.get("protocols") or []
        out["repo_method_fields"] = sum(bool((p.get("description") or "").strip()) for p in protos)
    elif acc.startswith("ST"):
        a = get(f"https://www.metabolomicsworkbench.org/rest/study/study_id/{acc}/analysis")
        if isinstance(a, dict) and "analysis_id" in a:
            a = {"1": a}
        ids = [v["analysis_id"] for v in a.values() if isinstance(v, dict) and "analysis_id" in v] if isinstance(a, dict) else []
        if ids:
            t = get(f"https://www.metabolomicsworkbench.org/rest/study/analysis_id/{ids[0]}/mwtab/txt", as_json=False)
            if isinstance(t, str):
                m = re.findall(r"^(?:PR|ST):(?:PUBLICATIONS|DOI)\s+(.+)$", t, re.M)
                out["repo_paper"] = ";".join(x.strip()[:80] for x in m)
                keys = ["CH:COLUMN_NAME", "CH:FLOW_GRADIENT", "MS:ION_MODE", "MS:INSTRUMENT_NAME", "MS:MS_COMMENTS",
                        "CH:CHROMATOGRAPHY_TYPE"]
                out["repo_method_fields"] = sum(bool(re.search(rf"^{k}\s+\S", t, re.M)) for k in keys)
        out["n_analyses"] = len(ids)
    elif acc.startswith("MSV"):
        d = get(f"https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/{acc}")
        pubs = d.get("publications", []) if isinstance(d, dict) else []
        out["repo_paper"] = ";".join(p.get("value") or "" for p in pubs
                                     if "no associated" not in p.get("name", "").lower())
        out["repo_method_fields"] = int(bool(isinstance(d, dict) and d.get("summary")))
    return out


def one(row):
    acc = row["ATTRIBUTE_DatasetAccession"]
    n, n_oa, ids = europepmc(acc)
    r = repo_info(acc, row["source"])
    time.sleep(0.5)
    return dict(accession=acc, source=row["source"], year=row["year"], fluid=row["fluid"],
                instrument=row["instrument"], files=row["files"],
                repo_paper=r["repo_paper"], repo_method_fields=r["repo_method_fields"],
                n_analyses=r.get("n_analyses", ""), epmc_hits=n, epmc_open_access=n_oa, epmc_ids=ids)


def main():
    rows = list(csv.DictReader(open(sys.argv[1])))
    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(one, rows))
    with open(sys.argv[2], "w", newline="") as fh:
        w = csv.DictWriter(fh, list(res[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(res)


if __name__ == "__main__":
    main()
