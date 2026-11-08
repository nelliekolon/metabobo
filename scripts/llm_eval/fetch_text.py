#!/usr/bin/env python3
"""Eval step 2: fetch the text each model condition reads (see SPEC.md).

usage: python3 fetch_text.py <papers.tsv> <out_dir>
papers.tsv: study, doi, pmcid, paper_role (original / companion / none), note.
Writes <out_dir>/<study>/:
  fulltext.txt    open-access paper body from Europe PMC (no references), if the study has a paper
  methods.txt     the paper's methods sections only
  repository.txt  the repository's own description, protocols and structured fields
  sources.json    URLs fetched and what was found
Paper text is kept out of git (licences vary); only this script and papers.tsv are tracked.
"""
import csv, json, os, re, sys, time, urllib.request
import xml.etree.ElementTree as ET

UA = {"User-Agent": "metabobo-eval/0.1 (academic research script)"}
METHODS = re.compile(r"method|material|experimental|procedure|sample (preparation|collection)|"
                     r"lc[- ]?ms|mass spectrometr|metabolom|data (acquisition|analysis)", re.I)
NOT_METHODS = re.compile(r"result|discussion|introduction|conclusion|acknowledg|author|funding|"
                         r"competing|availability of|supplementary|reference", re.I)


def get(url, as_json=False):
    for t in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
                b = r.read()
                return json.loads(b) if as_json else b.decode("utf-8", "replace")
        except Exception as e:
            err = e; time.sleep(3 * (t + 1))
    raise err


def text_of(el):
    """Readable text of an element: paragraphs on their own lines, tables kept as plain rows."""
    out = []
    for node in el.iter():
        tag = node.tag.split("}")[-1]
        if tag in ("title",) and node.text:
            out.append("\n## " + "".join(node.itertext()).strip())
        elif tag == "p":
            out.append("".join(node.itertext()).strip())
        elif tag == "tr":
            out.append(" | ".join("".join(c.itertext()).strip() for c in node))
    txt = "\n".join(x for x in out if x)
    return re.sub(r"[ \t]+", " ", txt)


def top_sections(body):
    return [s for s in body if s.tag.split("}")[-1] == "sec"]


def paper(pmcid):
    # Europe PMC first; some recent articles 500 there, so fall back to NCBI PMC (same JATS XML)
    urls = [f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML",
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid.removeprefix('PMC')}"]
    for url in urls:
        try:
            root = ET.fromstring(get(url))
            if root.find(".//body") is not None:
                break
        except Exception as e:
            print(f"  {url}: {e}", flush=True)
    else:
        raise ValueError(f"no full text for {pmcid}")
    body = root.find(".//body")
    if body is None:
        raise ValueError("no body in full text XML")
    full = text_of(body)
    meth = []
    for sec in top_sections(body):
        title = "".join((sec.find("title").itertext() if sec.find("title") is not None else [])).strip()
        stype = sec.get("sec-type", "")
        if "method" in stype or (METHODS.search(title) and not NOT_METHODS.search(title)):
            meth.append(text_of(sec))
    # Data in Brief and similar papers put methods under "Experimental design, materials and methods"
    # inside nested sections; fall back to any nested section with a methods-like title
    if not meth:
        for sec in body.iter():
            if sec.tag.split("}")[-1] != "sec":
                continue
            t = sec.find("title")
            title = "".join(t.itertext()).strip() if t is not None else ""
            if METHODS.search(title) and not NOT_METHODS.search(title):
                meth.append(text_of(sec))
    return url, full, "\n\n".join(meth)


def repository(study):
    if study.startswith("MTBLS"):
        url = f"https://www.ebi.ac.uk/metabolights/ws/studies/public/study/{study}"
        c = get(url, as_json=True).get("content", {})
        parts = [f"# {c.get('title', '')}", c.get("description", "")]
        for p in c.get("protocols") or []:
            parts.append(f"\n## Protocol: {p.get('name', '')}\n{p.get('description', '')}")
        for a in c.get("assays") or []:
            parts.append(f"\n## Assay: technology {a.get('technology', '')}; platform {a.get('platform', '')}; "
                         f"measurement {a.get('measurement', '')}")
        return url, "\n".join(x for x in parts if x)
    if study.startswith("MSV"):
        url = f"https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/{study}"
        d = get(url, as_json=True)
        parts = [f"# {d.get('title', '')}", d.get("summary", "")]
        parts.append("Instruments: " + ", ".join(i.get("name", "") for i in d.get("instruments", [])))
        parts.append("Keywords: " + ", ".join(k.get("value", "") for k in d.get("keywords", []) if k.get("value")))
        return url, "\n".join(parts)
    if study.startswith("ST"):
        url = f"https://www.metabolomicsworkbench.org/rest/study/study_id/{study}/analysis"
        a = get(url, as_json=True); a = {"1": a} if "analysis_id" in a else a
        parts, urls = [], [url]
        for v in a.values():
            u = f"https://www.metabolomicsworkbench.org/rest/study/analysis_id/{v['analysis_id']}/mwtab/txt"
            urls.append(u); parts.append(get(u))
        return ";".join(urls), "\n\n".join(parts)
    raise ValueError(study)


def main():
    papers, out_dir = sys.argv[1], sys.argv[2]
    for r in csv.DictReader(open(papers), delimiter="\t"):
        s = r["study"]; d = os.path.join(out_dir, s); os.makedirs(d, exist_ok=True)
        src = {"study": s, "paper_role": r["paper_role"]}
        rurl, rtxt = repository(s)
        open(os.path.join(d, "repository.txt"), "w").write(rtxt)
        src["repository_url"] = rurl; src["repository_chars"] = len(rtxt)
        if r["pmcid"]:
            purl, full, meth = paper(r["pmcid"])
            open(os.path.join(d, "fulltext.txt"), "w").write(full)
            open(os.path.join(d, "methods.txt"), "w").write(meth)
            src.update(paper_url=purl, doi=r["doi"], fulltext_chars=len(full), methods_chars=len(meth))
        json.dump(src, open(os.path.join(d, "sources.json"), "w"), indent=1)
        print(src, flush=True)
        time.sleep(1)


if __name__ == "__main__":
    main()
