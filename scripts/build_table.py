#!/usr/bin/env python3
"""Step 2: align retention times across runs, link features, export a sample x feature table.

usage: python3 build_table.py <feature_dir> <out_dir> [--no-align] [--samples samples.tsv]
Outputs: consensus.consensusXML, feature_table.tsv (features x samples, intensities),
         feature_info.tsv (mz, rt, how many samples detected in).
"""
import sys, os, glob
import pyopenms as oms
import pandas as pd

feat_dir, out_dir = sys.argv[1], sys.argv[2]
ALIGN = "--no-align" not in sys.argv
os.makedirs(out_dir, exist_ok=True)
files = sorted(glob.glob(os.path.join(feat_dir, "*.featureXML")))
# --samples <samples.tsv>: keep only runs marked include=1 (see make_samples.py)
if "--samples" in sys.argv:
    sheet = pd.read_csv(sys.argv[sys.argv.index("--samples") + 1], sep="\t", dtype=str)
    keep = set(sheet.loc[sheet.include == "1", "sample"])
    files = [f for f in files if os.path.basename(f).replace(".featureXML", "") in keep]
    if missing := keep - {os.path.basename(f).replace(".featureXML", "") for f in files}:
        print(f"warning: {len(missing)} included samples have no featureXML, e.g. {sorted(missing)[:3]}")
maps = []
for f in files:
    fm = oms.FeatureMap(); oms.FeatureXMLFile().load(f, fm); maps.append(fm)
names = [os.path.basename(f).replace(".featureXML", "") for f in files]

# RT alignment against the run with the most features
ref = max(range(len(maps)), key=lambda i: maps[i].size())
aligner = oms.MapAlignmentAlgorithmPoseClustering()
p = aligner.getDefaults()
p.setValue("max_num_peaks_considered", 1000)
p.setValue("pairfinder:distance_MZ:max_difference", 10.0)
p.setValue("pairfinder:distance_MZ:unit", "ppm")
p.setValue("pairfinder:distance_RT:max_difference", 30.0)
aligner.setParameters(p)
aligner.setReference(maps[ref])
shifts = []
for i, fm in enumerate(maps):
    if i == ref or not ALIGN:
        continue
    trafo = oms.TransformationDescription()
    aligner.align(fm, trafo)
    tr = oms.MapAlignmentTransformer()
    tr.transformRetentionTimes(fm, trafo, True)

# Link features across runs (QT clustering)
grouper = oms.FeatureGroupingAlgorithmQT()
gp = grouper.getDefaults()
gp.setValue("distance_RT:max_difference", 15.0)
gp.setValue("distance_MZ:max_difference", 10.0)
gp.setValue("distance_MZ:unit", "ppm")
grouper.setParameters(gp)
cmap = oms.ConsensusMap()
headers = cmap.getColumnHeaders()
for i, (fm, f) in enumerate(zip(maps, files)):
    h = oms.ColumnHeader(); h.filename = f; h.label = names[i]; h.size = fm.size(); headers[i] = h
cmap.setColumnHeaders(headers)
grouper.group(maps, cmap)
cmap.setUniqueIds()
oms.ConsensusXMLFile().store(os.path.join(out_dir, "consensus.consensusXML"), cmap)

rows, info = [], []
for k, cf in enumerate(cmap):
    vals = {names[h.getMapIndex()]: h.getIntensity() for h in cf.getFeatureList()}
    fid = f"F{k:05d}_{cf.getMZ():.4f}_{cf.getRT():.1f}"
    rows.append(pd.Series(vals, name=fid))
    info.append({"feature": fid, "mz": cf.getMZ(), "rt_s": cf.getRT(), "n_detected": len(vals)})
tab = pd.DataFrame(rows).reindex(columns=names)
tab.index.name = "feature"
tab.to_csv(os.path.join(out_dir, "feature_table.tsv"), sep="\t")
pd.DataFrame(info).to_csv(os.path.join(out_dir, "feature_info.tsv"), sep="\t", index=False)
nd = pd.DataFrame(info).n_detected
print(f"{len(names)} samples, {len(tab)} consensus features")
for t in [1, 0.5, 0.8]:
    thr = max(1, int(round(t * len(names)))) if t < 1 else len(names)
    print(f"  detected in >= {thr} samples: {(nd >= thr).sum()}")
