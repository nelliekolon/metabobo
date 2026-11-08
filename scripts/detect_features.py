#!/usr/bin/env python3
"""Step 1: untargeted feature detection on centroided mzML files with OpenMS.

usage: python3 detect_features.py <in_dir> <out_dir> [n_files]
Writes one .featureXML per input file into out_dir; skips files already done.
Fixed parameters (keep identical across every study so results are comparable).
"""
import sys, os, time, glob
import pyopenms as oms

PARAMS = {
    "mtd": {"mass_error_ppm": 5.0, "noise_threshold_int": 3000.0, "chrom_peak_snr": 3.0,
            "min_trace_length": 3.0},
    "epd": {"width_filtering": "fixed", "min_fwhm": 1.0, "max_fwhm": 60.0},
    "ffm": {"isotope_filtering_model": "none", "remove_single_traces": "false",
            "mz_scoring_by_elements": "false", "report_convex_hulls": "false",
            "local_rt_range": 10.0, "local_mz_range": 6.5, "charge_lower_bound": 1,
            "charge_upper_bound": 3},
}

def setp(algo, d):
    p = algo.getDefaults()
    for k, v in d.items():
        p.setValue(k, v)
    algo.setParameters(p)

def detect(path, out):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(path, exp)
    exp.setSpectra([s for s in exp.getSpectra() if s.getMSLevel() == 1])
    exp.sortSpectra(True)
    traces = []
    mtd = oms.MassTraceDetection(); setp(mtd, PARAMS["mtd"]); mtd.run(exp, traces, 0)
    split = []
    epd = oms.ElutionPeakDetection(); setp(epd, PARAMS["epd"]); epd.detectPeaks(traces, split)
    fm = oms.FeatureMap(); chroms = []
    ffm = oms.FeatureFindingMetabo(); setp(ffm, PARAMS["ffm"]); ffm.run(split, fm, chroms)
    fm.setUniqueIds()
    fm.setPrimaryMSRunPath([path.encode()])
    oms.FeatureXMLFile().store(out, fm)
    return len(traces), fm.size()

if __name__ == "__main__":
    in_dir, out_dir = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else None
    os.makedirs(out_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(in_dir, "*.mzML")))[:n]
    for i, f in enumerate(files, 1):
        out = os.path.join(out_dir, os.path.basename(f).replace(".mzML", ".featureXML"))
        if os.path.exists(out):
            continue
        t = time.time()
        nt, nf = detect(f, out)
        print(f"[{i}/{len(files)}] {os.path.basename(f)}: {nt} traces, {nf} features, {time.time()-t:.0f}s", flush=True)
