#!/usr/bin/env python3
"""Step 0b (MetaboLights only): convert Thermo .raw files to centroided MS1 mzML.

usage: python3 convert_raw.py <raw_dir> <mzml_dir> [--samples samples.tsv]
Uses ThermoRawFileParser (tools/ThermoRawFileParser, v2.0.0-dev osx-arm64 build).
Output matches the MassIVE files the pipeline already reads: mzML, MS1 only, vendor
centroiding ON (detect_features.py expects centroided spectra).
Skips files already converted; with --samples, only converts runs with include=1.
"""
import csv, glob, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
TRFP = os.path.join(HERE, "..", "tools", "ThermoRawFileParser", "ThermoRawFileParser")


def main():
    raw_dir, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    raws = sorted(glob.glob(os.path.join(raw_dir, "*.raw")))
    if "--samples" in sys.argv:
        with open(sys.argv[sys.argv.index("--samples") + 1]) as fh:
            keep = {r["sample"] for r in csv.DictReader(fh, delimiter="\t") if r["include"] == "1"}
        raws = [r for r in raws if os.path.basename(r)[:-4] in keep]
    for i, raw in enumerate(raws, 1):
        name = os.path.basename(raw)[:-4]
        out = os.path.join(out_dir, name + ".mzML")
        if os.path.exists(out):
            continue
        t = time.time()
        # -f=1 mzML (not indexed), -L=1 MS1 only; peak picking is on unless -p is given
        res = subprocess.run([TRFP, f"-i={raw}", f"-o={out_dir}", "-f=1", "-L=1"],
                             capture_output=True, text=True)
        if res.returncode != 0 or not os.path.exists(out):
            print(f"[{i}/{len(raws)}] FAILED {name}: {(res.stderr or res.stdout).strip()[-300:]}", flush=True)
            continue
        with open(out, "rb") as fh:            # same truncation check as ftp_download.py
            fh.seek(max(0, os.path.getsize(out) - 64))
            if b"mzML>" not in fh.read():
                os.remove(out)
                print(f"[{i}/{len(raws)}] FAILED {name}: truncated mzML", flush=True)
                continue
        print(f"[{i}/{len(raws)}] ok {name} ({os.path.getsize(out)/1e6:.0f} MB, {time.time()-t:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
