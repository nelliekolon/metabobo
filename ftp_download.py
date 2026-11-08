#!/usr/bin/env python3
"""Download files in a Pan-ReDU USI list straight from the repositories' own servers.

  MassIVE / GNPS (MSV...)  -> https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.<dataset>/<path>
  MetaboLights (MTBLS...)  -> ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/<dataset>/<path>

MassIVE uses HTTPS because its FTP server refuses the ccms_peak symlinks most USIs go through.

Usage:
  python3 ftp_download.py <usi_list.tsv> <output_folder> [--only DATASET ...]

Files land in <output_folder>/<dataset>/<filename>. Files already there are skipped,
so it is safe to rerun after an interruption. A log is written to <output_folder>/ftp_log.tsv.
Uses only the Python standard library.
"""
import csv, ftplib, os, shutil, sys, time, urllib.error, urllib.parse, urllib.request

MASSIVE_HOST = "massive.ucsd.edu"
MASSIVE_URL = f"https://{MASSIVE_HOST}/ProteoSAFe/DownloadResultFile?forceDownload=true&file="
EBI_HOST = "ftp.ebi.ac.uk"
EBI_ROOT = "/pub/databases/metabolights/studies/public"
MAX_ATTEMPTS = 8


def parse_usi(usi):
    # mzspec:<dataset>:<path>[:scan:N]
    parts = usi.split(":", 2)
    dataset, path = parts[1], parts[2]
    if ":scan:" in path:
        path = path.split(":scan:")[0]
    return dataset, path


def connect(host):
    ftp = ftplib.FTP(host, timeout=120)
    ftp.login()
    ftp.set_pasv(True)
    return ftp


def fetch(ftp, remote, local):
    # Resumes from any .part left by a stalled attempt or an earlier run
    tmp = local + ".part"
    ftp.voidcmd("TYPE I")
    total = ftp.size(remote)
    offset = os.path.getsize(tmp) if os.path.exists(tmp) else 0
    if offset > total:
        os.remove(tmp)
        offset = 0
    if offset < total:
        try:
            with open(tmp, "ab") as fh:
                ftp.retrbinary(f"RETR {remote}", fh.write, blocksize=1 << 20, rest=offset or None)
        except ftplib.error_perm:
            raise
        except Exception as e:
            size = os.path.getsize(tmp)
            raise IOError(f"{e} (started at {offset/1e6:.0f}, reached {size/1e6:.0f} of {total/1e6:.0f} MB)") from e
    size = os.path.getsize(tmp)
    if size != total:
        raise IOError(f"incomplete: got {size} of {total} bytes")
    os.replace(tmp, local)


def fetch_https(dataset, path, local):
    url = MASSIVE_URL + urllib.parse.quote(f"f.{dataset}/{path}", safe="/")
    tmp = local + ".part"
    with urllib.request.urlopen(url, timeout=120) as resp, open(tmp, "wb") as fh:
        expected = resp.headers.get("Content-Length")
        shutil.copyfileobj(resp, fh, 1 << 20)
    size = os.path.getsize(tmp)
    if expected and size != int(expected):
        raise IOError(f"truncated: got {size} of {expected} bytes")
    # MassIVE usually sends no Content-Length, and a partial file would be skipped on every rerun
    if path.endswith(".mzML"):
        with open(tmp, "rb") as fh:
            fh.seek(max(0, size - 64))
            if b"mzML>" not in fh.read():
                raise IOError(f"truncated: no closing mzML tag after {size} bytes")
    os.replace(tmp, local)


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        sys.exit(__doc__)
    usi_file, out_dir = args[0], args[1]
    only = set(args[args.index("--only") + 1:]) if "--only" in args else None

    with open(usi_file) as fh:
        usis = [row["usi"] for row in csv.DictReader(fh, delimiter="\t") if row.get("usi")]

    by_dataset = {}
    for u in usis:
        ds, path = parse_usi(u)
        if only and ds not in only:
            continue
        by_dataset.setdefault(ds, []).append(path)

    log_path = os.path.join(out_dir, "ftp_log.tsv")
    os.makedirs(out_dir, exist_ok=True)
    new_log = not os.path.exists(log_path)
    log = open(log_path, "a")
    if new_log:
        log.write("dataset\tpath\tstatus\tbytes\n")

    for ds, paths in by_dataset.items():
        ds_dir = os.path.join(out_dir, ds)
        os.makedirs(ds_dir, exist_ok=True)
        if ds.startswith("MSV"):
            host = MASSIVE_HOST
        elif ds.startswith("MTBLS"):
            host = EBI_HOST
        else:
            print(f"Skipping {ds}: only MassIVE and MetaboLights are handled here")
            continue

        print(f"\n== {ds}: {len(paths)} files from {host}")
        ftp = None

        for i, path in enumerate(paths, 1):
            local = os.path.join(ds_dir, os.path.basename(path))
            if os.path.exists(local) and os.path.getsize(local) > 10000:
                log.write(f"{ds}\t{path}\tEXISTS\t{os.path.getsize(local)}\n")
                continue
            status = "FAILED"
            for attempt in range(MAX_ATTEMPTS):
                try:
                    if host == MASSIVE_HOST:
                        fetch_https(ds, path, local)
                    else:
                        ftp = ftp or connect(host)
                        fetch(ftp, f"{EBI_ROOT}/{ds}/{path}", local)
                    size = os.path.getsize(local)
                    print(f"[{i}/{len(paths)}] ok   {os.path.basename(path)} ({size/1e6:.0f} MB)")
                    log.write(f"{ds}\t{path}\tOK\t{size}\n")
                    break
                except ftplib.error_perm as e:
                    print(f"[{i}/{len(paths)}] MISSING {path}: {e}")
                    log.write(f"{ds}\t{path}\tNOT_FOUND\t0\n")
                    break
                except Exception as e:
                    print(f"[{i}/{len(paths)}] retry {attempt+1} {os.path.basename(path)}: {e}")
                    # MassIVE answers a missing file with HTTP 500, which could also be a passing error
                    if isinstance(e, urllib.error.HTTPError) and e.code in (404, 500):
                        status = "NOT_FOUND"
                    is_429 = isinstance(e, urllib.error.HTTPError) and e.code == 429
                    time.sleep(60 if is_429 else 10)
                    try:
                        ftp.close()
                    except Exception:
                        pass
                    ftp = None
            else:
                print(f"[{i}/{len(paths)}] {status} {path}")
                log.write(f"{ds}\t{path}\t{status}\t0\n")
            if host == MASSIVE_HOST and os.path.exists(local + ".part"):
                os.remove(local + ".part")
            log.flush()
        try:
            ftp.quit()
        except Exception:
            pass
    log.close()
    print(f"\nDone. Log: {log_path}")


if __name__ == "__main__":
    main()
