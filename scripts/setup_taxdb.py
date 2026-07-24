#!/usr/bin/env python3
"""
setup_taxdb.py — build a PINNED ete3 NCBI-taxonomy database for reproducible
name/taxid resolution in profile2cami.py.

Why pin: ete3's default NCBITaxa() silently downloads whatever taxdump is current
that day, so results drift over time and older builds can miss renamed phyla
(Pseudomonadota, Bacillota, …). This script instead builds the ete3 SQLite from a
*specific* NCBI taxdump you choose, records its URL + SHA-256 + timestamp, and
checks that the updated phylum names resolve — so the benchmark is reproducible
and the FOCUS output maps cleanly.

Pinning knob (--url / $TAXDUMP_URL), most-pinned first:
  1. A dated archive dump (recommended, fully reproducible):
       https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_<YYYY-MM-DD>.zip
     (list available dates: browse .../pub/taxonomy/taxdump_archive/)
  2. The current dump (pins by content hash only, date = whatever is live today):
       https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz   [default]

Both .zip (dated archive) and .tar.gz (current/new_taxdump) inputs are accepted.
The built SQLite + the written taxdb.version file are the artifacts to keep/commit
alongside the analysis so anyone can reproduce the mapping.

Usage:
  python3 setup_taxdb.py --url <dump-url> -o ./taxdb/taxa.sqlite
  export PROFILE2CAMI_TAXDB=$PWD/taxdb/taxa.sqlite   # profile2cami.py picks it up
"""

import argparse
import datetime
import hashlib
import io
import json
import os
import shutil
import sys
import tarfile
import tempfile
import urllib.request
import zipfile

DEFAULT_URL = os.environ.get(
    "TAXDUMP_URL",
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")

# Updated phylum names that must resolve for FOCUS output to map cleanly.
CHECK_NAMES = ["Pseudomonadota", "Bacillota", "Actinomycetota", "Bacteroidota",
               "Campylobacterota", "Cyanobacteriota", "Thermodesulfobacteriota",
               "Bdellovibrionota", "Synergistota", "Planctomycetota"]

NEEDED_DMP = ["nodes.dmp", "names.dmp"]          # minimum ete3 needs
EXTRA_DMP = ["merged.dmp", "delnodes.dmp"]       # include if present


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_source(src, dest):
    """Return (path_to_dump, last_modified). `src` may be a LOCAL path (used in
    place, no copy/download) or an http(s)/ftp URL (downloaded to dest)."""
    if os.path.exists(src):
        mtime = datetime.datetime.utcfromtimestamp(os.path.getmtime(src))
        print("[setup_taxdb] using local dump %s (mtime %sZ)" % (src, mtime.isoformat()),
              file=sys.stderr)
        return src, mtime.isoformat() + "Z"
    print("[setup_taxdb] downloading %s" % src, file=sys.stderr)
    req = urllib.request.Request(src, headers={"User-Agent": "setup_taxdb/1.0"})
    with urllib.request.urlopen(req) as resp, open(dest, "wb") as out:
        last_mod = resp.headers.get("Last-Modified")
        shutil.copyfileobj(resp, out)
    return dest, last_mod


def normalize_to_targz(src, targz):
    """Ensure a taxdump.tar.gz whose members include nodes.dmp/names.dmp at the
    archive root (what ete3 expects). Accepts a .tar.gz or a dated .zip."""
    if src.endswith(".zip") or zipfile.is_zipfile(src):
        with tempfile.TemporaryDirectory() as tmp:
            with zipfile.ZipFile(src) as zf:
                members = {os.path.basename(n): n for n in zf.namelist()}
                for dmp in NEEDED_DMP:
                    if dmp not in members:
                        raise SystemExit("ERROR: %s not found inside %s" % (dmp, src))
                with tarfile.open(targz, "w:gz") as tar:
                    for dmp in NEEDED_DMP + EXTRA_DMP:
                        if dmp in members:
                            zf.extract(members[dmp], tmp)
                            tar.add(os.path.join(tmp, members[dmp]), arcname=dmp)
        return
    # already a tar.gz: verify it carries the needed members; repack flat if nested
    with tarfile.open(src) as tin:
        names = {os.path.basename(m.name): m for m in tin.getmembers()}
        for dmp in NEEDED_DMP:
            if dmp not in names:
                raise SystemExit("ERROR: %s not found inside %s" % (dmp, src))
        flat = all(names[d].name == d for d in NEEDED_DMP)
    if flat:
        shutil.copyfile(src, targz)
    else:                                        # repack members at root
        with tarfile.open(src) as tin, tarfile.open(targz, "w:gz") as tout:
            for dmp in NEEDED_DMP + EXTRA_DMP:
                if dmp in names:
                    m = names[dmp]
                    f = tin.extractfile(m)
                    info = tarfile.TarInfo(name=dmp)
                    data = f.read()
                    info.size = len(data)
                    tout.addfile(info, io.BytesIO(data))


def build_sqlite(targz, dbfile):
    try:
        from ete3 import NCBITaxa
        import ete3
    except ImportError:
        raise SystemExit("ERROR: ete3 is required (pip install ete3).")
    if os.path.exists(dbfile):
        os.remove(dbfile)
    os.makedirs(os.path.dirname(os.path.abspath(dbfile)), exist_ok=True)
    print("[setup_taxdb] building ete3 SQLite at %s" % dbfile, file=sys.stderr)
    ncbi = NCBITaxa(dbfile=dbfile, taxdump_file=targz)
    return ncbi, getattr(ete3, "__version__", "unknown")


def verify(ncbi):
    got = ncbi.get_name_translator(CHECK_NAMES)
    resolved = {n: (got[n][0] if n in got else None) for n in CHECK_NAMES}
    n_ok = sum(1 for v in resolved.values() if v)
    for n, tid in resolved.items():
        print("  %-28s %s" % (n, ("taxid %s" % tid) if tid else "NOT FOUND"),
              file=sys.stderr)
    print("[setup_taxdb] updated-phylum check: %d/%d resolved"
          % (n_ok, len(CHECK_NAMES)), file=sys.stderr)
    return resolved, n_ok


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-u", "--url", default=DEFAULT_URL, metavar="SRC",
                    help="taxdump to pin: a LOCAL path (used in place, no download) "
                         "or an http(s)/ftp URL (.tar.gz or dated .zip). "
                         "Default: $TAXDUMP_URL or the current NCBI taxdump.tar.gz")
    ap.add_argument("-o", "--output", default="./taxdb/taxa.sqlite",
                    help="path for the pinned ete3 SQLite (default: ./taxdb/taxa.sqlite)")
    ap.add_argument("--cache", default=None,
                    help="dir to keep the downloaded dump (default: alongside output)")
    args = ap.parse_args()

    outdir = os.path.dirname(os.path.abspath(args.output))
    cache = args.cache or outdir
    os.makedirs(cache, exist_ok=True)

    fname = os.path.basename(args.url.split("?")[0]) or "taxdump.download"
    raw, last_mod = fetch_source(args.url, os.path.join(cache, fname))
    raw_sha = sha256(raw)
    print("[setup_taxdb] source sha256=%s  Last-Modified=%s"
          % (raw_sha, last_mod), file=sys.stderr)

    targz = os.path.join(cache, "taxdump.pinned.tar.gz")
    normalize_to_targz(raw, targz)

    ncbi, ete3_version = build_sqlite(targz, args.output)
    resolved, n_ok = verify(ncbi)

    version = {
        "source_url": args.url,
        "source_sha256": raw_sha,
        "source_last_modified": last_mod,
        "built_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "sqlite_path": os.path.abspath(args.output),
        "sqlite_sha256": sha256(args.output),
        "ete3_version": ete3_version,
        "updated_phylum_check": {n: resolved[n] for n in CHECK_NAMES},
        "updated_phylum_resolved": "%d/%d" % (n_ok, len(CHECK_NAMES)),
    }
    vpath = os.path.join(outdir, "taxdb.version")
    with open(vpath, "w") as fh:
        json.dump(version, fh, indent=2)
    print("[setup_taxdb] wrote %s" % vpath, file=sys.stderr)
    print("\nPin recorded. To use it everywhere:\n"
          "  export PROFILE2CAMI_TAXDB=%s" % os.path.abspath(args.output),
          file=sys.stderr)
    if n_ok < len(CHECK_NAMES):
        print("[setup_taxdb] WARNING: some updated phylum names did not resolve — "
              "choose a newer dated dump via --url.", file=sys.stderr)


if __name__ == "__main__":
    main()
