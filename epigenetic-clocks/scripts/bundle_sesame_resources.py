#!/usr/bin/env python3
"""Bundle only cached human-array resources; this command never downloads data.

Prepare the source cache first with preprocess_idat.R --init-cache --cache PATH.
Uses Rscript from the active R library to resolve the matching sesameData IDs.
"""
import argparse
import csv
import hashlib
import io
import json
import shutil
import sqlite3
import subprocess
import tarfile
import tempfile
from pathlib import Path


def bundle(cache, output, rscript="Rscript"):
    cache, output = Path(cache).resolve(), Path(output).resolve()
    rcode = '''
    d <- sesameData::sesameDataList()
    keep <- d$Title %in% c("idatSignature", paste0(c("HM450", "EPIC", "EPICv2", "MSA"), ".address")) |
        grepl("^KYCG[.](HM450|EPIC|EPICv2|MSA)[.]Mask[.]", d$Title)
    d <- as.data.frame(d[keep, c("EHID", "Title")])
    d$sesame <- as.character(packageVersion("sesame"))
    d$sesameData <- as.character(packageVersion("sesameData"))
    write.csv(d, stdout(), row.names=FALSE)
    '''
    result = subprocess.run([rscript, "-e", rcode], text=True, capture_output=True, check=True)
    resources = list(csv.DictReader(io.StringIO(result.stdout)))
    basic = {"idatSignature", *(f"{p}.address" for p in ("HM450", "EPIC", "EPICv2", "MSA"))}
    if not resources or not basic.issubset({r["Title"] for r in resources}):
        raise ValueError("Installed sesameData lacks required human-array resources")
    archive = output / "sesame-reference-cache.tar.gz"
    manifest_path = output / "sesame-resources.json"
    if archive.exists() or manifest_path.exists():
        raise ValueError("Bundle outputs already exist; use a fresh output directory")
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="sesame-bundle-") as temporary:
        stage = Path(temporary) / "cache"
        stage.mkdir()
        with sqlite3.connect(f"file:{cache / 'BiocFileCache.sqlite'}?mode=ro", uri=True) as src:
            with sqlite3.connect(stage / "BiocFileCache.sqlite") as dst:
                src.backup(dst)
        manifest = {"sesame": resources[0]["sesame"], "sesameData": resources[0]["sesameData"],
                    "source": "https://bioconductor.org/packages/sesameData/",
                    "resources": []}
        with sqlite3.connect(stage / "BiocFileCache.sqlite") as db:
            rows = db.execute("SELECT rid,rname,rpath FROM resource").fetchall()
            selected = []
            for title in ["experimenthub.sqlite3", "experimenthub.index.rds"] + [r["EHID"] for r in resources]:
                matches = [r for r in rows if r[1] == title or r[1].startswith(title + " : ")]
                if len(matches) != 1:
                    raise ValueError(f"Missing or ambiguous cached resource {title}; initialize the cache first")
                rid, rname, rpath = matches[0]
                source = Path(rpath) if Path(rpath).is_absolute() else cache / rpath
                if not source.is_file():
                    raise ValueError(f"Cached resource file missing: {title}")
                destination = stage / source.name
                shutil.copy2(source, destination)
                db.execute("UPDATE resource SET rpath=? WHERE rid=?", (source.name, rid))
                selected.append(rid)
                label = next((r["Title"] for r in resources if r["EHID"] == title), title)
                manifest["resources"].append({"id": title, "title": label, "file": source.name,
                    "bytes": destination.stat().st_size,
                    "sha256": hashlib.sha256(destination.read_bytes()).hexdigest()})
            placeholders = ",".join("?" for _ in selected)
            db.execute(f"DELETE FROM resource WHERE rid NOT IN ({placeholders})", selected)
            db.commit()
            db.execute("VACUUM")
        with (stage / "bundle_versions.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["package", "version"])
            writer.writerows((p, manifest[p]) for p in ("sesame", "sesameData"))
        with tarfile.open(archive, "w:gz") as tar:
            for path in sorted(stage.iterdir()):
                # Store public reference files only; no absolute source paths or user files.
                info = tar.gettarinfo(str(path), arcname=path.name)
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                with path.open("rb") as handle:
                    tar.addfile(info, handle)
        manifest["archive_bytes"] = archive.stat().st_size
        manifest["archive_sha256"] = hashlib.sha256(archive.read_bytes()).hexdigest()
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Bundled {len(resources)} SeSAMe resources in {archive} ({archive.stat().st_size / 2**20:.1f} MiB)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    bundle(args.cache, args.output_dir, args.rscript)
