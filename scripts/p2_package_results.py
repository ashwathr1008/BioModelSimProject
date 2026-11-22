"""Package the completed X fits and their inputs as verified GitHub release assets."""

import argparse
import hashlib
import json
import shutil
from datetime import date
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

from filelock import FileLock

from p2_fit import input_path, load_config, manifest_path, repo_path, run_identity, sha256
from p2_run_all import configured_pcas, make_jobs


ARTIFACTS = ("run.json", "training.csv", "fitted.h5ad", "model/model.pt")
MANIFEST_PATH = "results/p2/artifact_manifest.json"


def file_record(path, root):
    path = path.resolve()
    relative = path.relative_to(root).as_posix()
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"Missing or empty artifact: {relative}")
    return {"path": relative, "bytes": path.stat().st_size, "sha256": sha256(path)}


def describe_assets(cfg, config, root, release_date):
    if configured_pcas(cfg) != ["X"]:
        raise ValueError("This release packages only the selected X experiment")
    fits, inputs = [], []
    jobs = make_jobs(cfg["sampling"]["n_repeats"], ["X"])
    for name, condition, repeat, pca in jobs:
        pca = "shipped" if condition == "tutorial" else pca
        folder = root / "results/p2" / name
        record = json.loads((folder / "run.json").read_text(encoding="utf-8"))
        if record.get("status") != "complete":
            raise ValueError(f"Incomplete fit: {name}")
        if record.get("identity") != run_identity(cfg, condition, repeat, pca):
            raise ValueError(f"Input or scientific settings mismatch: {name}")
        fits.extend(folder / artifact for artifact in ARTIFACTS)
        inputs.append(input_path(cfg, condition, repeat, pca))
        if condition != "tutorial":
            inputs.append(manifest_path(cfg, condition, repeat))
    raw = repo_path(cfg["data"]["source_file"])
    if sha256(raw) != cfg["data"]["expected_sha256"]:
        raise ValueError("Source data checksum does not match the configuration")
    archives = {}
    for kind, paths in (("fits", fits), ("inputs", inputs)):
        records = [file_record(path, root) for path in sorted(set(paths))]
        archives[f"regvelo-p2-x-{kind}-{release_date}.zip"] = {
            "kind": kind,
            "file_count": len(records),
            "uncompressed_bytes": sum(item["bytes"] for item in records),
            "files": records,
        }
    return {
        "schema_version": 1,
        "release_date": release_date,
        "configuration": {
            "path": config.relative_to(root).as_posix(),
            "sha256": sha256(config),
        },
        "pca_variants": ["X"],
        "main_fit_count": len(jobs) - 1,
        "tutorial_fit_count": 1,
        "extraction": "Extract into the repository root to restore repository-relative paths.",
        "archives": archives,
    }


def zip_entry(path, release_date):
    day = date.fromisoformat(release_date)
    entry = ZipInfo(path, (day.year, day.month, day.day, 0, 0, 0))
    entry.compress_type = ZIP_DEFLATED
    entry.create_system = 3
    entry.external_attr = 0o100644 << 16
    entry._compresslevel = 1
    return entry


def verify_archive(path, records, manifest_bytes):
    expected = {item["path"] for item in records} | {MANIFEST_PATH}
    with ZipFile(path) as archive:
        if len(archive.namelist()) != len(expected) or set(archive.namelist()) != expected:
            raise ValueError(f"Unexpected members in {path.name}")
        if archive.read(MANIFEST_PATH) != manifest_bytes:
            raise ValueError(f"Manifest differs in {path.name}")
        # Reading every member checks its CRC as well as the recorded SHA-256.
        for item in records:
            digest = hashlib.sha256()
            with archive.open(item["path"]) as member:
                for block in iter(lambda: member.read(4 << 20), b""):
                    digest.update(block)
            if (archive.getinfo(item["path"]).file_size != item["bytes"]
                    or digest.hexdigest() != item["sha256"]):
                raise ValueError(f"Archive content mismatch: {item['path']}")


def create_archive(path, records, manifest_bytes, root, release_date):
    if path.exists():
        verify_archive(path, records, manifest_bytes)
        print(f"Verified existing archive: {path.name}", flush=True)
        return
    temporary = path.with_suffix(".zip.partial")
    if temporary.exists():
        raise FileExistsError(f"Preserving incomplete prior archive: {temporary}")
    with ZipFile(temporary, "x", compression=ZIP_DEFLATED, compresslevel=1,
                 allowZip64=True) as archive:
        archive.writestr(zip_entry(MANIFEST_PATH, release_date), manifest_bytes)
        for item in records:
            source = root / item["path"]
            with source.open("rb") as stream, archive.open(
                    zip_entry(item["path"], release_date), "w", force_zip64=True) as member:
                shutil.copyfileobj(stream, member, length=4 << 20)
            print(f"Packed {item['path']}", flush=True)
    verify_archive(temporary, records, manifest_bytes)
    if temporary.stat().st_size >= 2 * 2**30:
        raise ValueError(f"Archive exceeds GitHub's per-asset size limit: {temporary}")
    temporary.rename(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=repo_path("config/p2_x.yaml"))
    parser.add_argument("--release-date", default=date.today().isoformat())
    parser.add_argument("--output-dir", type=Path, default=repo_path("results/p2/release_assets"))
    args = parser.parse_args()
    date.fromisoformat(args.release_date)
    root = repo_path(".").resolve()
    config = args.config.resolve()
    output = args.output_dir.resolve()
    output.relative_to(root)
    config.relative_to(root)
    output.mkdir(parents=True, exist_ok=True)
    with FileLock(output / ".package.lock", timeout=0):
        manifest = describe_assets(load_config(config), config, root, args.release_date)
        manifest_bytes = (json.dumps(manifest, indent=2) + "\n").encode("utf-8")
        manifest_path = root / MANIFEST_PATH
        if manifest_path.exists() and manifest_path.read_bytes() != manifest_bytes:
            raise ValueError("Existing artifact manifest differs; preserve the published release")
        manifest_path.write_bytes(manifest_bytes)
        checksums, summaries = [], []
        for name, asset in manifest["archives"].items():
            path = output / name
            create_archive(path, asset["files"], manifest_bytes, root, args.release_date)
            digest = sha256(path)
            checksums.append(f"{digest}  {name}")
            summaries.append({"name": name, "bytes": path.stat().st_size, "sha256": digest})
        checksums.append(f"{sha256(manifest_path)}  artifact_manifest.json")
        (output / "artifact_manifest.json").write_bytes(manifest_bytes)
        checksum_bytes = ("\n".join(checksums) + "\n").encode("utf-8")
        (output / "SHA256SUMS").write_bytes(checksum_bytes)
        (root / "results/p2/release_checksums.txt").write_bytes(checksum_bytes)
        (output / "package_summary.json").write_text(
            json.dumps(summaries, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(summaries, indent=2), flush=True)


if __name__ == "__main__":
    main()
