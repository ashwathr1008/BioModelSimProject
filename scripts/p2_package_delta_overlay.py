"""Package completed local fits and updated Slurm code over the uploaded base."""

import hashlib
import json
import tarfile
import uuid
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TRANSFER = ROOT / "results/p2/delta_transfer"
ARTIFACTS = ("run.json", "training.csv", "fitted.h5ad", "model/model.pt")


def digest(stream):
    value = hashlib.sha256()
    for chunk in iter(lambda: stream.read(2**20), b""):
        value.update(chunk)
    return value.hexdigest()


def checksum(path):
    with path.open("rb") as stream:
        return digest(stream)


def main():
    base = json.loads((TRANSFER / "deployment.json").read_text())
    identifier = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    destination = TRANSFER / identifier
    destination.mkdir(parents=True, exist_ok=False)
    receipt_path = destination / "imported_local_full0.json"
    run = json.loads((ROOT / "results/p2/X/full_r0/run.json").read_text())
    assert run["status"] == "complete" and run["pca_input"] == "X" and run["repeat"] == 0
    import_files = [{"path": f"results/p2/X/full_r0/{name}",
                     "sha256": checksum(ROOT / "results/p2/X/full_r0" / name)} for name in ARTIFACTS]
    receipt_path.write_text(json.dumps({
        "schema_version": 1, "run": "X/full_r0", "source_gpu": run["gpu"], "files": import_files,
    }, indent=2) + "\n")
    sources = {f"results/p2/{run_name}/{name}": ROOT / f"results/p2/{run_name}/{name}"
               for run_name in ("tutorial", "X/full_r0") for name in ARTIFACTS}
    for pattern in ("*.py", "*.sh", "README.md"):
        for path in (ROOT / "scripts/slurm").glob(pattern):
            sources[path.relative_to(ROOT).as_posix()] = path
    receipt_member = "results/p2/slurm/imported_local_full0.json"
    sources[receipt_member] = receipt_path
    records = []
    for name, path in sorted(sources.items()):
        assert path.is_file() and not path.is_symlink() and path.stat().st_size
        assert not name.startswith("/") and ".." not in Path(name).parts
        records.append({"path": name, "size_bytes": path.stat().st_size, "sha256": checksum(path)})
    manifest = destination / "delta_overlay_manifest.json"
    manifest.write_text(json.dumps({"created_utc": identifier, "base_sha256": base["archive_sha256"],
                                   "imported_runs": ["tutorial", "X/full_r0"],
                                   "remaining_fit_tasks": list(range(1, 18)), "files": records}, indent=2) + "\n")
    checksums = destination / "OVERLAY_SHA256SUMS"
    checksums.write_text("".join(f"{item['sha256']}  {item['path']}\n" for item in records)
                         + f"{checksum(manifest)}  {manifest.name}\n", newline="\n")
    sources[manifest.name] = manifest
    sources[checksums.name] = checksums
    expected = {name: checksum(path) for name, path in sources.items()}
    archive = destination / f"regvelo-overlay-{identifier}.tar"
    def metadata(info):
        info.uid = info.gid = 0
        info.uname = info.gname = ""
        info.mode = 0o755 if info.name.endswith(".sh") else 0o644
        return info
    with tarfile.open(archive, "w", format=tarfile.PAX_FORMAT) as bundle:
        for name, path in sorted(sources.items()):
            bundle.add(path, arcname=name, recursive=False, filter=metadata)
    with tarfile.open(archive) as bundle:
        assert len(bundle.getmembers()) == len(expected) and set(bundle.getnames()) == set(expected)
        for member in bundle.getmembers():
            assert member.isfile()
            with bundle.extractfile(member) as stream:
                assert digest(stream) == expected[member.name], member.name
    plan = {
        "account": "bjvz-delta-gpu", "partition": "gpuH200x8", "remaining_fits": 17,
        "base_remote_archive": base["remote_archive"], "base_sha256": base["archive_sha256"],
        "overlay_local_archive": str(archive), "overlay_sha256": checksum(archive),
        "overlay_size_bytes": archive.stat().st_size,
        "overlay_remote_archive": f"/u/araghuram/{archive.name}",
        "remote_project": f"/projects/bjvz/araghuram/BioModelSimProject-{identifier}",
        "import_receipt_member": receipt_member, "import_receipt_sha256": checksum(receipt_path),
        "minimum_free_gib_before_extract": 6,
        "minimum_free_gib_before_fit": 2,
        "verified_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    (destination / "verification.json").write_text(json.dumps(plan, indent=2) + "\n")
    temporary = TRANSFER / "resume_plan.json.tmp"
    temporary.write_text(json.dumps(plan, indent=2) + "\n")
    temporary.replace(TRANSFER / "resume_plan.json")
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()
