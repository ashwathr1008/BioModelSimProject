"""Create and verify an immutable, credential-free Delta deployment bundle."""

import hashlib
import json
import tarfile
import uuid
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TRANSFER = ROOT / "results/p2/delta_transfer"


def digest_stream(stream):
    digest = hashlib.sha256()
    for chunk in iter(lambda: stream.read(2**20), b""):
        digest.update(chunk)
    return digest.hexdigest()


def digest_file(path):
    with path.open("rb") as stream:
        return digest_stream(stream)


def main():
    seed = json.loads((TRANSFER / "20261121T180806Z_69f697b8/delta_transfer_manifest.json").read_text())
    members = {entry["path"] for entry in seed["files"]}
    for pattern in ("*.sh", "*.py", "README.md"):
        members.update(path.relative_to(ROOT).as_posix()
                       for path in (ROOT / "scripts/slurm").glob(pattern))
    identifier = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    destination = TRANSFER / identifier
    destination.mkdir(parents=True, exist_ok=False)
    records = []
    for name in sorted(members):
        path = ROOT / name
        if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(ROOT):
            raise ValueError(f"Unsafe or missing source: {name}")
        records.append({"path": name, "size_bytes": path.stat().st_size, "sha256": digest_file(path)})
    prepared = [record for record in records if record["path"].endswith(".h5ad")]
    assert len(prepared) == 18
    manifest_path = destination / "delta_transfer_manifest.json"
    manifest_path.write_text(json.dumps({
        "created_utc": identifier, "purpose": "18 RegVelo fits on Delta H200; both PCA variants",
        "payload_file_count": len(records), "prepared_input_count": len(prepared),
        "files": records,
        "excluded": ["credentials", "SSH files", "virtual environments", ".git",
                     "raw input", "local fit outputs", "local queue state"],
    }, indent=2) + "\n")
    checksums_path = destination / "SHA256SUMS"
    checksums_path.write_text("".join(f"{record['sha256']}  {record['path']}\n" for record in records)
                              + f"{digest_file(manifest_path)}  delta_transfer_manifest.json\n", newline="\n")
    archive = destination / f"regvelo-p2-{identifier}.tar"
    def metadata(info):
        info.uid = info.gid = 0
        info.uname = info.gname = ""
        info.mode = 0o755 if info.name.endswith(".sh") else 0o644
        return info
    with tarfile.open(archive, "w", format=tarfile.PAX_FORMAT) as bundle:
        for record in records:
            bundle.add(ROOT / record["path"], arcname=record["path"], recursive=False, filter=metadata)
        for path in (manifest_path, checksums_path):
            bundle.add(path, arcname=path.name, recursive=False, filter=metadata)
    expected = {record["path"]: record["sha256"] for record in records}
    expected.update({path.name: digest_file(path) for path in (manifest_path, checksums_path)})
    with tarfile.open(archive) as bundle:
        assert len(bundle.getmembers()) == len(expected)
        assert set(bundle.getnames()) == set(expected)
        for member in bundle.getmembers():
            assert member.isfile() and not member.name.startswith("/") and ".." not in Path(member.name).parts
            with bundle.extractfile(member) as stream:
                assert digest_stream(stream) == expected[member.name], member.name
    plan = {
        "account": "bjvz-delta-gpu", "partition": "gpuH200x8",
        "archive_path": str(archive), "archive_sha256": digest_file(archive),
        "archive_size_bytes": archive.stat().st_size, "verified_members": len(expected),
        "remote_archive": f"/u/araghuram/{archive.name}",
        "remote_project": f"/work/hdd/bjvz/araghuram/BioModelSimProject-{identifier}",
        "prepared_inputs": 18, "pilot_fits": 2, "remaining_fits": 16,
        "remaining_concurrency": 4, "gpus_per_fit": 1, "walltime_per_fit": "04:00:00",
        "verified_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    (destination / "verification.json").write_text(json.dumps(plan, indent=2) + "\n")
    (destination / (archive.name + ".sha256")).write_text(f"{plan['archive_sha256']}  {archive.name}\n")
    temporary = TRANSFER / "deployment.json.tmp"
    temporary.write_text(json.dumps(plan, indent=2) + "\n")
    temporary.replace(TRANSFER / "deployment.json")
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()
