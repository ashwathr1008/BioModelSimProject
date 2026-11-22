#!/usr/bin/env bash
# Run once before submission, using a verified Linux Python 3.12 interpreter.
set -euo pipefail

PROJECT_DIR="${1:?Usage: bash setup_environment.sh /absolute/project /python3.12 [/absolute/venv]}"
PYTHON_BIN="${2:?Supply the Python 3.12 interpreter available on Delta}"
P2_ENV="${3:-$PROJECT_DIR/.venv-delta}"
PROJECT_DIR="$(cd -- "$PROJECT_DIR" && pwd -P)"
[[ "$P2_ENV" = /* ]] || { echo 'The environment path must be absolute' >&2; exit 2; }
[[ "$(uname -s)" = Linux && "$(uname -m)" = x86_64 ]] || {
    echo 'This environment targets Linux x86_64 NVIDIA nodes.' >&2; exit 2;
}
"$PYTHON_BIN" -c 'import sys; assert sys.version_info[:2] == (3, 12), sys.version'
mkdir -p -- "$(dirname -- "$P2_ENV")" "$PROJECT_DIR/results/p2/slurm"
exec 9>"${P2_ENV}.setup.lock"
flock --exclusive 9
if [[ ! -x "$P2_ENV/bin/python" ]]; then
    "$PYTHON_BIN" -m venv "$P2_ENV"
fi
PYTHON="$P2_ENV/bin/python"
"$PYTHON" -c 'import sys; assert sys.version_info[:2] == (3, 12) and sys.prefix != sys.base_prefix'
rm -f -- "$P2_ENV/.p2_environment_ready.json"
"$PYTHON" -m pip install --upgrade pip
"$PYTHON" -m pip install --no-compile 'torch==2.5.1' --index-url https://download.pytorch.org/whl/cu124
"$PYTHON" -m pip install --no-compile -r "$PROJECT_DIR/requirements-p2.txt"
"$PYTHON" -m pip check
"$PYTHON" -m pip freeze --all > "$PROJECT_DIR/results/p2/slurm/environment-lock-linux.txt"
"$PYTHON" - "$PROJECT_DIR" "$P2_ENV" <<'PY'
import hashlib, json, platform, sys
from importlib.metadata import version
from pathlib import Path
import cellrank, regvelo, scvelo, scvi, torch
project, environment = map(Path, sys.argv[1:])
assert torch.version.cuda == "12.4", torch.version.cuda
record = {
    "python": platform.python_version(), "platform": platform.platform(),
    "interpreter": sys.executable, "cuda_build": torch.version.cuda,
    "requirements_sha256": {
        name: hashlib.sha256((project / name).read_bytes()).hexdigest()
        for name in ("requirements-p1.txt", "requirements-p2.txt")
    },
    "packages": {name: version(name) for name in
                 ("regvelo", "scvi-tools", "torch", "torchode", "lightning", "numpy", "scipy")},
    "gpu_test": "pending an allocated GPU node",
}
(environment / ".p2_environment_ready.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2))
PY
echo 'Environment installation/import checks passed; run a GPU pilot before the remaining array.'
