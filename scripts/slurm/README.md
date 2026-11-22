# Historical Delta deployment workflow

The final selected X-only experiment completed locally on November 22. All Delta
jobs were cancelled before they ran. Use [the current fitting guide](../../docs/p2_model_fitting.md)
for the completed outputs and team handoff. The workflow below preserves the
earlier, expanded both-PCA deployment attempt; it is not required for the final
results and does not implement the selected X-only queue.

The Windows deployment and packaging helpers refer to Ashwath's original allocation,
project directories, and ignored session receipts. They are historical support code,
not portable launch commands. A future deployment must supply its own account,
paths, quota checks, and manifests, and must include `p2_run_all.py` along with
`p2_compare.py` and `p2_verify_environment.py`, which import its scope selection.

These scripts were prepared for 18 experimental fits (two PCA variants, three conditions,
three repeats). They exclude the separate tutorial. Each
array task requests one GPU and writes its own existing `results/p2/{pca}/{run}`
directory. The frozen scientific settings, full-batch training, 30 posterior
samples, and 2,048-cell ODE checkpointing are preserved.

## Before submission

On Delta, use `accounts`, `quota`, `sinfo`, and `scontrol show partition` to verify
your GPU allocation, remaining storage, partition access, and current scheduler
limits. The account is an allocation code ending in `-delta-gpu`, not a username.
Current official documentation lists `gpuH200x8` with H200 141 GB GPUs. Request one
GPU on a shared node; do not request an entire eight-GPU node for one fit. H200 has
a 3.0 charging factor, so verify allocation balance before submitting the array.
The defaults here are eight CPU cores, 32 GB host RAM, and a four-hour wall limit
per fit. They are starting resource requests, not measured Delta runtimes.

Find an available Python 3.12 interpreter with the module system or an existing
user installation. Do not assume a particular module name is installed. Run the
environment setup once in a location where installation is allowed by site
policy, before any GPU array task:

```bash
bash scripts/slurm/setup_environment.sh "$PWD" /verified/path/to/python3.12
```

The script installs Linux CUDA 12.4 PyTorch and the project requirements into
`.venv-delta`, verifies imports and package consistency, and records the Linux
environment. It uses an exclusive setup lock; fitting jobs hold shared locks so
installation cannot modify their environment. Array tasks never install packages.
The Windows `.venv` and Windows package freeze must not be used as a Linux environment.

The deployment bootstrap accepts the verified `/projects/bjvz/araghuram` location
through its canonical filesystem path, as well as the original work-HDD location.
It keeps Python, the environment, and package caches under `/u/araghuram`; the
default environment is `/u/araghuram/.local/share/regvelo/p2-py312-cu124`.
`results/p2/slurm/environment_path.txt` records the environment used so resumed
monitoring uses the same interpreter.

Fits run from a fresh node-local working directory under `${TMPDIR:-/tmp}` while
reading prepared inputs and writing collected model/results files at their
absolute project paths. The array checks for at least 2 GiB of filesystem free
space before each fit; `P2_MIN_FREE_GIB` changes that threshold. This filesystem
capacity check supplements the deployment's separate project-quota check.

## Submit the full campaign

After setup and live allocation checks, submit both stages in one authenticated
session. The first array runs the two full-data pilots; the second contains the
remaining 16 fits, with at most four running together:

```bash
export P2_ACCOUNT='VERIFIED_ALLOCATION-delta-gpu'
export P2_PARTITION='gpuH200x8'
P2_PRINT_ONLY=1 bash scripts/slurm/submit_campaign.sh
bash scripts/slurm/submit_campaign.sh
```

To reuse the completed laptop `X/full_r0`, transfer its four artifacts and an
independently hashed import receipt, then explicitly enable reuse:

```bash
export P2_REUSE_LOCAL_FULL0=1
export P2_IMPORT_COMPLETED='/absolute/path/to/imported_full0.json'
export P2_IMPORT_COMPLETED_SHA256='VERIFIED_RECEIPT_SHA256'
bash scripts/slurm/submit_campaign.sh
```

The receipt is JSON with `schema_version: 1`, `run: "X/full_r0"`, `source_gpu`
matching the run record, and `files` containing exactly four objects with `path`
and `sha256`: `results/p2/X/full_r0/run.json`, `training.csv`, `fitted.h5ad`, and
`model/model.pt` (each path includes the same run-directory prefix). Before any
submission, the controller verifies the pinned receipt, all four artifact hashes,
scientific input/settings identity, finite outputs, and RTX 4090 Laptop GPU
provenance. It freezes the receipt into the campaign directory and records the
validation in `import_validation.json`.

Reuse mode submits task `1` as the only GPU pilot, then tasks `2-17`; task `0`
is never submitted. Only the explicitly authorized `X/full_r0` may have laptop
GPU provenance. Every newly allocated GPU and all other experimental fit records
must match the expected H200. The final comparison marker records the import
receipt hash and `mixed_hardware: true`. Include the completed tutorial's same
four artifacts as well to enable the full three-baseline comparison.

The second array has an `afterok` dependency on the entire pilot array. A failed
pilot cancels dependent tasks through `--kill-on-invalid-dep=yes`. Each remaining
task also checks both pilot records, input/settings identities, output files,
finite training losses and cell times/stages, and the H200 GPU name before fitting.
This check does not select models by the size or sign of a biological effect.
The allocated GPU name is checked before every fit. `P2_EXPECT_GPU` can explicitly
change the expected GPU when using a different verified partition.

Campaign receipts under `results/p2/slurm/campaign.*` contain `pilots.job_id`,
`remainder.job_id`, `jobs.tsv`, and a frozen configuration. Each array also has its
own submission receipt. If submission stops partway, inspect those receipts
before resubmitting; an earlier array may already have been accepted by Slurm.
`P2_CAMPAIGN_CONCURRENCY` overrides the default four concurrent remaining fits.

After each remaining fit, a shared file lock protects a completion check. The
last task to observe all 18 fits complete validates the artifacts and runs
`p2_compare.py` once within its existing allocation. It writes
`results/p2/slurm/comparison_complete.json` only after comparison succeeds, with
the configuration SHA-256 and checked runs. A comparison failure fails that
array task and leaves the success marker absent. Missing tutorial outputs remain
reported as unavailable until the laptop tutorial is transferred or comparison
is rerun locally with all outputs.

## Manual submission and inspection

Set the verified GPU account and partition, then submit the two full-data pilots:

```bash
export P2_ACCOUNT='VERIFIED_ALLOCATION-delta-gpu'
export P2_PARTITION='gpuH200x8'
P2_PRINT_ONLY=1 P2_ARRAY='0-1%2' bash scripts/slurm/submit_array.sh
P2_ARRAY='0-1%2' bash scripts/slurm/submit_array.sh
```

`P2_PRINT_ONLY=1` writes a local submission snapshot and prints the command without
calling `sbatch`. Every submission copies the configuration into its receipt
directory and records its SHA-256. Each job checks the snapshot before training.
The allocated GPU is checked on the compute node, not the login node.

After both pilots complete and their outputs are checked, submit tasks `2-17`.
Choose the `%N` concurrency limit from the live allocation/QOS limits and desired
allocation use. For example, `%2` runs at most two fits at once:

```bash
P2_ARRAY='2-17%2' bash scripts/slurm/submit_array.sh
squeue -u "$USER"
```

`P2_CPUS`, `P2_MEMORY`, `P2_TIME`, `P2_ENV`, `PROJECT_DIR`, and `P2_SOURCE_CONFIG`
can override the documented defaults. Use absolute project/environment paths.
`P2_DEPENDENCY` passes an explicit dependency such as `afterok:12345` to `sbatch`.
`P2_REQUIRE_PILOTS=1` enables pilot validation, and
`P2_COMPARE_WHEN_COMPLETE=1` enables the guarded final comparison for manual arrays.
Changing scheduler resources does not change model settings. Logs are
`results/p2/slurm/<array-job>_<task>.out` and `.err`; submission receipts and config
snapshots are in `results/p2/slurm/submission.*`.

Task mapping repeats the same six entries per repeat:

| Task offsets | PCA | Condition |
|---|---|---|
| 0, 6, 12 | X | full |
| 1, 7, 13 | shipped_pca | full |
| 2, 8, 14 | X | targeted |
| 3, 9, 15 | shipped_pca | targeted |
| 4, 10, 16 | X | random |
| 5, 11, 17 | shipped_pca | random |

Tasks 0–5 use repeat 0, 6–11 repeat 1, and 12–17 repeat 2. Existing completed fits
are skipped only if the fitting script's input/setting identity checks pass.
Failed fits restart from the beginning; the script does not resume within an epoch.

Run comparisons once after a completed batch, never simultaneously from all array
tasks, because the comparison script writes shared summary files:

```bash
.venv-delta/bin/python scripts/p2_compare.py --config config/p2_both_pca.yaml
```

## Transfer inventory

Preserve this repository-relative structure on Delta:

- Code: both requirements files, `config/experiment.yaml`,
  `config/p2_both_pca.yaml`, `src/p1_data/`, `scripts/p2_fit.py`,
  `scripts/p2_memory.py`, `scripts/p2_compare.py`, and `scripts/slurm/`.
- Inputs: the nine H5ADs in each of `data/prepared_p2/` and
  `data/prepared_p2_shipped_pca/`, plus all nine manifests in
  `results/p2/preprocessing/manifests/`.
- Provenance: `results/p2/input_validation.json` and
  `results/p2/memory_validation.json`.
- Optional: the original raw H5AD and `scripts/p2_verify_environment.py` for repeating
  all input/GPU checks remotely; the completed laptop tutorial's checkpoint,
  `fitted.h5ad`, `training.csv`, and `run.json` for the three-baseline comparison.

Do not transfer `.venv`, `.venv-py310-original`, local queue/keep-awake state, or
credentials. Hash the selected source, manifests, and prepared files before and
after transfer. Raw data are not required for the 18 prepared-input fits.

Official references: [GPU architecture](https://docs.ncsa.illinois.edu/systems/delta/en/latest/user_guide/architecture.html),
[running jobs](https://docs.ncsa.illinois.edu/systems/delta/en/latest/user_guide/running_jobs.html),
[job accounting](https://docs.ncsa.illinois.edu/systems/delta/en/latest/user_guide/job_accounting.html).
