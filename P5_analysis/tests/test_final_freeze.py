from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from verify_final_package import verify


def test_final_freeze_integrity():
    assert not verify()


def test_canonical_final_outputs_present():
    paths = (
        "PROJECT_FINAL_GUIDE.md",
        "outputs/reports/findings.md",
        "outputs/tables/baseline_fate_by_repeat.csv",
        "outputs/tables/cbc_matched_boundary_cells.csv",
        "outputs/tables/perturbation_effects_existing_p4.csv",
        "outputs/tables/cross_layer_per_cell_association.csv",
        "outputs/figures/fig_baseline_fate_tvd_by_repeat.png",
    )
    for rel in paths:
        assert (ROOT / rel).is_file(), rel
