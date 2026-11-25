"""Regression tests for published headline numbers in PROJECT_FINAL_GUIDE.md.

These tests use the frozen output CSVs, NOT re-trained models.
"""
from pathlib import Path

import pandas as pd

P = Path(__file__).resolve().parents[1] / 'outputs' / 'tables'


def test_primary_fate_mean_contrast_all_repeats():
    d = pd.read_csv(P / 'baseline_fate_contrast_by_repeat.csv').sort_values('repeat')
    assert d.repeat.tolist() == [0, 1, 2]
    assert d.n_cells.tolist() == [8234, 8235, 8243]
    assert (d.mean_targeted_minus_random > 0).all()
    assert abs(d.iloc[0].mean_targeted_minus_random - 0.017848) < 1e-5
    # The mean and median need not agree; specifically R0 median is < 0.
    assert d.iloc[0].median_targeted_minus_random < 0


def test_myoB1_matched_cbc_and_coverage_claims():
    c = pd.read_csv(P / 'cbc_matched_boundary_cells.csv')
    q = c[(c.source == 'MyoB1') & (c.target == 'MyoC1')]
    assert q.shape[0] == 6
    for r, n in enumerate([31, 30, 22]):
        t = q[(q.repeat == r) & (q.condition == 'targeted')].iloc[0]
        random = q[(q.repeat == r) & (q.condition == 'random')].iloc[0]
        assert t.n_paired_boundary_cells == n == random.n_paired_boundary_cells
        assert t.mean_abs_cbc_change > random.mean_abs_cbc_change
    cover = pd.read_csv(P / 'cbc_coverage_and_direction.csv')
    z = cover[(cover.source == 'MyoB1') & (cover.target == 'MyoC1')]
    assert z[z.condition == 'full'].sort_values('repeat').n_valid.tolist() == [120, 120, 120]
    assert z[z.condition == 'targeted'].sort_values('repeat').n_valid.tolist() == [44, 48, 37]


def test_p4_main_response_contrasts_mixed_signs():
    d = pd.read_csv(P / 'perturbation_effects_existing_p4.csv')
    for tf, endpoint in [('MSC', 'PAX7+MyoProg1'), ('MYOG', 'MYH3+MyoC')]:
        x = d[(d.regulator == tf) & (d.endpoint == endpoint)].sort_values('repeat')
        assert len(x) == 3
        vals = x.targeted_minus_random_mae.to_numpy()
        assert (vals > 0).any() and (vals < 0).any()


def test_genuine_per_cell_cross_layer_output_complete():
    d = pd.read_csv(P / 'cross_layer_per_cell_association.csv')
    assert len(d) == 36
    assert set(d.regulator) == {'MSC', 'MYOG'}
    assert set(d.condition) == {'targeted', 'random'}
    assert set(d.repeat) == {0, 1, 2}
    assert d.n_cells.isin([8234, 8235, 8243]).all()
    assert d.spearman_rho.notna().all()
