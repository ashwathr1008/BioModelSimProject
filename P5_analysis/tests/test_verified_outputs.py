"""Regression tests for source-derived P5 aggregate tables."""
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

D=Path(__file__).resolve().parents[1]/'outputs'/'verified_from_github'

class VerifiedRegressionTests(unittest.TestCase):
    def test_baseline_tvd_three_repeats(self):
        s=pd.read_csv(D/'baseline_fate_by_repeat.csv')
        p=s.pivot(index='repeat',columns='condition',values='mean_tvd')
        self.assertEqual(p.shape,(3,2))
        self.assertTrue((p.targeted>p.random).all())
        np.testing.assert_allclose((p.targeted-p.random).to_numpy(),[.01784794,.02237048,.00758477],atol=2e-7)

    def test_cbc_common_boundary_counts(self):
        a=pd.read_csv(D/'cbc_matched_boundary_verified.csv')
        p=a[(a.source=='MyoB1')&(a.target=='MyoC1')].pivot(index='repeat',columns='condition',values='mean_absolute_change')
        self.assertTrue((p.targeted>p.random).all())
        self.assertEqual(a[(a.source=='MyoB1')&(a.target=='MyoC1')&(a.condition=='targeted')].sort_values('repeat').n_common_three.tolist(),[31,30,22])

    def test_p4_contrasts_have_mixed_signs(self):
        a=pd.read_csv(D/'p4_pooled_verified.csv')
        self.assertEqual(len(a),18)
        for tf,endpoint in [('MSC','PAX7+MyoProg1'),('MYOG','MYH3+MyoC')]:
            x=a[(a.regulator==tf)&(a.endpoint==endpoint)].targeted_minus_random_mae.to_numpy()
            self.assertGreater(x.max(),0)
            self.assertLess(x.min(),0)

    def test_cross_layer_aggregate_matches_p4(self):
        x=pd.read_csv(D/'cross_layer_marginal_errors_verified.csv')
        self.assertEqual(len(x),36)
        self.assertTrue(np.isfinite(x[['mean_abs_baseline_fate_error','mean_abs_perturbation_effect_error']].to_numpy()).all())
        self.assertEqual(x.n_cells.min(),8234)

if __name__=='__main__': unittest.main()
