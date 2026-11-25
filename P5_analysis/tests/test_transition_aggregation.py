"""Independent of RegVelo and CellRank; validates the 10x10 reduction math."""
import sys
import unittest
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from transition_aggregation import aggregate_by_celltype

class KernelAggregationTest(unittest.TestCase):
    def test_all_and_common_destinations(self):
        K=csr_matrix(np.array([[.2,.3,.5],[0,.5,.5],[.6,.1,.3]]))
        labels=['A','A','B']
        # Cell0 and cell2 shared across all models; cell1 is a model-specific A destination.
        a,b,valid,counts=aggregate_by_celltype(K,labels,[0,2],[1,0,1],['A','B'])
        np.testing.assert_allclose(a,[[.5,.5],[.7,.3]])
        np.testing.assert_allclose(b,[[.2/.7,.5/.7],[.6/.9,.3/.9]])
        self.assertTrue(valid.all())
        np.testing.assert_equal(counts,[[2,1],[1,1]])
        np.testing.assert_allclose(a.sum(axis=1),1)
        np.testing.assert_allclose(b.sum(axis=1),1)

    def test_bad_kernel_detected(self):
        with self.assertRaises(ValueError):
            aggregate_by_celltype(np.eye(2)*.2,['A','B'],[0,1],[1,1],['A','B'])

if __name__=='__main__': unittest.main()
