"""Full seed reference must use each repeat's *three-way* evaluation cohort."""
import sys
import tempfile
import unittest
from pathlib import Path
import pandas as pd
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from run_analysis import full_seed_variability, ENDPOINTS

class SeedReferenceTest(unittest.TestCase):
    def test_seed_variability_cohort_alignment(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name)
            out=root/'P5_analysis'/'outputs'
            (out/'tables').mkdir(parents=True)
            inp=root/'results/p3/fate';inp.mkdir(parents=True)
            for r in range(3):
                for c in ['full','targeted','random']:
                    rows=[]
                    for cid,typ in [('c1','MyoB1'),('c2','MyoB2'),('c3','MyoB1'),('c4','PAX7+MyoProg1')]:
                        if c=='targeted' and cid=='c3':continue
                        if c=='random' and cid=='c3':continue
                        p=[.6,.3,.1]
                        if c=='full':p=[.6-.05*r,.3+.05*r,.1]
                        rows.append(dict(cell_id=cid,**dict(zip(ENDPOINTS,p)),cell_type=typ,condition=c,repeat=r,adj_stage_num=6.1))
                    pd.DataFrame(rows).to_csv(inp/f'{c}_r{r}_fate.csv',index=False)
            summary=full_seed_variability(root,out)
            self.assertEqual(len(summary),27)
            one=summary[(summary.evaluation_repeat==0)&(summary.reference_comparison=='full_r0_vs_full_r1')&(summary.source_group=='all_shared_nonterminal')].iloc[0]
            self.assertEqual(one.n_cells,2)
            self.assertAlmostEqual(one.mean_tvd,.05)

if __name__=='__main__':unittest.main()
