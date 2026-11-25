"""Small deterministic fixture test: use `python -m unittest discover -s tests` from P5_analysis."""
import sys
import tempfile
import unittest
from pathlib import Path
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from run_analysis import paired_baseline_fate

class BaselineFateTest(unittest.TestCase):
    def test_shared_nonterminal_and_total_variation(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);(root/'results/p3/fate').mkdir(parents=True)
            out=root/'P5_analysis/outputs';(out/'tables').mkdir(parents=True);(out/'figures').mkdir()
            # 2 nonterminal cells + 1 protected endpoint; Random removes one different cell.
            for r in range(3):
                base={
                    'c1':([.4,.5,.1],'MyoB1',5.6),
                    'c2':([.2,.7,.1],'MyoB2',6.1),
                    'c3':([1.,0.,0.],'PAX7+MyoProg1',6.5),
                    'c4':([.3,.6,.1],'MyoC1',7.2),
                }
                for condition in ('full','targeted','random'):
                    rows=[]
                    for id,(p,typ,stage) in base.items():
                        if condition=='targeted' and id=='c4':continue
                        if condition=='random' and id=='c4':continue
                        q=p.copy()
                        if condition=='targeted' and id=='c1':q=[.6,.3,.1]
                        if condition=='random' and id=='c1':q=[.45,.45,.1]
                        rows.append([id,*q,typ,condition,r,0.3,stage])
                    pd.DataFrame(rows,columns=['cell_id','PAX7+MyoProg1','MYH3+MyoC','MYL3+MyoC','cell_type','condition','repeat','latent_time','adj_stage_num']).to_csv(root/f'results/p3/fate/{condition}_r{r}_fate.csv',index=False)
            x=paired_baseline_fate(root,out)
            self.assertEqual(len(x),6)
            self.assertEqual(x.iloc[0].n_cells,2)
            self.assertAlmostEqual(x[(x.repeat==0)&(x.condition=='targeted')].mean_tvd.iloc[0],.1)
            self.assertAlmostEqual(x[(x.repeat==0)&(x.condition=='random')].mean_tvd.iloc[0],.025)
            self.assertTrue((out/'tables/baseline_fate_contrast_by_repeat.csv').is_file())

if __name__=='__main__': unittest.main()
