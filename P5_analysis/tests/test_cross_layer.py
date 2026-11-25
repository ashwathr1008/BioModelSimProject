"""Deterministic synthetic P4 fixture for cross-layer joining and per-cell association."""
import tempfile
import unittest
from pathlib import Path
import sys
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from run_analysis import analyze_cross_layer, ENDPOINTS

class CrossLayerFixtureTest(unittest.TestCase):
    def test_join_associations_and_output(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name)
            out=root/'P5_analysis'/'outputs'
            (out/'tables').mkdir(parents=True)
            (out/'figures').mkdir(parents=True)
            base={"c1":[.6,.3,.1],"c2":[.3,.6,.1],"c3":[.4,.2,.4],"c4":[.2,.3,.5]}
            ct={"c1":"MyoB1","c2":"MyoB1","c3":"MyoB2","c4":"MyoB2"}
            for repeat in range(3):
                for regulator in ['MSC','MYOG']:
                    for c in ['full','targeted','random']:
                        records=[]
                        offset={'full':0,'targeted':.02,'random':.01}[c]
                        for idx,(cid,prob) in enumerate(base.items()):
                            # Ensure perturbation probabilities remain in [0,1].
                            baseline=np.array(prob)+np.array([offset,-offset,0])*(idx+1)
                            effect=(idx+1)*{'full':.012,'targeted':.021,'random':.016}[c]
                            if regulator=='MYOG':effect*=.65
                            new=baseline+np.array([effect,-effect,0])
                            for j,endpoint in enumerate(ENDPOINTS):
                                records.append(dict(cell_id=cid,cell_type=ct[cid],condition=c,
                                    repeat=repeat,regulator=regulator,endpoint=endpoint,
                                    p_unperturbed=baseline[j],p_perturbed=new[j],delta_p=new[j]-baseline[j]))
                        path=root/'results'/'p4'/f'{c}_r{repeat}'/f'{regulator}.csv.gz'
                        path.parent.mkdir(parents=True,exist_ok=True)
                        pd.DataFrame(records).to_csv(path,index=False,compression='gzip')
            result=analyze_cross_layer(root,out)
            self.assertEqual(len(result),3*2*3*2)
            association=pd.read_csv(out/'tables'/'cross_layer_per_cell_association.csv')
            self.assertEqual(len(association),36)
            self.assertTrue(association.spearman_rho.dropna().between(-1,1).all())
            self.assertTrue((association.top_decile_overlap.between(0,1)).all())
            self.assertEqual(result.n_cells.unique().tolist(),[4])
            self.assertTrue((out/'figures'/'fig_cross_layer_error.png').is_file())

if __name__=='__main__':unittest.main()
