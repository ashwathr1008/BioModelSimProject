"""Independently verify the saved main tables, AUCs and paired contrasts."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

def main():
    root = Path(__file__).resolve().parents[1]
    p4 = root / 'results/p4'
    manifest = pd.read_csv(p4 / 'run_manifest.csv')
    stats = pd.read_csv(p4 / 'shared_cell_effects_by_repeat.csv')
    rank = pd.read_csv(p4 / 'depletion_by_run.csv')
    anchors = json.loads((root / 'results/p3/terminal_anchors/terminal_anchors.json').read_text())
    assert len(manifest) == 18 and not manifest.duplicated(['condition','repeat','regulator']).any()
    total = 0
    joint_counts = {}
    max_delta_mass = 0
    for (rep,tf), block in manifest.groupby(['repeat','regulator']):
        wide = []
        for run in block.itertuples(index=False):
            f = root/run.path
            assert hashlib.sha256(f.read_bytes()).hexdigest() == run.sha256
            df = pd.read_csv(f)
            total += len(df)
            obs = df.pivot(index='cell_id', columns='endpoint', values='p_unperturbed')
            alt = df.pivot(index='cell_id', columns='endpoint', values='p_perturbed')
            ref = pd.read_csv(root / 'results/p3/fate' / f'{run.condition}_r{rep}_fate.csv',index_col='cell_id').loc[obs.index]
            np.testing.assert_allclose(ref[obs.columns], obs, rtol=0, atol=1e-12)
            assert alt.min().min() >= -1e-6 and alt.max().max() <= 1+1e-6
            np.testing.assert_allclose(alt.sum(axis=1),1,atol=1e-3,rtol=0)
            for ep,cells in anchors.items():
                np.testing.assert_allclose(alt.loc[cells,ep],1,atol=1e-3,rtol=0)
                v = roc_auc_score([1]*len(obs) + [0]*len(alt), np.r_[obs[ep],alt[ep]])
                expected = rank.loc[(rank.condition==run.condition) & (rank['repeat']==rep) & (rank.regulator==tf) & (rank.endpoint==ep),'depletion_auc'].item()
                assert abs(v-expected)<1e-12
            delta = alt - obs
            max_delta_mass = max(max_delta_mass, abs(delta.sum(axis=1)).max())
            wide.append((run.condition,delta.loc[~ref.cell_type.isin(anchors)]))
        matrices=dict(wide)
        common=sorted(set.intersection(*(set(t.index) for t in matrices.values())))
        joint_counts[int(rep)]=len(common)
        targ = (matrices['targeted'].loc[common] - matrices['full'].loc[common]).abs().mean()
        rand = (matrices['random'].loc[common] - matrices['full'].loc[common]).abs().mean()
        primary = stats.loc[(stats.source_group=='all_shared_nonterminal') & (stats['repeat']==rep) & (stats.regulator==tf)].set_index('endpoint')
        np.testing.assert_allclose(targ-rand, primary.loc[targ.index,'targeted_minus_random_mae'],atol=1e-12)
    report={'main_analyses':len(manifest),'cell_regulator_endpoint_rows':total,'shared_nonterminal_cells_by_repeat':joint_counts,'max_abs_sum_of_delta_across_endpoints':max_delta_mass,'validated':['File digests','Published unperturbed probabilities by cell ID','Perturbed fate bounds and normalization','106 anchor identities and fates','ROC AUC against rank-based depletion','Paired contrast re-derived from wide cell-by-endpoint matrices']}
    (p4/'qc/main_outputs.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report)


if __name__ == "__main__":
    main()
