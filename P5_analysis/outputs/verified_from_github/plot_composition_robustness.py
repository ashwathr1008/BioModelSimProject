from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parent
x=pd.read_csv(p/'baseline_fate_by_celltype.csv')
records=[]
for (r,c),d in x.groupby(['repeat','condition'],sort=True):
    for grp,part in [('all_nonterminal_types',d),('excluding_MyoB1_sources',d[d.cell_type!='MyoB1'])]:
        records.append(dict(repeat=r,condition=c,source_cohort=grp,n_source_types=len(part),n_cells=int(part.n_cells.sum()),
                            cell_weighted_mean_tvd=np.average(part.mean_tvd,weights=part.n_cells),equal_celltype_mean_tvd=part.mean_tvd.mean()))
df=pd.DataFrame(records)
df.to_csv(p/'fate_composition_robustness_verified.csv',index=False)
fig,axs=plt.subplots(1,2,figsize=(11,4.5),sharey=True)
for ax,(name,grp) in zip(axs,[('All seven nonterminal source types','all_nonterminal_types'),('Excluding remaining MyoB1 sources','excluding_MyoB1_sources')]):
    z=df[df.source_cohort==grp]
    for cond in ['targeted','random']:
        sub=z[z.condition==cond].sort_values('repeat')
        ax.plot(sub.repeat,sub.cell_weighted_mean_tvd,marker='o',label=cond)
    ax.set_xticks([0,1,2],['R0','R1','R2']);ax.set_title(name);ax.set_ylabel('Mean TVD vs Full')
    ax.legend()
fig.suptitle('Sampling-effect robustness: is the TVD contrast limited to residual MyoB1 cells?')
fig.tight_layout();fig.savefig(p/'figures'/'verified_excluding_myoB1.png',dpi=220);plt.close(fig)
print(df.to_string(index=False))
