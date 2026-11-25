from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
here=Path(__file__).resolve().parent
b=pd.read_csv(here/'baseline_fate_by_repeat.csv')
s=pd.read_csv(here/'full_training_seed_variability_verified.csv')
fig,axes=plt.subplots(1,2,figsize=(12,4.5),sharex=True)
for ax,group in zip(axes,['all_shared_nonterminal','MyoB1']):
    if group=='all_shared_nonterminal':f=b
    else:
        f=pd.read_csv(here/'baseline_fate_by_celltype.csv').query('cell_type == "MyoB1"')
    for c in ['targeted','random']:
        q=f[f.condition==c].sort_values('repeat')
        ax.plot(q.repeat,q.mean_tvd,marker='o',label=c)
    seeds=s[s.source_group==group]
    for comparison in seeds.reference_comparison.unique():
        x=seeds[seeds.reference_comparison==comparison].sort_values('evaluation_repeat')
        ax.plot(x.evaluation_repeat,x.mean_tvd,ls=':',lw=1.7,marker='x',alpha=.75,label=comparison)
    ax.set_title('All shared nonterminal cells' if group=='all_shared_nonterminal' else 'Remaining MyoB1 cells')
    ax.set_xticks([0,1,2],['R0','R1','R2'])
    ax.set_ylabel('Mean fate TVD')
    ax.legend(fontsize=8)
fig.suptitle('Sampling contrast versus cross-seed Full-model variation\n(each panel uses the identical evaluation cohort within a repeat)', fontsize=11)
fig.tight_layout()
out=here/'figures'/'verified_full_seed_variability.png'
fig.savefig(out,dpi=220)
plt.close(fig)
print(out)
