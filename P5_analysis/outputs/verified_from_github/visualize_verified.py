"""Plot static summary tables computed from full connected-GitHub CSVs."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
charts=ROOT/'figures'
charts.mkdir(exist_ok=True)
base=pd.read_csv(ROOT/'baseline_fate_by_repeat.csv')
cell=pd.read_csv(ROOT/'baseline_fate_by_celltype.csv')
p4=pd.read_csv(ROOT/'p4_key_contrasts.csv')
a=base.pivot(index='repeat',columns='condition',values='mean_tvd')
fig,ax=plt.subplots(figsize=(7.8,4.6));xx=np.arange(3)
ax.bar(xx-.18,a.targeted.to_numpy(),width=.36,label='Targeted')
ax.bar(xx+.18,a.random.to_numpy(),width=.36,label='Random')
for i,(t,r) in enumerate(zip(a.targeted,a.random)):
    ax.text(i,t+.0015,f'+{(t-r):.3f}',ha='center',fontsize=9)
ax.set_xticks(xx,['R0','R1','R2']);ax.set_ylabel('Mean TVD vs Full')
ax.set_title('Baseline fate: Targeted > Random in all three repeats')
ax.legend();fig.tight_layout();fig.savefig(charts/'verified_baseline_fate_TVD.png',dpi=220);plt.close(fig)

h=cell.pivot(index=['repeat','cell_type'],columns='condition',values='mean_tvd')
d=(h.targeted-h.random).unstack('repeat').sort_index();arr=d.to_numpy();lim=np.max(np.abs(arr))
fig,ax=plt.subplots(figsize=(8.5,5));im=ax.imshow(arr,cmap='coolwarm',vmin=-lim,vmax=lim,aspect='auto')
ax.set_xticks(range(3),['R0','R1','R2']);ax.set_yticks(range(len(d.index)),d.index)
for i in range(len(d.index)):
 for j in range(3):ax.text(j,i,f'{arr[i,j]:+.3f}',ha='center',va='center',fontsize=8,color='black')
ax.set_title('Cell-type sensitivity: Targeted TVD minus Random TVD')
fig.colorbar(im,ax=ax,label='Mean TVD contrast');fig.tight_layout()
fig.savefig(charts/'verified_fate_celltype_contrast.png',dpi=220);plt.close(fig)

fig,axes=plt.subplots(2,1,figsize=(8,6.3),sharex=True)
axes[0].bar(xx,(a.targeted-a.random).to_numpy(),width=.5)
axes[0].axhline(0,color='black',linewidth=.8);axes[0].set_ylabel('Mean fate TVD contrast')
axes[0].set_title('Baseline fate: Targeted minus Random')
for tf,ep in [('MSC','PAX7+MyoProg1'),('MYOG','MYH3+MyoC')]:
 q=p4[(p4.regulator==tf)&(p4.endpoint==ep)].sort_values('repeat')
 axes[1].plot(xx,100*q.targeted_minus_random_mae,marker='o',label=f'{tf} -> {ep}')
axes[1].axhline(0,color='black',linewidth=.8);axes[1].set_ylabel('TF effect discrepancy contrast (pp)')
axes[1].set_title('P4 perturbation: no consistent targeted excess')
axes[1].legend(fontsize=8);axes[1].set_xticks(xx,['R0','R1','R2'])
fig.tight_layout();fig.savefig(charts/'verified_baseline_vs_P4_contrast.png',dpi=220);plt.close(fig)
