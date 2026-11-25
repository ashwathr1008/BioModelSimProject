#!/usr/bin/env python3
"""Reproduce the verified repository-derived P5 summary figures.

Inputs: pre-aggregated tables verified directly against GitHub blobs on 2026-11-25.
These are not substitutes for original per-cell files; run_analysis.py handles those.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
SRC=HERE/'outputs'/'verified_from_github'
FIG=SRC/'figures'
FIG.mkdir(exist_ok=True,parents=True)

baseline=pd.read_csv(SRC/'baseline_fate_by_repeat.csv')
cbc=pd.read_csv(SRC/'cbc_matched_boundary_verified.csv')
p4=pd.read_csv(SRC/'p4_pooled_verified.csv')
stages=pd.read_csv(SRC/'baseline_fate_stage_bins_verified.csv')
celltypes=pd.read_csv(SRC/'baseline_fate_by_celltype.csv')
assert set(baseline.repeat)=={0,1,2} and len(baseline)==6
assert set(cbc.repeat)=={0,1,2} and len(cbc)==24
assert len(p4)==18 and len(stages)==12

# Cross-layer marginal *mean absolute error* comparison on the same P4 shared
# nonterminal cohort, with each regulator-endpoint focus represented separately.
# This does NOT infer per-cell correlation between the baseline and perturbation error.
endpoints=['PAX7+MyoProg1','MYH3+MyoC','MYL3+MyoC']
rows=[]
for _,r in p4.iterrows():
    b=baseline[(baseline.repeat==r['repeat'])&(baseline.condition=='targeted')]
    assert len(b)==1 and b.iloc[0].n_cells==r.n_cells
    for condition in ('targeted','random'):
        br=baseline[(baseline.repeat==r['repeat'])&(baseline.condition==condition)].iloc[0]
        err=float(br['mean_abs_'+r.endpoint])
        effect_err=float(r['mae_'+condition])
        rows.append(dict(repeat=int(r['repeat']),condition=condition,regulator=r.regulator,endpoint=r.endpoint,
                    n_cells=int(r.n_cells),mean_abs_baseline_fate_error=err,
                    mean_abs_perturbation_effect_error=effect_err,
                    effect_error_minus_baseline_error=effect_err-err,
                    provenance='P3 per-cell fate aggregation + P4 published pooled matched-cell summary'))
full=pd.DataFrame(rows)
full.to_csv(SRC/'cross_layer_marginal_errors_verified.csv',index=False)

# Paired CBC contrasts on the exact SAME evaluable boundary-cell subset in each repeat.
indexes=['repeat','source','target']
p=cbc.pivot(index=indexes,columns='condition',values=['mean_absolute_change','mean_signed_change','n_common_three'])
compact=[]
for idx in p.index:
    r,src,tgt=idx
    compact.append(dict(repeat=int(r),source=src,target=tgt,
        n_common_three=int(p.loc[idx,('n_common_three','targeted')]),
        targeted_abs_change=float(p.loc[idx,('mean_absolute_change','targeted')]),
        random_abs_change=float(p.loc[idx,('mean_absolute_change','random')]),
        targeted_minus_random=float(p.loc[idx,('mean_absolute_change','targeted')]-p.loc[idx,('mean_absolute_change','random')])))
pd.DataFrame(compact).to_csv(SRC/'cbc_matched_contrasts_verified.csv',index=False)

# Repeats summary primary
pvt=baseline.pivot(index='repeat',columns='condition',values='mean_tvd')
pvt['targeted_minus_random']=pvt.targeted-pvt.random
assert (pvt.targeted_minus_random>0).all()
pvt.to_csv(SRC/'baseline_fate_repeat_contrasts_verified.csv')

# Baseline & P4 contrast for headline comparisons (separate metrics).
fig,axes=plt.subplots(1,2,figsize=(11,4.5),sharex=True)
x=np.arange(3)
for ax,(tf,ep) in zip(axes,[('MSC','PAX7+MyoProg1'),('MYOG','MYH3+MyoC')]):
    sub=full[(full.regulator==tf)&(full.endpoint==ep)]
    for j,c in enumerate(('targeted','random')):
        d=sub[sub.condition==c].sort_values('repeat')
        ax.plot(x,d.mean_abs_baseline_fate_error.to_numpy(),marker='o',ls=('-' if c=='targeted' else '--'),label=f'{c}: baseline')
        ax.plot(x,d.mean_abs_perturbation_effect_error.to_numpy(),marker='s',ls=('-' if c=='targeted' else '--'),label=f'{c}: TF effect')
    ax.set_xticks(x,['R0','R1','R2']);ax.set_title(f'{tf} → {ep}')
    ax.set_ylabel('Mean absolute discrepancy vs Full (probability units)')
    ax.legend(fontsize=8)
fig.suptitle('Cross-layer sensitivity on shared nonterminal cells',fontsize=12)
fig.tight_layout()
fig.savefig(FIG/'verified_cross_layer_marginal_error.png',dpi=210)
plt.close(fig)

# CBC scatter: targeted vs random (four paths, three repeats)
fig,axs=plt.subplots(2,2,figsize=(10,8),sharex=True,sharey=True)
path_names=[('PAX3+MyoProg','MyoB1'),('MyoB1','MyoC1'),('PAX3+PAX7+MyoProg','MyoB2'),('MyoB2','MyoC2')]
for ax,(src,dst) in zip(axs.flat,path_names):
    q=cbc[(cbc.source==src)&(cbc.target==dst)]
    for c in ('targeted','random'):
        d=q[q.condition==c].sort_values('repeat')
        ax.plot(d.repeat,d.mean_absolute_change,marker='o',label=c)
    ax.set_title(f'{src} → {dst}',fontsize=10)
    ax.set_xticks([0,1,2],['R0','R1','R2']);ax.set_ylim(0,.16)
    ax.legend(fontsize=8)
    ax.set_ylabel('Mean abs paired CBC change vs Full')
fig.suptitle('P3: matched boundary-cell CBC sensitivity')
fig.tight_layout();fig.savefig(FIG/'verified_matched_cbc.png',dpi=220);plt.close(fig)

# Stage stratification
fig,axs=plt.subplots(1,2,figsize=(10,4),sharey=True)
for ax,stage in zip(axs,('early','late')):
    s=stages[stages.stage_bin==stage]
    for c in ('targeted','random'):
        z=s[s.condition==c].sort_values('repeat')
        ax.plot(z.repeat,z.mean_tvd,marker='o',label=c)
    ax.set_title('PCW ≤ 6.1' if stage=='early' else 'PCW > 6.1')
    ax.set_xticks([0,1,2],['R0','R1','R2']);ax.set_ylabel('Mean baseline TVD');ax.legend()
fig.suptitle('Baseline fate sensitivity by measured developmental stage')
fig.tight_layout();fig.savefig(FIG/'verified_stage_stratified_baseline.png',dpi=220);plt.close(fig)

# summary figure: all 3 kinds of contrasts (separate scales, linked by sign)
sub=pd.DataFrame(compact)
fig,axs=plt.subplots(1,3,figsize=(14,4.5))
axs[0].bar(x,pvt.targeted_minus_random);axs[0].axhline(0,lw=.8)
axs[0].set_xticks(x,['R0','R1','R2']);axs[0].set_title('Baseline fate TVD: T − R')
for ax,(src,dst) in zip(axs[1:],[('PAX3+MyoProg','MyoB1'),('MyoB1','MyoC1')]):
    q=sub[(sub.source==src)&(sub.target==dst)].sort_values('repeat')
    ax.bar(x,q.targeted_minus_random);ax.axhline(0,lw=.8)
    ax.set_xticks(x,['R0','R1','R2']);ax.set_title(f'Paired CBC: {src} → {dst}',fontsize=9)
fig.suptitle('Effect of targeted vs random thinning (positive = targeted more sensitive)')
fig.tight_layout();fig.savefig(FIG/'verified_cross_readout_contrasts.png',dpi=220);plt.close(fig)

summary=SRC/'verified_completion_report.md'
lines=['# P5: verified analytical update','',
       'Data provenance: `results/p3/fate/*.csv`, `results/p3/cell_level_cbc/all_cell_cbc.csv`, and `results/p4/shared_cell_effects_by_repeat.csv` fetched from the connected GitHub repository. Aggregated input tables in this folder were derived from complete text blobs (not synthetic examples).','',
       '## Main findings','',
       '**1. Unperturbed fate stability.** In all three computational blocks, targeted removal produced larger mean TVD than equal-sized random removal on the same nonterminal cells.','',
       '| Repeat | Targeted TVD | Random TVD | T−R |','|---|---:|---:|---:|']
for k,r in pvt.iterrows():lines.append(f'| R{k} | {r.targeted:.4f} | {r.random:.4f} | {r.targeted_minus_random:+.4f} |')
lines+=['','**2. Local trajectory CBC.** When CBC is compared on a common evaluable boundary-cell subset, both MyoB1-linked transitions have larger mean absolute changes under Targeted than Random in each of three repeats. This can be obscured by nearly unchanged unpaired CBC mean values. Boundary eligibility and size change substantially.','',
        '| Transition | R0 T/R | R1 T/R | R2 T/R |','|---|---:|---:|---:|']
for src,dst in path_names:
    q=sub[(sub.source==src)&(sub.target==dst)].sort_values('repeat')
    vals=[]
    for _,z in q.iterrows():vals.append(f'{z.targeted_abs_change:.3f} / {z.random_abs_change:.3f}')
    lines.append('| '+src+' → '+dst+' | '+' | '.join(vals)+' |')
lines+=['','**3. Stage stratification.** Within each repeat, targeted thinning causes larger mean baseline TVD than random thinning among early cells (PCW ≤6.1). This is descriptive; the stage effect is confounded with MyoB1 identity in the removal design.','',
         '**4. TF-effect stability.** P4 contrasts for the two headline axes have mixed signs over repeats. Cross-layer mean absolute discrepancies are combined in `cross_layer_marginal_errors_verified.csv` on the same shared nonterminal cohort. These marginal means cannot establish a *per-cell association* between baseline error and TF-effect error. For the per-cell joint relationship, the original 18 compressed P4 result tables are needed.','',
         '**5. 10×10 transition matrices.** Not computed because none of the nine `results/p2/X/*/fitted.h5ad` files exist in this runtime and the RegVelo/CellRank environment is unavailable. The necessary model archive is 1,429,028,483 bytes in the P2 GitHub release. The tested Python implementation in `run_analysis.py --matrix` will compute these after installation and extraction. Do not substitute CBC matrices or invented transition probabilities.','',
         '## Constraints','',
         '- All three repeats come from one biological dataset; they are not independent biological replicates. No population-level significance claims.','- Full is a computational reference, not ground truth.','- GRN prior and gene panel were fixed using original data information.','- Stage and cell identity are confounded under the present targeted vs random removal design.','- P4 manipulates TF outgoing regulatory weights in silico; it does not reproduce experimental knockout effects.']
summary.write_text('\n'.join(lines)+'\n')
print('Reports, tables, and figures prepared:',len(list(FIG.glob('*.png'))),'figures,',len(list(SRC.glob('*.csv'))),'CSV files')
print(summary)
