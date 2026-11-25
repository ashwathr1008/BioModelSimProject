#!/usr/bin/env python3
"""RegVelo P5: paired cell-level fate, CBC, TF effects and optional transition analysis.

Run from anywhere: python P5_analysis/run_analysis.py [--project-root PATH] [--matrix]
Inputs are read-only. All generated files go into P5_analysis/outputs/.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ENDPOINTS = ('PAX7+MyoProg1', 'MYH3+MyoC', 'MYL3+MyoC')
CONDS = ('full', 'targeted', 'random')
REPEATS = range(3)
LOG = logging.getLogger('P5')

def locate_project(explicit: str | None, script: Path) -> Path:
    if explicit:
        p = Path(explicit).expanduser().resolve()
        if not (p/'config'/'experiment.yaml').is_file():
            raise FileNotFoundError(f'No config/experiment.yaml in {p}. Unzip the GitHub repository first.')
        return p
    candidates = [script.parent, script.parent.parent, *script.parents, *script.parent.parent.glob('*')]
    for p in candidates:
        if (p/'config'/'experiment.yaml').is_file():
            return p.resolve()
    raise FileNotFoundError('Could not locate extracted repository. Use --project-root /path/to/BioModelSimProject')

def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f'Required input missing: {path}. The unextracted ZIP cannot be read directly.')
    return pd.read_csv(path, **kwargs)

def read_fate(root: Path, condition: str, repeat: int) -> pd.DataFrame:
    p = root/'results'/'p3'/'fate'/f'{condition}_r{repeat}_fate.csv'
    d = read_csv(p, index_col='cell_id')
    if not d.index.is_unique:
        raise ValueError(f'Duplicate cell IDs in {p}')
    if not set(ENDPOINTS).issubset(d.columns) or not {'cell_type','condition','repeat','adj_stage_num'}.issubset(d.columns):
        raise ValueError(f'Missing fate columns in {p}')
    if not ((d.condition == condition).all() and (d.repeat.astype(int) == repeat).all()):
        raise ValueError(f'Incorrect condition or repeat in {p}')
    a = d.loc[:,list(ENDPOINTS)].to_numpy(dtype=float)
    if not np.isfinite(a).all() or (a < -1e-5).any() or (a > 1.00001).any() or not np.allclose(a.sum(axis=1), 1, atol=1e-3):
        raise ValueError(f'Invalid fate probabilities in {p}')
    return d

def aggregate_fate(d: pd.DataFrame, repeat: int, cond: str, group_by: str, name: str):
    rec = {'repeat':repeat, 'condition':cond, 'grouping':group_by, 'group':str(name), 'n_cells':len(d),
           'mean_tvd':d.tvd.mean(), 'median_tvd':d.tvd.median(), 'p90_tvd':d.tvd.quantile(.9),
           'fate_flip_rate':d.fate_flip.mean()}
    for end in ENDPOINTS:
        rec['mean_signed_'+end] = d['signed_'+end].mean()
        rec['mean_abs_'+end] = d['abs_'+end].mean()
    return rec

def paired_baseline_fate(root: Path, out: Path) -> pd.DataFrame:
    summary, cell_rows, cohort_rows = [],[],[]
    for r in REPEATS:
        ds = {c:read_fate(root,c,r) for c in CONDS}
        common = ds['full'].index.intersection(ds['targeted'].index).intersection(ds['random'].index)
        meta = ds['full'].loc[common,['cell_type','adj_stage_num']]
        for c in ('targeted','random'):
            if not ds[c].loc[common,'cell_type'].equals(meta.cell_type):
                raise ValueError(f'Cell type disagreement in repeat {r}, condition {c}')
        eligible = meta.index[~meta.cell_type.isin(ENDPOINTS)]
        cohort_rows.append({'repeat':r,'n_shared_all':len(common),'n_shared_nonterminal':len(eligible),
                            'n_endpoint_excluded':len(common)-len(eligible)})
        base = ds['full'].loc[eligible,list(ENDPOINTS)].to_numpy(dtype=float)
        types = meta.loc[eligible,'cell_type']
        stages = meta.loc[eligible,'adj_stage_num']
        base_max = base.argmax(axis=1)
        for c in ('targeted','random'):
            q = ds[c].loc[eligible,list(ENDPOINTS)].to_numpy(dtype=float)
            delta = q-base
            row = pd.DataFrame({'cell_id':eligible,'repeat':r,'condition':c,
                                'cell_type':types.to_numpy(), 'stage':stages.to_numpy(),
                                'tvd':np.abs(delta).sum(axis=1)/2,
                                'fate_flip':(q.argmax(axis=1)!=base_max).astype(int)})
            for j,end in enumerate(ENDPOINTS):
                row['signed_'+end]=delta[:,j]
                row['abs_'+end]=np.abs(delta[:,j])
            cell_rows.append(row)
            summary.append(aggregate_fate(row,r,c,'all_shared_nonterminal','all'))
            for typ,sub in row.groupby('cell_type',sort=True):
                summary.append(aggregate_fate(sub,r,c,'cell_type',typ))
            for stage,sub in row.groupby('stage',sort=True):
                summary.append(aggregate_fate(sub,r,c,'stage',stage))
        LOG.info('Fate repeat %d: shared all %d, nonterminal %d',r,len(common),len(eligible))
    cells = pd.concat(cell_rows, ignore_index=True)
    cells.to_csv(out/'tables'/'baseline_fate_per_cell.csv.gz',index=False,compression='gzip')
    pd.DataFrame(cohort_rows).to_csv(out/'tables'/'comparison_cohorts.csv',index=False)
    summaries = pd.DataFrame(summary)
    summaries.to_csv(out/'tables'/'baseline_fate_all_groups.csv',index=False)
    overall = summaries[summaries.grouping=='all_shared_nonterminal'].copy()
    overall.to_csv(out/'tables'/'baseline_fate_by_repeat.csv',index=False)
    for group in ('cell_type','stage'):
        summaries[summaries.grouping==group].to_csv(out/'tables'/f'baseline_fate_by_{group}.csv',index=False)
    piv = cells.pivot(index=['repeat','cell_id'],columns='condition',values='tvd')
    if piv.isna().any().any():
        raise ValueError('Missing matched cells in paired TVD comparison')
    diff = piv['targeted']-piv['random']
    contrasts = diff.groupby(level='repeat').agg(n_cells='count',mean_targeted_minus_random='mean',median_targeted_minus_random='median')
    contrasts.to_csv(out/'tables'/'baseline_fate_contrast_by_repeat.csv')
    make_fate_figures(summaries,overall,out)
    return overall

def full_seed_variability(root: Path, out: Path) -> pd.DataFrame:
    """Full-vs-Full training-seed reference on each matched evaluation cohort.

    This is a descriptive *different-seed* sensitivity reference, not a null
    distribution for sampling and not a substitute for repeated biological data.
    Full_r0/r1/r2 have the same input cells and preprocessing.
    """
    full = {r: read_fate(root, 'full', r) for r in REPEATS}
    records = []
    for r in REPEATS:
        targeted = read_fate(root, 'targeted', r)
        random = read_fate(root, 'random', r)
        ids = full[r].index.intersection(targeted.index).intersection(random.index)
        types = full[r].loc[ids, 'cell_type']
        ids = types.index[~types.isin(ENDPOINTS)]
        for a,b in [(0,1),(0,2),(1,2)]:
            for group in ('all_shared_nonterminal','MyoB1','MyoB2'):
                use = ids if group == 'all_shared_nonterminal' else ids[full[r].loc[ids,'cell_type']==group]
                if not len(use):continue
                p = full[a].loc[use,list(ENDPOINTS)].to_numpy(float)
                q = full[b].loc[use,list(ENDPOINTS)].to_numpy(float)
                scores = 0.5*np.abs(q-p).sum(axis=1)
                flips = (q.argmax(axis=1)!=p.argmax(axis=1)).mean()
                records.append(dict(evaluation_repeat=r,reference_comparison=f'full_r{a}_vs_full_r{b}',
                    source_group=group,n_cells=len(use),mean_tvd=scores.mean(),
                    median_tvd=np.median(scores),p90_tvd=np.quantile(scores,.9),flip_rate=flips))
    df=pd.DataFrame(records)
    df.to_csv(out/'tables'/'full_training_seed_variability.csv',index=False)
    return df

def composition_robustness(out: Path) -> pd.DataFrame:
    """Audit whether large cell types or the depleted group drive headline TVD."""
    d=pd.read_csv(out/'tables'/'baseline_fate_by_cell_type.csv')
    records=[]
    for (r,cond),part in d.groupby(['repeat','condition'],sort=True):
        for label,group in [('all_nonterminal_types',part),('excluding_MyoB1_sources',part[part.group!='MyoB1'])]:
            if len(group)==0:continue
            records.append({'repeat':r,'condition':cond,'source_cohort':label,
                'n_source_types':len(group),'n_cells':int(group.n_cells.sum()),
                'cell_weighted_mean_tvd':np.average(group.mean_tvd,weights=group.n_cells),
                'equal_celltype_mean_tvd':group.mean_tvd.mean()})
    df=pd.DataFrame(records)
    df.to_csv(out/'tables'/'fate_composition_robustness.csv',index=False)
    return df

def make_fate_figures(all_groups: pd.DataFrame, overall: pd.DataFrame,out: Path):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,4.5))
    x=np.arange(3)
    t=overall[overall.condition=='targeted'].sort_values('repeat').mean_tvd.to_numpy()
    r=overall[overall.condition=='random'].sort_values('repeat').mean_tvd.to_numpy()
    ax.bar(x-.18,t,width=.36,label='Targeted')
    ax.bar(x+.18,r,width=.36,label='Random')
    ax.set_xticks(x,['R0','R1','R2']);ax.set_ylabel('Mean total variation distance vs Full')
    ax.set_title('Unperturbed fate prediction stability (shared nonterminal cells)')
    ax.legend();fig.tight_layout();fig.savefig(out/'figures'/'fig_baseline_fate_tvd_by_repeat.png',dpi=200);plt.close(fig)
    rows=all_groups[all_groups.grouping=='cell_type']
    p=rows.pivot_table(index=['repeat','group'],columns='condition',values='mean_tvd')
    diff=(p.targeted-p.random).unstack('repeat').sort_index()
    fig,ax=plt.subplots(figsize=(8,5))
    im=ax.imshow(diff.to_numpy(),cmap='coolwarm',vmin=-max(.001,np.max(np.abs(diff.to_numpy()))),vmax=max(.001,np.max(np.abs(diff.to_numpy()))),aspect='auto')
    ax.set_yticks(range(len(diff.index)),diff.index);ax.set_xticks(range(3),['R0','R1','R2'])
    ax.set_title('Cell-type fate sensitivity: Targeted minus Random')
    fig.colorbar(im,ax=ax,label='Difference in mean TVD');fig.tight_layout()
    fig.savefig(out/'figures'/'fig_fate_tvd_celltype_contrasts.png',dpi=200);plt.close(fig)

def analyze_cbc(root:Path,out:Path):
    p=root/'results'/'p3'/'summary'/'all_cbc.csv'
    s=read_csv(p)
    s['boundary_fraction_of_source']=s['n_valid']/s['n_source']
    s.to_csv(out/'tables'/'cbc_coverage_and_direction.csv',index=False)
    keys=['repeat','source','target','pathway']
    cb=s.pivot(index=keys,columns='condition',values=['mean_cbc','n_valid','boundary_fraction_of_source'])
    contrasts=pd.DataFrame(index=cb.index)
    for v in ['mean_cbc','n_valid','boundary_fraction_of_source']:
        for c in ['targeted','random']:
            contrasts[f'{v}_{c}_minus_full']=cb[(v,c)]-cb[(v,'full')]
        contrasts[f'{v}_targeted_minus_random']=cb[(v,'targeted')]-cb[(v,'random')]
    contrasts.reset_index().to_csv(out/'tables'/'cbc_summary_contrasts.csv',index=False)
    raw=root/'results'/'p3'/'cell_level_cbc'/'all_cell_cbc.csv'
    if raw.is_file():
        d=read_csv(raw)
        rec=[]
        for (r,src,tgt),group in d.groupby(['repeat','source','target'],sort=True):
            m={}
            for c in CONDS:
                one=group[(group.condition==c)&group.cbc.notna()&np.isfinite(group.cbc)]
                if one.cell_id.duplicated().any(): raise ValueError('Duplicated cell-level CBC')
                m[c]=one.set_index('cell_id').cbc
            shared=m['full'].index.intersection(m['targeted'].index).intersection(m['random'].index)
            for c in ('targeted','random'):
                err=m[c].loc[shared]-m['full'].loc[shared]
                rec.append({'repeat':r,'source':src,'target':tgt,'condition':c,
                            'n_paired_boundary_cells':len(shared), 'mean_signed_cbc_change':err.mean() if len(err) else np.nan,
                            'mean_abs_cbc_change':err.abs().mean() if len(err) else np.nan})
        pd.DataFrame(rec).to_csv(out/'tables'/'cbc_matched_boundary_cells.csv',index=False)
    if (root/'results'/'p2'/'comparison'/'paired_latent_time.csv').is_file():
        d=read_csv(root/'results'/'p2'/'comparison'/'paired_latent_time.csv')
        d.to_csv(out/'tables'/'latent_time_pairwise_existing.csv',index=False)
    return s

def analyze_perturbation(root:Path,out:Path):
    path=root/'results'/'p4'/'shared_cell_effects_by_repeat.csv'
    if not path.is_file():
        LOG.warning('P4 summary absent: %s',path);return None
    d=read_csv(path)
    pooled=d[d.source_group=='all_shared_nonterminal'].copy()
    if len(pooled)!=3*2*3:raise ValueError(f'Unexpected P4 pooled summary rows {len(pooled)}')
    for field in ['mean_full','mean_targeted','mean_random','mae_targeted','mae_random','targeted_minus_random_mae']:
        pooled[field+'_pp']=100*pooled[field]
    pooled.to_csv(out/'tables'/'perturbation_effects_existing_p4.csv',index=False)
    pertype=d[d.source_group!='all_shared_nonterminal']
    pertype.to_csv(out/'tables'/'perturbation_by_source_celltype_existing_p4.csv',index=False)
    return pooled

def analyze_cross_layer(root:Path,out:Path):
    """Read P4's original per-cell .csv.gz for cross-layer comparison, no re-fitting."""
    missing=[root/'results'/'p4'/f'{c}_r{r}'/f'{tf}.csv.gz'
             for r in REPEATS for tf in ('MSC','MYOG') for c in CONDS
             if not (root/'results'/'p4'/f'{c}_r{r}'/f'{tf}.csv.gz').is_file()]
    if missing:
        LOG.warning('Per-cell cross-layer NOT run: %d of 18 P4 gzipped CSVs missing.',len(missing))
        return None
    rows=[]
    for r in REPEATS:
        for tf in ('MSC','MYOG'):
            paths={c:root/'results'/'p4'/f'{c}_r{r}'/f'{tf}.csv.gz' for c in CONDS}
            d={}
            for c in CONDS:
                raw = read_csv(paths[c])
                if raw.duplicated(['cell_id','endpoint']).any():
                    raise ValueError(f'Duplicate P4 cell/endpoint rows in {paths[c]}')
                d[c]=raw.set_index(['cell_id','endpoint'])  # uniqueness checked above; pandas 4 compatible
            shared=d['full'].index.intersection(d['targeted'].index).intersection(d['random'].index)
            typ=d['full'].loc[shared,'cell_type']
            eligible=shared[~typ.isin(ENDPOINTS)]
            for c in ('targeted','random'):
                q=d[c].loc[eligible];f=d['full'].loc[eligible]
                if not typ.loc[eligible].equals(q.cell_type):raise ValueError('P4 cell type mismatch')
                if not np.allclose(q.p_perturbed-q.p_unperturbed,q.delta_p,atol=1e-9):
                    raise ValueError(f'P4 delta_p inconsistency: {r} {tf} {c}')
                base_error=(q.p_unperturbed-f.p_unperturbed).abs()
                effect_error=(q.delta_p-f.delta_p).abs()
                table=pd.DataFrame({'repeat':r,'regulator':tf,'condition':c,
                    'cell_id':eligible.get_level_values('cell_id'),
                    'endpoint':eligible.get_level_values('endpoint'),'cell_type':typ.loc[eligible].to_numpy(),
                    'baseline_abs_error':base_error.to_numpy(),'perturb_effect_abs_error':effect_error.to_numpy()})
                rows.append(table)
    if not rows:return None
    a=pd.concat(rows,ignore_index=True)
    groupcols=['repeat','regulator','endpoint','condition']
    summary=a.groupby(groupcols,as_index=False).agg(n_cells=('baseline_abs_error','size'),
                  mean_baseline_error=('baseline_abs_error','mean'),
                  mean_effect_error=('perturb_effect_abs_error','mean'),
                  median_baseline_error=('baseline_abs_error','median'),
                  median_effect_error=('perturb_effect_abs_error','median'))
    summary['effect_minus_baseline']=summary.mean_effect_error-summary.mean_baseline_error
    summary.to_csv(out/'tables'/'cross_layer_baseline_vs_perturbation.csv',index=False)
    # Genuine per-cell association, unlike the previous marginal-mean comparison.
    association=[]
    for key,frame in a.groupby(groupcols,sort=True):
        x=frame.baseline_abs_error
        y=frame.perturb_effect_abs_error
        correlation=x.corr(y,method='spearman') if x.nunique()>1 and y.nunique()>1 else np.nan
        qx=x.quantile(.9); qy=y.quantile(.9)
        # Share of baseline top-decile cells also among effect top-decile cells.
        overlap=float(((x>=qx)&(y>=qy)).sum() / (x>=qx).sum()) if (x>=qx).any() else np.nan
        association.append(dict(zip(groupcols,key),n_cells=len(frame),spearman_rho=correlation,
            top_decile_overlap=overlap,median_baseline_error=x.median(),median_effect_error=y.median()))
    pd.DataFrame(association).to_csv(out/'tables'/'cross_layer_per_cell_association.csv',index=False)
    bytype=a.groupby(groupcols+['cell_type'],as_index=False)[['baseline_abs_error','perturb_effect_abs_error']].mean()
    bytype.to_csv(out/'tables'/'cross_layer_by_celltype.csv',index=False)
    import matplotlib.pyplot as plt
    focus=summary[((summary.regulator=='MYOG')&(summary.endpoint=='MYH3+MyoC'))|((summary.regulator=='MSC')&(summary.endpoint=='PAX7+MyoProg1'))]
    fig,axes=plt.subplots(1,2,figsize=(11,4),sharey=True)
    for ax,tf in zip(axes,('MSC','MYOG')):
        z=focus[focus.regulator==tf]
        labels=[f'R{r}-{c[0].upper()}' for r in REPEATS for c in ('targeted','random')]
        z=z.set_index(['repeat','condition']).loc[[(r,c) for r in REPEATS for c in ('targeted','random')]].reset_index()
        x=np.arange(len(z));ax.bar(x-.2,z.mean_baseline_error,width=.4,label='Baseline fate error')
        ax.bar(x+.2,z.mean_effect_error,width=.4,label='TF effect error')
        ax.set_xticks(x,labels,rotation=45);ax.set_title(tf)
        if tf=='MSC':ax.set_ylabel('Mean absolute error vs Full')
        ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'figures'/'fig_cross_layer_error.png',dpi=200);plt.close(fig)
    return summary

def transition_matrix(root:Path,out:Path):
    """Optional; requires non-Git P2 fitted.h5ad release and RegVelo/CellRank environment."""
    missing=[root/'results'/'p2'/'X'/f'{c}_r{r}'/'fitted.h5ad' for r in REPEATS for c in CONDS if not (root/'results'/'p2'/'X'/f'{c}_r{r}'/'fitted.h5ad').is_file()]
    if missing:
        LOG.warning('Transition matrices NOT run: %d fitted.h5ad missing (not in GitHub ZIP). See README.',len(missing))
        return False
    try:
        import anndata as ad
        from scipy.sparse import issparse
        sys.path.insert(0,str(root/'src'))
        from p3_analysis.fate_inference import build_kernel
    except ImportError as e:
        LOG.warning('Transition matrix requires pinned P2/CellRank environment: %s',e);return False
    rows=[];matrices={}; availability=[];diagnostics=[]
    for r in REPEATS:
        adatas={c:ad.read_h5ad(root/'results'/'p2'/'X'/f'{c}_r{r}'/'fitted.h5ad') for c in CONDS}
        common=adatas['full'].obs_names.intersection(adatas['targeted'].obs_names).intersection(adatas['random'].obs_names)
        classes=sorted(set(adatas['full'].obs['leiden_R_celltype'].astype(str)))
        for c,adata in adatas.items():
            typ=adata.obs['leiden_R_celltype'].astype(str)
            labels=typ.to_numpy()
            K=build_kernel(adata,n_jobs=1).transition_matrix.tocsr()
            inds=adata.obs_names.get_indexer(common)
            if (inds<0).any():raise ValueError('Missing common source cells')
            type_ids={x:np.flatnonzero(labels==x) for x in classes}
            present=np.zeros(adata.n_obs,bool);present[inds]=True
            dest_ids={x:np.flatnonzero((labels==x)&present) for x in classes}
            for x in classes:
                availability.append({'repeat':r,'condition':c,'cell_type':x,'n_available':len(type_ids[x]),'n_shared':len(dest_ids[x])})
            from transition_aggregation import aggregate_by_celltype
            common_indices=adata.obs_names.get_indexer(common)
            common_labels=labels[common_indices]
            type_target,type_shared,shared_valid,_=aggregate_by_celltype(
                K,labels,common_indices,present,classes
            )
            for mode,sub in [('all_available_targets',type_target),('common_target_restricted',type_shared)]:
                valid=(np.ones(len(common_indices),bool) if mode=='all_available_targets'
                       else shared_valid)
                for source in classes:
                    mask=(common_labels==source)&valid
                    if not mask.any():continue
                    row=sub[mask].mean(axis=0)
                    diagnostics.append({'repeat':r,'condition':c,'mode':mode,'source':source,
                                        'n_source_common':int((common_labels==source).sum()),'n_source_valid':int(mask.sum())})
                    for target,v in zip(classes,row):
                        rows.append({'repeat':r,'condition':c,'mode':mode,'source':source,'target':target,'probability':float(v)})
    d=pd.DataFrame(rows);d.to_csv(out/'tables'/'celltype_transition_matrices.csv',index=False)
    pd.DataFrame(availability).to_csv(out/'tables'/'transition_target_availability.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(out/'tables'/'transition_row_validity.csv',index=False)
    contrasts=[]
    for (r,mode),sub in d.groupby(['repeat','mode']):
        keys=['source','target']
        a=sub.pivot(index=keys,columns='condition',values='probability')
        if a.isna().any().any():raise ValueError('Incomplete 10x10 transition matrices')
        for c in ('targeted','random'):
            diff=a[c]-a['full']
            for (s,t),v in diff.items():contrasts.append({'repeat':r,'mode':mode,'condition':c,'source':s,'target':t,'delta_vs_full':v})
    contrasts=pd.DataFrame(contrasts)
    contrasts.to_csv(out/'tables'/'transition_matrix_contrasts.csv',index=False)
    import matplotlib.pyplot as plt
    # Plot 3-repeat averaged contrast (descriptive only)
    mode='all_available_targets'
    for c in ('targeted','random'):
        sub=contrasts[(contrasts.condition==c)&(contrasts['mode']==mode)]
        hm=sub.pivot_table(index='source',columns='target',values='delta_vs_full',aggfunc='mean')
        mat=hm.to_numpy();lim=max(np.max(np.abs(mat)),.0001)
        fig,ax=plt.subplots(figsize=(10,8));im=ax.imshow(mat,aspect='auto',cmap='coolwarm',vmin=-lim,vmax=lim)
        ax.set_xticks(range(len(hm.columns)),hm.columns,rotation=75,ha='right')
        ax.set_yticks(range(len(hm.index)),hm.index)
        ax.set_title(f'{c.capitalize()} – Full: cell-type transition probabilities')
        fig.colorbar(im,ax=ax,label='Mean Δ transition mass');fig.tight_layout()
        fig.savefig(out/'figures'/f'fig_transition_{c}_minus_full.png',dpi=180);plt.close(fig)
    return True

def write_report(out:Path, overall:pd.DataFrame,cbc:pd.DataFrame,p4:pd.DataFrame|None,matrix_complete:bool):
    v=overall.pivot(index='repeat',columns='condition',values='mean_tvd')
    v['targeted_minus_random']=v.targeted-v.random
    lines=['# P5 analysis: generated results','','**Estimand:** within-repeat common nonterminal cells; Full is a reference, not ground truth.','',
           '## Primary finding: unperturbed fate probability stability','',
           '| Repeat | Targeted mean TVD | Random mean TVD | T − R |','|---|---:|---:|---:|']
    for r,row in v.iterrows():lines.append(f'| {r} | {row.targeted:.5f} | {row.random:.5f} | {row.targeted_minus_random:+.5f} |')
    lines += ['','Across all three computational blocks, Targeted induced a larger mean fate-distribution shift on matched nonterminal cells. This does **not** establish a biological effect or isolate cell identity from developmental stage.','',
              '## Additional sensitivity checks', '']
    comp=pd.read_csv(out/'tables'/'fate_composition_robustness.csv')
    for r in REPEATS:
        q=comp[(comp.repeat==r)&(comp.source_cohort=='excluding_MyoB1_sources')].set_index('condition')
        lines.append(f'- R{r}, excluding MyoB1 as source: Targeted mean TVD={q.loc["targeted","cell_weighted_mean_tvd"]:.4f}, Random={q.loc["random","cell_weighted_mean_tvd"]:.4f}, n={int(q.loc["targeted","n_cells"]):,} identical cells.')
    lines+=['', 'The above reweighting/exclusion addresses whether the result is *only* present in remaining MyoB1 source cells, not the confounding of removed cell identities with developmental stage.','',
            '### Full-data repeat variation (unchanged cell input)','','These comparisons change training seed (and may include finite inference-sampling effects). Their dispersion is descriptive, not a formal sampling null distribution.','']
    seed=pd.read_csv(out/'tables'/'full_training_seed_variability.csv')
    lines.append('| Full comparison | Mean TVD, evaluation cohort R0 |')
    lines.append('|---|---:|')
    for key,q in seed[(seed.evaluation_repeat==0)&(seed.source_group=='all_shared_nonterminal')].groupby('reference_comparison'):
        lines.append(f'| {key} | {float(q.mean_tvd.iloc[0]):.5f} |')
    lines += ['', 'The Full-vs-Full variation overlaps the Targeted-vs-Full variability in magnitude; avoid interpreting the contrast as a direct, isolated effect of deleting MyoB1.', '',
              '## Existing P3 direction and coverage checks','',
              'P3 CBC means must be evaluated alongside source boundary counts; the evaluable MyoB1→MyoC1 boundary declines sharply under Targeted thinning. See the generated paired CBC tables.','',
              '## Existing P4 perturbation effects','']
    if p4 is not None:
        for tf,ep in [('MSC','PAX7+MyoProg1'),('MYOG','MYH3+MyoC')]:
            x=p4[(p4.regulator==tf)&(p4.endpoint==ep)].sort_values('repeat')
            lines.append(f'- {tf} → {ep}: targeted-minus-random discrepancy (percentage points): '+', '.join(f'R{int(r["repeat"])}={r["targeted_minus_random_mae_pp"]:+.2f}' for _,r in x.iterrows())+'.')
        lines.append('No uniformly larger Targeted discrepancy in P4 across the three computational blocks.')
    lines+=['','## Optional whole-network transition analysis','',
            ('The 10×10 transition matrix analysis was completed. See tables/transition_matrix_contrasts.csv.' if matrix_complete else 'Not executed. P2 fitted.h5ad files and CellRank/RegVelo were not available; the ordinary GitHub ZIP excludes these large assets.'),
            '', '## Direct per-cell cross-layer comparison','',
            ('Computed from all 18 original P4 compressed tables; see tables/cross_layer_per_cell_association.csv.' if (out/'tables'/'cross_layer_per_cell_association.csv').is_file() else 'Not executed in this environment. Requires the 18 P4 per-cell .csv.gz files.'),
            '', '## Critical constraints','',
            '- Three computational repeats are not biological replicates. No cell-level p-values are interpreted as biological inference.',
            '- The fixed 606-gene panel and prior GRN include information extracted from the original full study.',
            '- Deleting MyoB1 also disproportionately removes earlier PCW stages; the current random control is not stage matched.',
            '- P3 CellRank fate probabilities are conditional on fixed terminal anchors; P4 blockades manipulate learned TF regulatory weights, not empirical knockouts.',
            '- Paired cells participated in model fitting; these results measure within-dataset prediction sensitivity, not held-out predictive accuracy.']
    (out/'reports'/'findings.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',help='Extracted BioModelSimProject repository root')
    parser.add_argument('--matrix',action='store_true',help='Also compute 10x10 transition matrices if P2 fitted H5ADs exist')
    args=parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO,format='%(levelname)s: %(message)s')
    src=Path(__file__).resolve().parent
    root=locate_project(args.project_root,Path(__file__).resolve())
    out=src/'outputs';(out/'tables').mkdir(parents=True,exist_ok=True);(out/'figures').mkdir(exist_ok=True);(out/'reports').mkdir(exist_ok=True)
    LOG.info('Project root: %s',root)
    baseline=paired_baseline_fate(root,out)
    composition_robustness(out)
    full_seed_variability(root,out)
    cbc=analyze_cbc(root,out)
    p4=analyze_perturbation(root,out)
    cross=analyze_cross_layer(root,out)
    matrix_done=transition_matrix(root,out) if args.matrix else False
    write_report(out,baseline,cbc,p4,matrix_done)
    LOG.info('Results: %s',out)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
