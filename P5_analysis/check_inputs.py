#!/usr/bin/env python3
"""Inventory original P1-P4 files required by P5, without modifying them."""
import argparse
import importlib.util
from pathlib import Path

CONDITIONS=('full','targeted','random')
REPEATS=range(3)

def audit(root:Path):
    groups={
       'P3 fate (required)': [root/'results/p3/fate'/f'{c}_r{r}_fate.csv' for r in REPEATS for c in CONDITIONS],
       'P3 CBC summaries (required)': [root/'results/p3/summary/all_cbc.csv'],
       'P3 cell-level CBC (recommended)': [root/'results/p3/cell_level_cbc/all_cell_cbc.csv'],
       'P4 pooled effects (required)': [root/'results/p4/shared_cell_effects_by_repeat.csv'],
       'P4 per-cell effects (cross-layer)': [root/'results/p4'/f'{c}_r{r}'/f'{tf}.csv.gz' for r in REPEATS for c in CONDITIONS for tf in ('MSC','MYOG')],
       'P2 fitted AnnData (10x10 matrix)': [root/'results/p2/X'/f'{c}_r{r}'/'fitted.h5ad' for r in REPEATS for c in CONDITIONS],
    }
    for group,paths in groups.items():
        found=sum(p.is_file() for p in paths)
        print(f'{group}: {found}/{len(paths)}')
        if found<len(paths):
            for p in paths:
                if not p.is_file(): print('  missing',p.relative_to(root))
    print('Optional 10x10 matrix libraries:')
    for name in ('anndata','cellrank','scvelo','regvelo','scanpy'):
        print(f'  {name}: {"installed" if importlib.util.find_spec(name) else "missing"}')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--project-root',required=True)
    args=parser.parse_args()
    root=Path(args.project_root).expanduser().resolve()
    if not (root/'config/experiment.yaml').is_file():
        parser.error(f'Expected BioModelSimProject root with config/experiment.yaml; got {root}')
    audit(root)
