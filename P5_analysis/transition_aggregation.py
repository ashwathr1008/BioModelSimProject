"""Pure numerical helper for aggregating CellRank cell-to-cell transition kernels.

No RegVelo or CellRank dependency: can be checked on small synthetic sparse kernels.
The transition matrix itself must still be computed from the fitted AnnData.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix, issparse


def aggregate_by_celltype(K, labels, common_indices, common_mask, classes):
    """Return (all_targets, common_targets_conditioned, valid_rows, target availability).

    K: (N,N) row-stochastic CellRank transition kernel.
    labels: N cell-type labels in AnnData order.
    common_indices: positions of cells present in every condition in THIS AnnData.
    common_mask: boolean length-N marking the same common cell IDs as destinations.
    classes: stable ordered set of cell types.

    Each output is aligned to common_indices; output columns are classes. The second
    matrix excludes destination cells absent in any condition and re-normalizes.
    This is a sensitivity diagnostic; it is NOT the actual original Markov chain.
    """
    if not issparse(K):K=csr_matrix(K)
    K=K.tocsr()
    labels=np.asarray(labels).astype(str)
    common_indices=np.asarray(common_indices,dtype=int)
    common_mask=np.asarray(common_mask,dtype=bool)
    n=len(labels)
    if K.shape!=(n,n) or common_mask.shape!=(n,):raise ValueError('Kernel/labels dimensions disagree')
    if len(set(classes))!=len(classes):raise ValueError('Duplicate cell types')
    if len(common_indices)==0 or np.any(common_indices<0) or np.any(common_indices>=n):
        raise ValueError('Empty or invalid shared source indices')
    if not common_mask[common_indices].all():raise ValueError('Shared source IDs must also be shared destinations')
    if not np.isfinite(K.data).all() or np.any(K.data < -1e-8):raise ValueError('Invalid transition entries')
    rs=np.asarray(K.sum(axis=1)).ravel()
    if not np.allclose(rs,1,atol=1e-5,rtol=0):raise ValueError('Kernel rows must sum to one')
    all_cols=[];shared_cols=[];counts=[]
    for cl in classes:
        destination=(labels==str(cl))
        shared=destination & common_mask
        all_cols.append(np.asarray(K[:,destination].sum(axis=1)).ravel()[common_indices])
        shared_cols.append(np.asarray(K[:,shared].sum(axis=1)).ravel()[common_indices])
        counts.append((int(destination.sum()),int(shared.sum())))
    all_targets=np.column_stack(all_cols)
    common_targets=np.column_stack(shared_cols)
    retained_mass=common_targets.sum(axis=1)
    valid=retained_mass>1e-12
    conditional=np.full_like(common_targets,np.nan,dtype=float)
    conditional[valid]=common_targets[valid]/retained_mass[valid,None]
    if not np.allclose(all_targets.sum(axis=1),1,atol=1e-5,rtol=0):
        raise ValueError('Aggregated probability mass must sum to one')
    return all_targets,conditional,valid,np.asarray(counts)
