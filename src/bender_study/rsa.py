"""Module for analysis of pattern similarity evidence of reactivation."""

from pathlib import Path
import random

import numpy as np
import scipy.stats as stats
import pandas as pd


def unique_ordered(a):
    _, idx = np.unique(a, return_index=True)
    return a[np.sort(idx)]


def perm_within(cat, n_perm):
    """Create permutation indices to scramble within category."""
    ucat = np.unique(cat)
    rand_ind = []
    for i in range(n_perm):
        perm_ind = np.zeros((len(cat),), dtype=int)
        for j in ucat:
            cat_ind = np.nonzero(cat == j)[0]
            perm_ind[cat_ind] = random.sample(list(cat_ind), len(cat_ind))
        rand_ind.append(perm_ind)
    return rand_ind


def perm_z(stat_perm):
    """Calculate a z-statistic from permutation results."""
    p = np.mean(stat_perm >= stat_perm[0])
    e = 1.0 / len(stat_perm)
    if p > (1 - e):
        p = 1 - e
    return stats.norm.ppf(1 - p)


def load_roi_rdms(rdm_dir, subjects, rois, analyses, clusters, suffix):
    """Load RDMs for multiple ROIs."""
    rdm_dir = Path(rdm_dir)
    rdms = {
        roi: [
            np.loadtxt(
                rdm_dir / f'{analysis}_{cluster}_{suffix}' / f'bender_{s}_rdm.txt'
            ) for s in subjects
        ] for analysis, cluster, roi in zip(analyses, clusters, rois)
    }
    return rdms


def rdm_reactivation_stats(subjects, rdms, dfs):
    """Statistics of similarity between pre-exposure and study patterns."""
    results_list = []
    for rdm, df in zip(rdms, dfs):
        # selector for self pairs
        pair_self = np.eye(df.shape[0], dtype=bool)

        # selector for within-category pairs
        category = df['category'].to_numpy()
        pair_within = category == category[:, np.newaxis]

        # transform to Fisher z of correlation
        z = np.arctanh(1 - rdm)

        # calculate statistics
        res = pd.Series(
            {
                'self': np.mean(z[pair_self]),
                'within': np.mean(z[~pair_self & pair_within]),
                'between': np.mean(z[~pair_within]),
            }
        )
        results_list.append(res)
    results = pd.DataFrame(results_list, index=subjects)
    results['item'] = results['self'] - results['within']
    results['category'] = results['within'] - results['between']
    return results


def reactivation_stats(subjects, roi_rdms, dfs):
    """Reactivation statistics for multiple ROIs."""
    stats_list = []
    for roi, rdms in roi_rdms.items():
        results = rdm_reactivation_stats(subjects, rdms, dfs)
        stats_list.append(results)
    stats = pd.concat(stats_list, keys=roi_rdms.keys())
    stats.index.rename(['roi', 'subject'], inplace=True)
    return stats
