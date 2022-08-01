"""Module for analysis of pattern similarity evidence of reactivation."""

from pathlib import Path
import random

import click
import numpy as np
import scipy.stats as stats
import scipy.io as sio
import scipy.spatial.distance as sd
import statsmodels.api as sm
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


def rdm_reactivation_stats_split(
    subjects, rdms, dfs, accuracy_filter, category_filter, min_count=2
):
    """Statistics of similarity between pre-exposure and study patterns."""
    results_list = []
    for rdm, df in zip(rdms, dfs):
        # get matching trials
        category_match = (df['AC'].fillna(0).to_numpy() == accuracy_filter)
        accuracy_match = (df['category'].to_numpy() == category_filter)
        match = category_match & accuracy_match

        # if no matching trials, statistics are undefined
        if np.count_nonzero(match) < min_count:
            res = pd.Series({'self': np.nan, 'within': np.nan})
            results_list.append(res)
            continue

        # filter the dataframe and RDM
        df = df.loc[match]
        rdm = rdm[np.ix_(match, match)]

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
            }
        )
        results_list.append(res)
    results = pd.DataFrame(results_list, index=subjects)
    results['item'] = results['self'] - results['within']
    return results


def reactivation_stats_split(subjects, roi_rdms, dfs):
    """Reactivation statistics by category and accuracy."""
    categories = dfs[0]['category'].unique()
    stats_list = []
    for roi, rdms in roi_rdms.items():
        for a in [1, 0]:
            for c in categories:
                results = rdm_reactivation_stats_split(
                    subjects, rdms, dfs, accuracy_filter=a, category_filter=c
                )
                results['roi'] = roi
                results['correct'] = a
                results['category'] = c
                stats_list.append(results)
    stats = pd.concat(stats_list)
    stats.index.rename('subject', inplace=True)
    return stats


def robust_slope(x, y, min_count=5):
    """Calculate slope based on robust regression"""
    if len(x) >= min_count:
        # ROI pair correlation for this accuracy and category bin
        exog = sm.add_constant(x)
        rlm = sm.RLM(y, exog, M=sm.robust.norms.TukeyBiweight())
        rlm_res = rlm.fit()
        slope = rlm_res.params[1]
    else:
        # correlation undefined if less than 3 included trials
        slope = np.nan
    return slope


def item_roi_correlation(subjects, roi1, roi2, roi_rdms, dfs, split_accuracy=False):
    """Correlation between item reactivation in two ROIs."""
    results_list = []
    for subject, rdm1, rdm2, df in zip(subjects, roi_rdms[roi1], roi_rdms[roi2], dfs):
        # Fisher z of correlation
        item1 = np.diag(np.arctanh(1 - rdm1))
        item2 = np.diag(np.arctanh(1 - rdm2))

        # category of the A item
        category = df['category'].to_numpy()
        categories = np.unique(category)

        if split_accuracy:
            # AC test performance for each group
            # no-response trials scored as incorrect
            correct = df['AC'].fillna(0).to_numpy()
            for a in [1, 0]:
                for c in categories:
                    include = (correct == a) & (category == c)
                    slope = robust_slope(item1[include], item2[include])

                    # package stats for this ROI pair and bin
                    r = pd.Series(
                        {
                            'subject': subject,
                            'roi1': roi1,
                            'roi2': roi2,
                            'correct': a,
                            'category': c,
                            'slope': slope,
                        }
                    )
                    results_list.append(r)
        else:
            for c in categories:
                include = category == c
                slope = robust_slope(item1[include], item2[include])
                r = pd.Series(
                    {
                        'subject': subject,
                        'roi1': roi1,
                        'roi2': roi2,
                        'category': c,
                        'slope': slope
                    }
                )
                results_list.append(r)
    item_corr = pd.DataFrame(results_list)
    return item_corr


def trial_model_correlation(model_rdm, neural_rdm):
    """Correlation with a model per trial."""
    n_items = model_rdm.shape[0]
    rho = np.empty(n_items)
    include = ~np.eye(n_items, dtype=bool)
    for i in range(n_items):
        # correlation between model distance and neural distance,
        # excluding self-similarity
        model_vec = model_rdm[include[i, :], i]
        neural_vec = neural_rdm[include[i, :], i]
        rho[i], _ = stats.spearmanr(model_vec, neural_vec)
    return rho


def subject_item_model_correlation(subject, item_rdm, integ_rdm, df, model):
    """Reactivation model correlation by AC accuracy."""
    # extract self-similarity as a measure of item reactivation
    r_self = np.diag(np.arctanh(1 - item_rdm))

    # create model RDM for AC semantics
    a_vectors = get_item_vectors(df['item1'], model)
    c_vectors = get_item_vectors(df['item3'], model)
    ac_vectors = a_vectors + c_vectors
    ac_rdm = sd.squareform(sd.pdist(ac_vectors, 'correlation'))

    # calculate model correlation for each trial
    r_ac = trial_model_correlation(ac_rdm, integ_rdm)

    results_list = []
    for a in [1, 0]:
        include = df['AC'].fillna(0).to_numpy() == a
        x = r_self[include]
        y = r_ac[include]
        slope = robust_slope(x, y)
        res = pd.Series(
            {
                'subject': subject,
                'correct': a,
                'slope': slope,
            }
        )
        results_list.append(res)
    results = pd.DataFrame(results_list)
    return results


def item_model_correlation(subjects, item_rdms, integ_rdms, dfs, model):
    """Correlation between item reactivation and model by AC accuracy."""
    results_list = []
    for item_roi, item_roi_rdms in item_rdms.items():
        for integ_roi, integ_roi_rdms in integ_rdms.items():
            inputs = zip(subjects, item_roi_rdms, integ_roi_rdms, dfs)
            for subject, item_rdm, integ_rdm, df in inputs:
                res = subject_item_model_correlation(
                    subject, item_rdm, integ_rdm, df, model
                )
                res['item_roi'] = item_roi
                res['integ_roi'] = integ_roi
                results_list.append(res)
    results = pd.concat(results_list)
    return results


def load_model_mat(model_file):
    """Read a model RDM from a MAT-file."""
    # load the MAT-file
    mat = sio.loadmat(model_file)
    if 'rdm' not in mat:
        raise ValueError(f'Model RDM file not in standard format: {model_file}')
    rdm = mat['rdm']
    if 'vectors' not in mat:
        vectors = None
    else:
        vectors = mat['vectors']

    # check that the rdm is sensible
    assert np.allclose(np.diag(rdm), 0)
    assert rdm.shape[0] == rdm.shape[1]

    # get items in list format
    f = mat['items']
    if len(f) == 1:
        items = [i[0] for i in f[0]]
    elif type(f[0]) == np.str_:
        items = [i.strip() for i in f]
    else:
        items = [i[0][0] for i in f]
    assert rdm.shape[0] == len(items)
    return items, vectors, rdm


def pool_index(trial_items, pool_items_list):
    """
    Get the index of each item in the full pool.

    Parameters
    ----------
    trial_items : pandas.Series
        The item presented on each trial.

    pool_items_list : list or numpy.ndarray
        List of items in the full pool.

    Returns
    -------
    item_index : pandas.Series
        Index of each item in the pool. Trials with items not in the
        pool will be NaN.
    """
    trial_items = pd.Series(trial_items)
    pool_map = dict(zip(pool_items_list, np.arange(len(pool_items_list))))
    item_index = trial_items.map(pool_map)
    return item_index


def get_item_vectors(items, model):
    """Get vectors for a set of items."""
    if not len(np.unique(model['items'])) == len(model['items']):
        raise ValueError(f'Items in model are not unique.')

    item_index = pool_index(items, model['items'])
    if np.any(np.isnan(item_index)):
        n = np.count_nonzero(np.isnan(item_index))
        raise ValueError(f'{n} item(s) not found in model.')
    item_vectors = model['vectors'][item_index, :]
    return item_vectors


@click.command()
@click.argument("mat_file", type=click.Path(exists=True))
@click.argument("npz_file", type=click.Path())
def convert_model(mat_file, npz_file):
    """Convert a model MAT-file to npz format."""
    items, vectors, rdm = load_model_mat(mat_file)
    np.savez(npz_file, items=items, vectors=vectors, rdm=rdm)
