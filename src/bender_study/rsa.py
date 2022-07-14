"""Module for analysis of pattern similarity evidence of reactivation."""

from pathlib import Path
import random

import numpy as np
import scipy.stats as stats
from scipy.spatial.distance import cdist, pdist, squareform
from mvpa2.measures.base import Measure
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


class SimItemReact(Measure):
    """Test for reactivation of item patterns.

    Calculate similarity of items between two phases, and contrast
    same-item and different-item, same-(sub)category similarity
    values. Estimate null by permuting phase 2 items within
    (sub)category.

    Items must be sorted the same way in both datasets.

    """

    def __init__(self, split, cat, n_perm, stat="react", output="full"):
        Measure.__init__(self)
        self.split_ind = tuple([np.nonzero(split == u)[0] for u in np.unique(split)])
        self.n_perm = n_perm
        self.stat = stat
        self.output = output

        # label pairs by item (items must be ordered the same way in
        # both datasets)
        same_item = np.eye(len(cat))

        # get indices within each category
        self.cat_ind = []
        self.w_ind = []
        self.b_ind = []
        for i, c in enumerate(np.unique(cat)):
            # get labels for this category
            cat_ind = np.nonzero(cat == c)[0]
            self.cat_ind.append(cat_ind)
            m_item = same_item[np.ix_(cat_ind, cat_ind)]

            # scramble items within category
            n_cat_item = len(cat_ind)
            rand_ind = [
                random.sample(range(n_cat_item), n_cat_item) for _ in range(n_perm)
            ]
            rand_ind.insert(0, list(range(n_cat_item)))

            # within- and between-item indices
            w_ind = []
            b_ind = []
            for ind in rand_ind:
                r_item = m_item[:, ind]
                w_ind.append(np.nonzero(r_item == 1))
                b_ind.append(np.nonzero(r_item == 0))
            self.w_ind.append(w_ind)
            self.b_ind.append(b_ind)

    def __call__(self, dataset):
        stat_perm_all = np.empty((len(self.cat_ind), self.n_perm + 1))
        stat_perm_all.fill(np.nan)
        samples1 = dataset.samples[self.split_ind[0]]
        samples2 = dataset.samples[self.split_ind[1]]
        for i, cat_ind in enumerate(self.cat_ind):
            if self.w_ind[i][0][0].size == 0 or self.b_ind[i][0][0].size == 0:
                continue

            rho = 1 - cdist(samples1[cat_ind], samples2[cat_ind], "correlation")

            for j in range(self.n_perm + 1):
                w = np.mean(rho[self.w_ind[i][j]])
                b = np.mean(rho[self.b_ind[i][j]])
                if self.stat == "react":
                    stat_perm_all[i, j] = w - b
                elif self.stat == "suppress":
                    stat_perm_all[i, j] = b - w
        if np.any(np.isnan(stat_perm_all)):
            raise ValueError("statistic is undefined.")
        stat_perm = np.mean(stat_perm_all, 0)

        if self.output == "full":
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)


class SimItemReactSME(Measure):
    """Test whether reactivation of items patterns predicts memory.

    Calculate similarity of items between two phases, and contrast
    same-item and different-item, same-(sub)category similarity
    values, separately for recalled and forgotten items. Pairs of
    items that include one recalled item and one forgotten item are
    excluded. Hypothesis is that patterns are more item-specific on
    subsequently recalled trials. Estimate null by permuting phase 2
    items within (sub)category.

    """

    def __init__(self, split, cat, correct, n_perm, stat="react_sme", output="full"):
        Measure.__init__(self)
        self.split_ind = tuple([np.nonzero(split == u)[0] for u in np.unique(split)])
        self.n_perm = n_perm
        self.stat = stat
        self.output = output

        # item and recall labels
        n_item = len(cat)
        item_mat = np.zeros((n_item, n_item))
        rec_mat = np.zeros((n_item, n_item))
        for i in range(n_item):
            for j in range(n_item):
                if i == j:
                    item_mat[i, j] = 1
                else:
                    item_mat[i, j] = 2
                if correct[i] == 1 and correct[j] == 1:
                    rec_mat[i, j] = 1
                elif correct[i] == 0 and correct[j] == 0:
                    rec_mat[i, j] = 2

        # get contrast indices for each category
        self.cat_ind = []
        self.wr_ind = []
        self.br_ind = []
        self.wf_ind = []
        self.bf_ind = []
        for i, c in enumerate(np.unique(cat)):
            # get labels for this category
            cat_ind = np.nonzero(cat == c)[0]
            self.cat_ind.append(cat_ind)
            m_item = item_mat[np.ix_(cat_ind, cat_ind)]
            m_rec = rec_mat[np.ix_(cat_ind, cat_ind)]

            # scramble items within category
            n_cat_item = len(cat_ind)
            rand_ind = [
                random.sample(range(n_cat_item), n_cat_item) for _ in range(n_perm)
            ]
            rand_ind.insert(0, list(range(n_cat_item)))

            # indices to calculate interaction: (wr - br) - (wf - bf)
            iwr = []
            ibr = []
            iwf = []
            ibf = []
            for ind in rand_ind:
                # scramble labels for both datasets (this preserves
                # item identity while breaking relationship to memory)
                r_item = m_item[np.ix_(ind, ind)]
                r_rec = m_rec[np.ix_(ind, ind)]
                wr = np.nonzero(np.logical_and(r_item == 1, r_rec == 1))
                br = np.nonzero(np.logical_and(r_item == 2, r_rec == 1))
                wf = np.nonzero(np.logical_and(r_item == 1, r_rec == 2))
                bf = np.nonzero(np.logical_and(r_item == 2, r_rec == 2))
                if (
                    wr[0].size == 0
                    or br[0].size == 0
                    or wf[0].size == 0
                    or bf[0].size == 0
                ):
                    raise ValueError("One or more conditions have no observations.")
                iwr.append(wr)
                ibr.append(br)
                iwf.append(wf)
                ibf.append(bf)
            self.wr_ind.append(iwr)
            self.br_ind.append(ibr)
            self.wf_ind.append(iwf)
            self.bf_ind.append(ibf)

    def __call__(self, dataset):
        stat_perm_all = np.empty((len(self.cat_ind), self.n_perm + 1))
        stat_perm_all.fill(np.nan)
        samples1 = dataset.samples[self.split_ind[0]]
        samples2 = dataset.samples[self.split_ind[1]]
        for i, cat_ind in enumerate(self.cat_ind):
            rho = 1 - cdist(samples1[cat_ind], samples2[cat_ind], "correlation")

            for j in range(self.n_perm + 1):
                wr = np.mean(rho[self.wr_ind[i][j]])
                br = np.mean(rho[self.br_ind[i][j]])
                wf = np.mean(rho[self.wf_ind[i][j]])
                bf = np.mean(rho[self.bf_ind[i][j]])

                if self.stat == "react_sme":
                    # more item reactivation -> better memory
                    stat_perm_all[i, j] = (wr - br) - (wf - bf)
                elif self.stat == "suppress_sme":
                    # more item suppression -> better memory
                    stat_perm_all[i, j] = (br - wr) - (bf - wf)
        stat_perm = np.nanmean(stat_perm_all, 0)
        if np.any(np.isnan(stat_perm)):
            raise ValueError("Undefined statistic.")

        if self.output == "full":
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)


class SimModelCond(Measure):
    def __init__(self, cond, model_rdms, n_perm, output="full"):
        """Contrast similarity between model and neural data by condition.

        Calculate model-neural similarity for each condition, for each
        model RDM, and calculate the interaction between model and
        condition. Pairs of items that are in different conditions are
        excluded. Hypothesis is that there is a greater difference
        between conditions 1 and 2 for model 1 than there is for model
        2. Estimate null by permuting items within condition.

        """

        Measure.__init__(self)
        self.n_perm = n_perm
        self.output = output
        self.cond = cond

        # scramble similarity values across items, within category
        # (breaking relationship to model within each condition;
        # should also break interaction)
        n_item = len(cond)
        rand_ind = perm_within(cond, n_perm)
        rand_ind.insert(0, np.arange(n_item))

        self.mat_ind = []
        self.model = []
        self.n_model = len(model_rdms)
        self.ucond = np.unique(cond)
        self.n_cond = len(self.ucond)
        self.include = []

        # get model and scrambled models for each condition
        for c in self.ucond:
            item_ind = cond == c
            sub_ind = np.ix_(item_ind, item_ind)
            self.mat_ind.append(item_ind)

            c_list = []
            c_include = []
            for rdm in model_rdms:
                # model matrix for this condition
                c_mat = rdm[sub_ind]

                # indices to vectorize matrix
                n_item = c_mat.shape[0]
                vec_ind = np.triu_indices(n_item, 1)

                # items/pairs to include for this model
                vec = c_mat[vec_ind]
                vec_include = np.logical_not(np.isnan(vec))
                item_include = np.sum(np.isnan(c_mat), 0) < (n_item - 1)
                all_include = np.zeros(len(cond), dtype=bool)
                all_include[item_ind] = item_include

                m_list = []
                for i in rand_ind:
                    r_ind = np.ix_(i, i)
                    rdm_rand = rdm[r_ind][sub_ind]
                    inc_rand = all_include[i][item_ind]
                    vec_rand = squareform(rdm_rand[np.ix_(inc_rand, inc_rand)])
                    m_list.append(stats.rankdata(vec_rand))
                # array with actual and permuted models
                c_list.append(np.asarray(m_list))
                c_include.append(vec_include)
            # [condition][model] list of model arrays/included vector indices
            self.model.append(c_list)
            self.include.append(c_include)

    def __call__(self, dataset):

        # calculate data RDM, convert to ranking vector
        rho = np.empty((self.n_perm + 1, self.n_cond, self.n_model))
        rho.fill(np.nan)
        for c, c_ind in enumerate(self.mat_ind):
            if np.count_nonzero(c_ind) < 3:
                continue

            # items corresponding to this condition
            rdm = pdist(dataset.samples[c_ind, :], "correlation")
            x = stats.rankdata(rdm)
            x = x.reshape((1, len(x)))
            for m, ymat in enumerate(self.model[c]):
                # calculate all correlations for this condition and
                # model
                include = self.include[c][m]
                rho[:, c, m] = 1 - cdist(x[:, include], ymat, "correlation")

        # calculate interaction: (C1M1 - C1M2) - (C2M1 - C2M2)
        stat_perm = (rho[:, 0, 0] - rho[:, 0, 1]) - (rho[:, 1, 0] - rho[:, 1, 1])
        if np.any(np.isnan(stat_perm)):
            raise ValueError("Undefined statistic.")

        if self.output == "full":
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)


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
                'item': np.mean(z[pair_self]),
                'within': np.mean(z[~pair_self & pair_within]),
                'between': np.mean(z[~pair_within]),
            }
        )
        results_list.append(res)
    results = pd.DataFrame(results_list, index=subjects)
    return results
