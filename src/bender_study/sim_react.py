"""Module for analysis of pattern similarity evidence of reactivation."""

import random
import numpy as np
import scipy.stats as stats
from scipy.spatial.distance import cdist
from mvpa2.measures.base import Measure


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
            perm_ind[cat_ind] = random.sample(cat_ind, len(cat_ind))
        rand_ind.append(perm_ind)
    return rand_ind

def perm_z(stat_perm):
    """Calculate a z-statistic from permutation results."""
    p = np.mean(stat_perm >= stat_perm[0])
    e = 1.0 / len(stat_perm)
    if p > (1 - e):
        p = 1 - e
    return stats.norm.ppf(1 - p)

class SimCatReact(Measure):
    """Test for reactivation of category patterns.

    Calculate similarity of items between two phases, and contrast
    within- and between-category similarity values. Estimate null by
    permuting phase 2 items.

    """
    
    def __init__(self, split, cat, n_perm, output='full'):
        Measure.__init__(self)
        self.split_ind = tuple([np.nonzero(split==u)[0]
                                for u in np.unique(split)])
        self.output = output
        self.n_perm = n_perm

        # permute items across categories
        n_item = len(cat)
        rand_ind = [random.sample(range(n_item), n_item)
                    for i in range(n_perm)]
        rand_ind.insert(0, range(n_item))

        # create a matrix labeling pairs:
        # 1 - same category
        # 2 - different category
        # will contrast 1 and 2.
        cat_mat = np.zeros((n_item, n_item))
        for i in range(n_item):
            for j in range(n_item):
                if cat[i] == cat[j]:
                    cat_mat[i,j] = 1
                else:
                    cat_mat[i,j] = 2

        # get within indices and between indices for each permutation,
        # to get the correct values of the matrix that cdist will
        # produce
        self.w_ind = []
        self.b_ind = []
        for ind in rand_ind:
            # scramble labels for one phase only (ignores item
            # relationships)
            r_cat_mat = cat_mat[:,ind]
            self.w_ind.append(np.nonzero(r_cat_mat==1))
            self.b_ind.append(np.nonzero(r_cat_mat==2))
    
    def __call__(self, dataset):
        rho = 1 - cdist(dataset[self.split_ind[0]],
                        dataset[self.split_ind[1]], 'correlation')

        stat_perm = np.empty((self.n_perm + 1,))
        stat_perm.fill(np.nan)
        for i in range(self.n_perm + 1):
            if (self.w_ind[i].size == 0 or self.b_ind[i].size == 0):
                continue
            
            within = np.mean(rho[self.w_ind[i]])
            between = np.mean(rho[self.b_ind[i]])
            stat_perm[i] = within - between

        if np.any(np.isnan(stat_perm)):
            raise ValueError('Undefined statistic.')

        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)

class SimCatReactSME(Measure):
    """Test whether reactivation of category patterns predicts memory.

    Calculate similarity of items between two phases, and contrast
    within- and between-category similarity values, separately for
    recalled and forgotten items. Pairs of items that include one
    recalled item and one forgotten item are excluded. Hypothesis is
    that there is a larger separation between categories on
    subsequently recalled trials. Estimate null by permuting phase 2
    items within category.

    """
    
    def __init__(self, split, cat, correct, n_perm, output='full'):
        Measure.__init__(self)
        self.split_ind = tuple([np.nonzero(split==u)[0]
                                for u in np.unique(split)])
        self.n_perm = n_perm
        self.output = output

        # scramble items within category
        rand_ind = perm_within(cat, n_perm)
        n_item = len(cat)
        rand_ind.insert(0, np.arange(n_item))

        # create a matrix with within- and between-category labels,
        # and a matrix with neither recalled and both recalled labels
        cat_mat = np.zeros((n_item, n_item))
        rec_mat = np.zeros((n_item, n_item))
        for i in range(n_item):
            for j in range(n_item):
                if cat[i] == cat[j]:
                    cat_mat[i,j] = 1
                else:
                    cat_mat[i,j] = 2
                if correct[i]==1 and correct[j]==1:
                    rec_mat[i,j] = 1
                elif correct[i]==0 and correct[j]==0:
                    rec_mat[i,j] = 2

        # indices to calculate interaction: (wr - br) - (wf - bf)
        self.wr_ind = []
        self.br_ind = []
        self.wf_ind = []
        self.bf_ind = []
        for ind in rand_ind:
            r_cat = cat_mat[:,ind]
            r_rec = rec_mat[:,ind]
            self.wr_ind.append(np.nonzero(np.logical_and(r_cat==1, r_rec==1)))
            self.br_ind.append(np.nonzero(np.logical_and(r_cat==2, r_rec==1)))
            self.wf_ind.append(np.nonzero(np.logical_and(r_cat==1, r_rec==2)))
            self.bf_ind.append(np.nonzero(np.logical_and(r_cat==2, r_rec==2)))
            
    def __call__(self, dataset):
        rho = 1 - cdist(dataset[self.split_ind[0]],
                        dataset[self.split_ind[1]], 'correlation')
        
        stat_perm = np.empty((self.n_perm + 1,))
        stat_perm.fill(np.nan)
        for i in range(self.n_perm + 1):
            if (self.wr_ind[i][0].size==0 or self.br_ind[i][0].size==0 or
                self.wf_ind[i][0].size==0 or self.bf_ind[i][0].size==0):
                continue
            
            wr = np.mean(rho[self.wr_ind[i]])
            br = np.mean(rho[self.br_ind[i]])
            wf = np.mean(rho[self.wf_ind[i]])
            bf = np.mean(rho[self.bf_ind[i]])
            stat_perm[i] = (wr - br) - (wf - bf)

        if np.any(np.isnan(stat_perm)):
            raise ValueError('Undefined statistic.')
            
        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)
                                                                    
class SimSubcatReact(Measure):
    """Test for reactivation of subcategory patterns.

    Calculate similarity of items between two phases, and for each
    category, contrast within- and between-subcategory similarity
    values (both are within-category). Average the contrast across
    categories. Estimate null by permuting phase 2 items within each
    category.

    """
    
    def __init__(self, split, cat, subcat, n_perm, output='full'):
        Measure.__init__(self)
        self.split_ind = tuple([np.nonzero(split==u)[0]
                                for u in np.unique(split)])
        self.output = output
        self.n_perm = n_perm

        # create label matrix
        n_item = len(cat)
        subcat_mat = np.zeros((n_item, n_item))
        for i in range(n_item):
            for j in range(n_item):
                if subcat[i] == subcat[j]:
                    subcat_mat[i,j] = 1
                else:
                    subcat_mat[i,j] = 2

        # get indices within each category
        self.cat_ind = []
        self.w_ind = []
        self.b_ind = []
        for i, c in enumerate(np.unique(cat)):
            # get labels for this category
            cat_ind = np.nonzero(cat == c)[0]
            self.cat_ind.append(cat_ind)
            m_subcat = subcat_mat[np.ix_(cat_ind, cat_ind)]

            # scramble items within category
            n_cat_item = len(cat_ind)
            rand_ind = [random.sample(range(n_cat_item), n_cat_item)
                        for i in range(n_perm)]
            rand_ind.insert(0, range(n_cat_item))

            # within- and between-subcategory indices
            w_ind = []
            b_ind = []
            for ind in rand_ind:
                r_subcat = m_subcat[:,ind]
                w_ind.append(np.nonzero(r_subcat==1))
                b_ind.append(np.nonzero(r_subcat==2))
            self.w_ind.append(w_ind)
            self.b_ind.append(b_ind)
            
    def __call__(self, dataset):
        stat_perm_all = np.empty((len(self.cat_ind), self.n_perm + 1))
        stat_perm_all.fill(np.nan)
        samples1 = dataset.samples[self.split_ind[0]]
        samples2 = dataset.samples[self.split_ind[1]]
        for i, cat_ind in enumerate(self.cat_ind):
            if (self.w_ind[i][0][0].size==0 or
                self.b_ind[i][0][0].size==0):
                continue

            rho = 1 - cdist(samples1[cat_ind], samples2[cat_ind],
                            'correlation')
            
            for j in range(self.n_perm + 1):
                within = np.mean(rho[self.w_ind[i][j]])
                between = np.mean(rho[self.b_ind[i][j]])
                stat_perm_all[i,j] = within - between
        stat_perm = np.nanmean(stat_perm_all, 0)
        if np.any(np.isnan(stat_perm)):
            raise ValueError('Undefined statistic.')
        
        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)
        
class SimSubcatReactSME(Measure):
    """Test whether reactivation of subcategory patterns predicts memory.
    
    Calculate similarity of items between two phases, and contrast
    within- and between-subcategory similarity values, separately for
    recalled and forgotten items. Pairs of items that include one
    recalled item and one forgotten item are excluded. Hypothesis is
    that there is a larger separation between subcategories on
    subsequently recalled trials. Contrast is averaged across
    categories. Estimate null by permuting phase 2 items within
    subcategory (this keeps the subcategory structure intact and
    maintains the number of remembered and forgotten items within
    subcategory constant).

    """

    def __init__(self, split, cat, subcat, correct, n_perm=100, output='full'):
        Measure.__init__(self)
        self.split_ind = tuple([np.nonzero(split==u)[0]
                                for u in np.unique(split)])
        self.n_perm = n_perm
        self.output = output

        # subcategory and recall labels
        n_item = len(cat)
        rec_mat = np.zeros((n_item, n_item))
        subcat_mat = np.zeros((n_item, n_item))
        for i in range(n_item):
            for j in range(n_item):
                if subcat[i] == subcat[j]:
                    subcat_mat[i,j] = 1
                else:
                    subcat_mat[i,j] = 2
                if correct[i]==1 and correct[j]==1:
                    rec_mat[i,j] = 1
                elif correct[i]==0 and correct[j]==0:
                    rec_mat[i,j] = 2

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
            m_subcat = subcat_mat[np.ix_(cat_ind, cat_ind)]
            m_rec = rec_mat[np.ix_(cat_ind, cat_ind)]

            # scramble items within subcategory
            rand_ind = perm_within(subcat[cat_ind], n_perm)
            rand_ind.insert(0, np.arange(len(cat_ind)))

            # indices to calculate interaction: (wr - br) - (wf - bf)
            iwr = []
            ibr = []
            iwf = []
            ibf = []
            for ind in rand_ind:
                r_subcat = m_subcat[:,ind]
                r_rec = m_rec[:,ind]
                iwr.append(np.nonzero(np.logical_and(r_subcat==1, r_rec==1)))
                ibr.append(np.nonzero(np.logical_and(r_subcat==2, r_rec==1)))
                iwf.append(np.nonzero(np.logical_and(r_subcat==1, r_rec==2)))
                ibf.append(np.nonzero(np.logical_and(r_subcat==2, r_rec==2)))
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
            if (self.wr_ind[i][0][0].size==0 or
                self.br_ind[i][0][0].size==0 or
                self.wf_ind[i][0][0].size==0 or
                self.bf_ind[i][0][0].size==0):
                continue
                
            rho = 1 - cdist(samples1[cat_ind], samples2[cat_ind],
                            'correlation')
            for j in range(self.n_perm + 1):
                wr = np.mean(rho[self.wr_ind[i][j]])
                br = np.mean(rho[self.br_ind[i][j]])
                wf = np.mean(rho[self.wf_ind[i][j]])
                bf = np.mean(rho[self.bf_ind[i][j]])
                stat_perm_all[i,j] = (wr - br) - (wf - bf)
        stat_perm = np.nanmean(stat_perm_all, 0)
        if np.any(np.isnan(stat_perm)):
            raise ValueError('Undefined statistic.')
                
        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)

class SimItemReact(Measure):
    """Test for reactivation of item patterns.

    Calculate similarity of items between two phases, and contrast
    same-item and different-item, same-(sub)category similarity
    values. Estimate null by permuting phase 2 items within
    (sub)category.

    Items must be sorted the same way in both datasets.

    """
    
    def __init__(self, split, cat, n_perm, stat='react', output='full'):
        Measure.__init__(self)
        self.split_ind = tuple([np.nonzero(split==u)[0]
                                for u in np.unique(split)])
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
            rand_ind = [random.sample(range(n_cat_item), n_cat_item)
                        for i in range(n_perm)]
            rand_ind.insert(0, range(n_cat_item))

            # within- and between-item indices
            w_ind = []
            b_ind = []
            for ind in rand_ind:
                r_item = m_item[:,ind]
                w_ind.append(np.nonzero(r_item==1))
                b_ind.append(np.nonzero(r_item==0))
            self.w_ind.append(w_ind)
            self.b_ind.append(b_ind)

    def __call__(self, dataset):
        stat_perm_all = np.empty((len(self.cat_ind), self.n_perm + 1))
        stat_perm_all.fill(np.nan)
        samples1 = dataset.samples[self.split_ind[0]]
        samples2 = dataset.samples[self.split_ind[1]]
        for i, cat_ind in enumerate(self.cat_ind):
            if (self.w_ind[i][0][0].size==0 or
                self.b_ind[i][0][0].size==0):
                continue

            rho = 1 - cdist(samples1[cat_ind], samples2[cat_ind],
                            'correlation')
            
            for j in range(self.n_perm + 1):
                w = np.mean(rho[self.w_ind[i][j]])
                b = np.mean(rho[self.b_ind[i][j]])
                if self.stat == 'react':
                    stat_perm_all[i,j] = w - b
                elif self.stat == 'suppress':
                    stat_perm_all[i,j] = b - w
        if np.any(np.isnan(stat_perm_all)):
            raise ValueError('statistic is undefined.')
        stat_perm = np.mean(stat_perm_all, 0)
        
        if self.output == 'full':
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
    
    def __init__(self, split, cat, correct, n_perm, stat='react_sme',
                 output='full'):
        Measure.__init__(self)
        self.split_ind = tuple([np.nonzero(split==u)[0]
                                for u in np.unique(split)])
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
                    item_mat[i,j] = 1
                else:
                    item_mat[i,j] = 2
                if correct[i]==1 and correct[j]==1:
                    rec_mat[i,j] = 1
                elif correct[i]==0 and correct[j]==0:
                    rec_mat[i,j] = 2

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
            rand_ind = [random.sample(range(n_cat_item), n_cat_item)
                        for i in range(n_perm)]
            rand_ind.insert(0, range(n_cat_item))

            # indices to calculate interaction: (wr - br) - (wf - bf)
            iwr = []
            ibr = []
            iwf = []
            ibf = []
            for ind in rand_ind:
                # scramble labels for both datasets (this preserves
                # item identity while breaking relationship to memory)
                r_item = m_item[np.ix_(ind,ind)]
                r_rec = m_rec[np.ix_(ind,ind)]
                wr = np.nonzero(np.logical_and(r_item==1, r_rec==1))
                br = np.nonzero(np.logical_and(r_item==2, r_rec==1))
                wf = np.nonzero(np.logical_and(r_item==1, r_rec==2))
                bf = np.nonzero(np.logical_and(r_item==2, r_rec==2))
                if (wr[0].size==0 or br[0].size==0 or
                    wf[0].size==0 or bf[0].size==0):
                    raise ValueError('One or more conditions have no observations.')
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
            rho = 1 - cdist(samples1[cat_ind], samples2[cat_ind],
                            'correlation')
            
            for j in range(self.n_perm + 1):
                wr = np.mean(rho[self.wr_ind[i][j]])
                br = np.mean(rho[self.br_ind[i][j]])
                wf = np.mean(rho[self.wf_ind[i][j]])
                bf = np.mean(rho[self.bf_ind[i][j]])

                if self.stat == 'react_sme':
                    # more item reactivation -> better memory
                    stat_perm_all[i,j] = (wr - br) - (wf - bf)
                elif self.stat == 'suppress_sme':
                    # more item suppression -> better memory
                    stat_perm_all[i,j] = (br - wr) - (bf - wf)
        stat_perm = np.nanmean(stat_perm_all, 0)
        if np.any(np.isnan(stat_perm)):
            raise ValueError('Undefined statistic.')
        
        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)
