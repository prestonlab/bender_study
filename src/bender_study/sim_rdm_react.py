import random
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats.stats import pearsonr
from mvpa2.measures.base import Measure

class SimRDMReact(Measure):

    def __init__(self, split_attr, n_item, n_perm=5000):
        Measure.__init__(self)
        self.split_attr = split_attr
        self.n_perm = n_perm

        self.rand_ind = []
        for i in range(n_perm):
            rand_ind = range(n_item)
            random.shuffle(rand_ind)
            self.rand_ind.append(rand_ind)

    def __call__(self, dataset):

        # split into the chunks we are comparing
        split = dataset.sa[self.split_attr].value
        usplit = np.unique(split)
        if len(usplit) != 2:
            raise ValueError('Split attribute must have two unique values.')
        ds1 = dataset[split == usplit[0],:]
        ds2 = dataset[split == usplit[1],:]

        # calculate pairwise similarity of stimuli in each phase
        rdm_vec1 = 1 - pdist(ds1.samples, 'correlation')
        rdm_vec2 = 1 - pdist(ds2.samples, 'correlation')
        rdm_mat2 = squareform(rdm_vec2)
        rho_perm = np.zeros(self.n_perm + 1)
        if np.any(np.isnan(rdm_vec1)) or np.any(np.isnan(rdm_vec2)):
            # probably no variation between voxels in searchlight,
            # probably because we're outside the brain. Just return
            # zeros for everything
            return tuple(rho_perm)

        # actual correlation
        rho_data = pearsonr(rdm_vec1, rdm_vec2)[0]
        rho_perm[0] = rho_data

        # permuted correlation
        for i, rand_ind in enumerate(self.rand_ind):
            rho_perm[i+1] = pearsonr(rdm_vec1, squareform(rdm_mat2[np.ix_(rand_ind, rand_ind)]))[0]
        return tuple(rho_perm)
