"""Module for analysis of pattern similarity correlations with models."""

import random
import numpy as np
from numpy import linalg
import scipy.stats as stats
from scipy.spatial.distance import cdist, pdist, squareform
from mvpa2.measures.base import Measure
import scipy.optimize as optim

from bender_study.sim_react import perm_within, perm_z, unique_ordered


class SimRDM(Measure):
    """
    Calculate similarity between data and a target RDM.
    
    Calculate an RDM for the data including all pairs of
    items. Correlate the data RDM with the target RDM. Hypothesis is
    that the correlation is greater than chance. Estimate null by
    creating permuted models with shuffled items.
    """
    
    def __init__(self, cat, model_rdm, n_perm, output='full', rand_ind=None):
        Measure.__init__(self)
        self.n_perm = n_perm
        self.output = output

        # actual and random item indices
        if rand_ind is None:
            rand_ind = perm_within(cat, n_perm)
        rand_ind.insert(0, np.arange(len(cat)))

        # create all model vectors, rank transformed
        rand_model = []
        for ind in rand_ind:
            r_model = squareform(model_rdm[np.ix_(ind, ind)])
            rand_model.append(stats.rankdata(r_model))
        self.model = np.asarray(rand_model)

    def __call__(self, dataset):
        # calculate data RDM
        data_rdm = stats.rankdata(pdist(dataset.samples, 'correlation'))
        xmat = data_rdm.reshape((1,len(data_rdm)))

        # calculate all data-model correlations
        stat_perm = 1 - cdist(xmat, self.model, 'correlation').squeeze()
        if np.any(np.isnan(stat_perm)):
            raise ValueError('statistic is undefined.')
        
        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)

class SimRDMPartial(Measure):

    def __init__(self, cat, model_rdm, control_rdms, n_perm,
                 fit='ls', output='full'):
        Measure.__init__(self)
        self.n_perm = n_perm
        self.output = output
        self.fit = fit

        rand_ind = perm_within(cat, n_perm)
        rand_ind.insert(0, np.arange(len(cat)))

        # control models
        control_vecs = np.asarray([stats.rankdata(squareform(rdm))
                                   for rdm in control_rdms]).T
        intercept = np.ones((control_vecs.shape[0],1))
        self.control_mat = np.hstack((control_vecs, intercept))

        # model residuals
        model_vec = stats.rankdata(squareform(model_rdm))
        if fit == 'ls':
            beta = linalg.lstsq(self.control_mat, model_vec)[0]
        elif fit == 'nnls':
            beta = optim.nnls(self.control_mat, model_vec)[0]
        else:
            raise ValueError('Unsupported fit type: {}'.format(fit))
        resid = model_vec - self.control_mat.dot(beta)
        self.resid = resid[:,None].T
        self.rand_ind = rand_ind

    def __call__(self, dataset):
        # calculate data RDM
        data_mat = squareform(stats.rankdata(pdist(dataset.samples,
                                                   'correlation')))

        # get residuals for randomized data
        data_resid = []
        for ind in self.rand_ind:
            data_vec = squareform(data_mat[np.ix_(ind,ind)])
            if self.fit == 'ls':
                beta = linalg.lstsq(self.control_mat, data_vec)[0]
            elif self.fit == 'nnls':
                beta = optim.nnls(self.control_mat, data_vec)[0]
            data_resid.append(data_vec - self.control_mat.dot(beta))
            
        # correlate with the residualized model
        stat_perm = 1 - cdist(np.asarray(data_resid), self.resid, 'correlation').squeeze()
        if np.any(np.isnan(stat_perm)):
            raise ValueError('statistic is undefined.')
        
        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)

class SimRDMPartial2(Measure):

    def __init__(self, cat, model_rdm, control_rdms, n_perm,
                 fit='ls', output='full'):
        Measure.__init__(self)
        self.n_perm = n_perm
        self.output = output
        self.fit = fit

        rand_ind = perm_within(cat, n_perm)
        rand_ind.insert(0, np.arange(len(cat)))

        # control models in real order
        control_vecs = np.asarray([stats.rankdata(squareform(rdm))
                                   for rdm in control_rdms]).T
        intercept = np.ones((control_vecs.shape[0],1))
        self.control_mat = np.hstack((control_vecs, intercept))

        # residualize the (random) model of interest
        resid = []
        for ind in rand_ind:
            rand_model = stats.rankdata(squareform(model_rdm[np.ix_(ind,ind)]))
            if fit == 'ls':
                beta = linalg.lstsq(self.control_mat, rand_model)[0]
            elif fit == 'nnls':
                beta = optim.nnls(self.control_mat, rand_model)[0]
            else:
                raise ValueError('Unsupported fit type: {}'.format(fit))
            res = rand_model - self.control_mat.dot(beta)
            resid.append(res)
        self.resid = np.asarray(resid)

    def __call__(self, dataset):
        if np.count_nonzero(dataset.fa.include) < 10:
            if self.output in ['full', 'roi']:
                return (tuple(np.zeros(self.n_perm + 1)), 0)
            else:
                return (0, 0)
        dataset = dataset[:,dataset.fa.include]
        
        # calculate data RDM
        data_vec = stats.rankdata(pdist(dataset.samples, 'correlation'))

        # data residuals
        if self.fit == 'ls':
            beta = linalg.lstsq(self.control_mat, data_vec)[0]
        elif self.fit == 'nnls':
            beta = optim.nnls(self.control_mat, data_vec)[0]
        data_resid = data_vec - self.control_mat.dot(beta)

        # correlate with the residualized (random) model
        xmat = data_resid.reshape((1,len(data_resid)))
        
        stat_perm = 1 - cdist(xmat, self.resid, 'correlation').squeeze()
        if np.any(np.isnan(stat_perm)):
            raise ValueError('statistic is undefined.')
        
        if self.output in ['full', 'roi']:
            return (tuple(stat_perm), 1)
        else:
            return (perm_z(stat_perm), 1)
        
class SimRDMSME(Measure):
    """Test whether correlation with an RDM predicts memory.

    Calculate pairwise similarity for correct and incorrect items, and
    correlation with an RDM for those two sets of items (without
    regard to category). Contrast correlation for correct
    vs. incorrect items. Hypothesis is that there is a greater
    correlation with the RDM on correct trials. Estimate null by
    scrambling items within category (this maintains any accuracy
    differences between category).

    """
    
    def __init__(self, cat, correct, model_rdm, n_perm, output='full'):
        Measure.__init__(self)
        self.n_perm = n_perm
        self.output = output

        # scramble within category, relative to correct and incorrect
        # labels (same as scrambling the labels within category)
        rand_ind = perm_within(cat, n_perm)
        rand_ind.insert(0, np.arange(len(cat)))

        self.mat_ind = []
        self.model = []
        uacc = np.unique(correct)
        self.n_acc = len(uacc)

        for a in uacc:
            a_rand = []
            item_ind = correct == a
            self.mat_ind.append(item_ind)
            for ind in rand_ind:
                r_ind = np.ix_(ind, ind)
                model_rand = model_rdm[r_ind][np.ix_(item_ind, item_ind)]
                a_rand.append(stats.rankdata(squareform(model_rand)))
            self.model.append(np.asarray(a_rand))

    def __call__(self, dataset):

        # calculate spearman correlation for each permutation
        # and accuracy condition
        rho = np.empty((self.n_perm+1, self.n_acc))
        rho.fill(np.nan)
        for a, a_ind in enumerate(self.mat_ind):
            if np.count_nonzero(a_ind) < 3:
                continue
            
            rdm = pdist(dataset.samples[a_ind,:], 'correlation')
            x = stats.rankdata(rdm)
            x = x.reshape((1,len(x)))
            ymat = self.model[a]
            rho[:,a] = 1 - cdist(x, ymat, 'correlation')

        # difference between correct and incorrect
        stat_perm = rho[:,1] - rho[:,0]
        if np.any(np.isnan(stat_perm)):
            raise ValueError('Statistic is undefined.')

        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)

class SimRDMCatSME(Measure):
    """Test whether correlation with an RDM predicts memory.
    
    For each category and accuracy condition, calculate an
    RDM. Correlate each of those RDMs with a model RDM. Contrast
    correlation for correct vs. incorrect items, and average that
    difference across category. Hypothesis is that there is a greater
    correlation with the RDM on correct trials. Estimate null by
    scrambling items within category.

    """

    def __init__(self, cat, correct, model_rdm, n_perm, output='full'):
        Measure.__init__(self)
        self.n_perm = n_perm
        self.output = output
        self.cat = cat
        self.correct = correct

        # scramble within category, relative to correct and incorrect
        # labels (same as scrambling the labels within category)
        rand_ind = perm_within(cat, n_perm)
        rand_ind.insert(0, np.arange(len(cat)))

        self.mat_ind = []
        self.model = []
        ucat = np.unique(cat)
        uacc = np.unique(correct)
        self.n_cat = len(ucat)
        self.n_acc = len(uacc)
        
        for c in ucat:
            c_rand_list = []
            c_mat_ind = []
            for a in uacc:
                # indices for this submatrix of the RDM
                item_ind = np.logical_and(cat == c, correct == a)
                c_mat_ind.append(item_ind)
                sub_ind = np.ix_(item_ind, item_ind)
                
                # permuted model dissimilarity vectors
                c_rand = []
                for i in rand_ind:
                    r_ind = np.ix_(i, i)
                    c_rand.append(stats.rankdata(squareform(model_rdm[r_ind][sub_ind])))
                c_rand_list.append(np.asarray(c_rand))
            self.mat_ind.append(c_mat_ind)
            self.model.append(c_rand_list)

    def __call__(self, dataset):

        # calculate spearman correlation for each permutation,
        # category, and accuracy condition
        rho = np.empty((self.n_perm+1, self.n_cat, self.n_acc))
        rho.fill(np.nan)
        for c, c_ind in enumerate(self.mat_ind):
            for a, a_ind in enumerate(c_ind):
                if np.count_nonzero(a_ind) < 3:
                    # minimum samples to have a defined data-model
                    # correlation is 3
                    continue
                
                rdm = pdist(dataset.samples[a_ind,:], 'correlation')
                x = stats.rankdata(rdm)
                x = x.reshape((1,len(x)))
                ymat = self.model[c][a]
                rho[:,c,a] = 1 - cdist(x, ymat, 'correlation')

        # difference between correct and incorrect, averaged over category
        stat_perm = np.nanmean(rho[:,:,1] - rho[:,:,0], 1)
        if np.any(np.isnan(stat_perm)):
            raise ValueError('Statistic is undefined.')
        
        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)

class SimRDMContrast(Measure):

    def __init__(self, cond, contrast, model_rdm, n_perm, output='full'):
        """Contrast similarity between model and neural data by condition."""
        
        Measure.__init__(self)
        self.n_perm = n_perm
        self.output = output
        self.cond = cond
        self.contrast = contrast

        # scramble similarity values across items
        n_item = len(cond)
        rand_ind = [random.sample(range(n_item), n_item) for i in range(n_perm)]

        self.mat_ind = []
        self.model = []
        self.ucond = np.unique(cond)
        self.n_cond = len(self.ucond)

        # get model and scrambled models for each condition
        for c in self.ucond:
            item_ind = cond == c
            sub_ind = np.ix_(item_ind, item_ind)

            c_list = []
            c_list.append(stats.rankdata(squareform(model_rdm[sub_ind])))
            for i in rand_ind:
                r_ind = np.ix_(i, i)
                c_list.append(stats.rankdata(squareform(model_rdm[r_ind][sub_ind])))
            self.mat_ind.append(item_ind)
            self.model.append(np.asarray(c_list))

    def __call__(self, dataset):

        data_problem = False
        rho = np.zeros((self.n_perm+1, self.n_cond))
        for c, c_ind in enumerate(self.mat_ind):
            # calculate this sub-matrix of the full RDM
            rdm = pdist(dataset.samples[c_ind,:], 'correlation')
            if np.any(np.isnan(rdm)):
                data_problem = True
                break
            x = stats.rankdata(rdm)
            x = x.reshape((1,len(x)))

            # get model RDM for this condition
            ymat = self.model[c]

            # calculate rank correlation
            rho[:,c] = 1 - cdist(x, ymat, 'correlation')

        if data_problem:
            if self.output == 'full':
                return tuple(np.zeros(self.n_perm+1))
            else:
                return np.float64(0)
            
        # calculate contrast for each model (contrast must be in the
        # order of the sorted conditions)
        stat_perm = np.dot(rho, self.contrast)
        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)

class SimModelCond(Measure):

    def __init__(self, cond, model_rdms, n_perm, output='full'):
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
                item_include = np.sum(np.isnan(c_mat), 0) < (n_item-1)
                all_include = np.zeros(len(cond), dtype=bool)
                all_include[item_ind] = item_include
                
                m_list = []
                for i in rand_ind:
                    r_ind = np.ix_(i, i)
                    rdm_rand = rdm[r_ind][sub_ind]
                    inc_rand = all_include[i][item_ind]
                    vec_rand = squareform(rdm_rand[np.ix_(inc_rand,inc_rand)])
                    m_list.append(stats.rankdata(vec_rand))
                # array with actual and permuted models
                c_list.append(np.asarray(m_list))
                c_include.append(vec_include)
            # [condition][model] list of model arrays/included vector indices
            self.model.append(c_list)
            self.include.append(c_include)

    def __call__(self, dataset):

        # calculate data RDM, convert to ranking vector
        rho = np.empty((self.n_perm+1, self.n_cond, self.n_model))
        rho.fill(np.nan)
        for c, c_ind in enumerate(self.mat_ind):
            if np.count_nonzero(c_ind) < 3:
                continue
            
            # items corresponding to this condition
            rdm = pdist(dataset.samples[c_ind,:], 'correlation')
            x = stats.rankdata(rdm)
            x = x.reshape((1,len(x)))
            for m, ymat in enumerate(self.model[c]):
                # calculate all correlations for this condition and
                # model
                include = self.include[c][m]
                rho[:,c,m] = 1 - cdist(x[:,include], ymat, 'correlation')

        # calculate interaction: (C1M1 - C1M2) - (C2M1 - C2M2)
        stat_perm = (rho[:,0,0]-rho[:,0,1]) - (rho[:,1,0]-rho[:,1,1])
        if np.any(np.isnan(stat_perm)):
            raise ValueError('Undefined statistic.')
        
        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)

class SimModelContrast(Measure):
    """Contrast different models.
    
    Calculate an RDM for neural data and compare to models, then
    contrast their match to the data. Contrast within different
    (optional) categories, with differences averaged across
    categories. Scrambling is done within category, and also within a
    separate grouping variable (e.g. subcategory). This allows for
    preserving subcategory structure in the permuted runs.

    """
    
    def __init__(self, cat, subcat, models, n_perm, contrast, output='full'):
        Measure.__init__(self)
        self.n_perm = n_perm
        self.n_model = len(models)
        self.output = output
        self.contrast = contrast

        # categories x models x permutations RDVs
        self.models = []

        ucat = np.unique(cat)
        for c in ucat:
            # subcat labels within category
            c_subcat = subcat[cat == c]
            
            # use indices relative to category items
            n_item = np.count_nonzero(cat == c)
            rand_ind = [range(n_item)]
            for j in range(n_perm):
                # scramble within subcat
                rand_ind.append(perm_within(c_subcat, n_perm))

            # create model list
            c_list = []
            for mat in models:
                c_mat = mat[np.ix_(cat==c, cat==c)]
                r_list = []
                for ind in rand_ind:
                    c_mat_rand = c_mat[np.ix_(ind, ind)]
                    r_list.append(stats.rankdata(squareform(c_mat_rand)))
                c_list.append(np.asarray(r_list))
            self.models.append(c_list)
        import pdb
        pdb.set_trace()

    def __call__(self, dataset):
        # calculate data RDM, convert to ranking vector
        n_cat = len(self.ucat)

        rho = np.zeros((self.n_perm+1, self.n_model, n_cat))
        for i, c in enumerate(self.ucat):
            # get RDV for this category
            data_rdm = pdist(dataset.samples[self.cat==c,:], 'correlation')
            if np.any(np.isnan(data_rdm)):
                raise ValueError('NaN in data RDM.')
            data_rdm = stats.rankdata(data_rdm)
            ymat = data_rdm.reshape((1,len(data_rdm)))

            for j in range(self.n_model):
                # compare the data to all permutations of this model
                xmat = self.models[i,j]
                rho[:,i,j] = 1 - cdist(xmat, ymat, 'correlation')

        import pdb
        pdb.set_trace()
        
class SimMaxContrastRDM(Measure):

    def __init__(self, cat, models, n_perm, contrast, output='full',
                 rand_ind_in=None):
        Measure.__init__(self)
        self.n_perm = n_perm
        self.output = output
        self.cat = cat
        self.ucat = np.unique(cat)
        self.contrast = contrast

        # dict of vis and sem with a list of vectorized models
        self.models = {}

        # set random indices (must be the same across models)
        ucat = np.unique(cat)
        if rand_ind_in is not None:
            # only have enough indices for one category, so use for
            # all included categories
            rand_ind = {}
            for i, c in enumerate(ucat):
                rand_ind[c] = rand_ind_in
                rand_ind[c].insert(0, range(len(rand_ind_in[0])))
        else:
            # scramble within each category
            rand_ind = {}
            for c in ucat:
                n_item = np.count_nonzero(cat == c)
                c_ind = [range(n_item)]
                for j in range(n_perm):
                    c_ind.append(random.sample(range(n_item), n_item))
                rand_ind[c] = c_ind

        # create model and random model dicts
        for model_type in models.keys():
            m_list = []
            for mat in models[model_type]:
                c_list = []
                for c in ucat:
                    # get model similarity for this category
                    c_mat = mat[np.ix_(cat==c, cat==c)]

                    # replace any missing data with the median within
                    # this category
                    vec_ind = np.triu_indices(c_mat.shape[0], 1)
                    c_mat[np.isnan(c_mat)] = np.nanmedian(c_mat[vec_ind])

                    r_list = []
                    for ind in rand_ind[c]:
                        c_mat_rand = c_mat[np.ix_(ind, ind)]
                        r_list.append(stats.rankdata(squareform(c_mat_rand)))
                    c_list.append(np.asarray(r_list))
                m_list.append(c_list)
            self.models[model_type] = m_list

    def __call__(self, dataset):
        # calculate data RDM, convert to ranking vector
        xmat = []
        n_cat = len(self.ucat)
        for c in self.ucat:
            data_rdm = pdist(dataset.samples[self.cat==c,:], 'correlation')
            if np.any(np.isnan(data_rdm)):
                raise ValueError('NaN in data RDM.')
            data_rdm = stats.rankdata(data_rdm)
            xmat.append(data_rdm.reshape((1,len(data_rdm))))

        rho = {}
        for model_type in self.models.keys():
            n_model = len(self.models[model_type])
            rho_type = np.zeros((self.n_perm+1, n_model, n_cat))
            for i in range(n_model):
                for j in range(n_cat):
                    # calculate model-data RDM correlation for this
                    # model (within this model type) and category
                    ymat = self.models[model_type][i][j]
                    rho_type[:,i,j] = 1 - cdist(xmat[j], ymat, 'correlation')

            # take max over all models within the model type
            rho[model_type] = np.max(rho_type, 1)

        if '-' in self.contrast:
            # this is a contrast of two model types
            mlist = self.contrast.split('-')
            stat_perm = np.mean(rho[mlist[0]] - rho[mlist[1]], 1)
        else:
            # just testing if one model type is significantly
            # correlated
            stat_perm = np.mean(rho[self.contrast], 1)
            
        # statistic is a difference of model types, averaged over
        # categories
        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)
            
class SimFitContrastRDM(Measure):
    """Fit multiple models to an RDM and contrast them."""
    
    def __init__(self, cat, models, n_perm, contrast,
                 fit='nnls', output='full'):
        Measure.__init__(self)
        self.n_perm = n_perm
        self.output = output
        self.cat = cat
        self.ucat = np.unique(cat)
        self.contrast = contrast
        self.fit = fit

        # because we're building a model of the full RDM in a given
        # sphere, will replace missing values instead of excluding
        # them, so models that do define those values can contribute
        self.designs = {}

        # scramble within each category, keeping same scrambles across
        # all models
        rand_ind = []
        for i, c in enumerate(self.ucat):
            # indices to scramble the within-category sub-matrix
            n_item = np.count_nonzero(cat == c)
            c_ind = [range(n_item)]
            for j in range(n_perm):
                c_ind.append(random.sample(range(n_item), n_item))
            rand_ind.append(c_ind)
            
        # create model and random model dicts
        for model_type in models.keys():
            # make a [cat][perm] nested list of design matrices
            n_model = len(models[model_type])
            
            # for each category and permutation, want a
            # design matrix for the nnls fit that is [features x models]
            c_list = []
            for i, c in enumerate(self.ucat):
                # prepare model RDMs for this category
                c_mats = []
                for mat in models[model_type]:
                    # get model similarity for this category
                    c_mat = mat[np.ix_(cat==c, cat==c)]

                    # replace any missing data with the median within
                    # this category
                    vec_ind = np.triu_indices(c_mat.shape[0], 1)
                    c_mat[np.isnan(c_mat)] = np.nanmedian(c_mat[vec_ind])
                    c_mats.append(c_mat)

                # prepare design matrices
                r_list = []
                n_item = c_mat.shape[0]
                n_feat = (n_item ** 2 - n_item) / 2
                for r_ind in rand_ind[i]:
                    design = np.ones((n_feat, n_model + 1))
                    for j, mat in enumerate(c_mats):
                        design[:,j+1] = squareform(mat[np.ix_(r_ind,r_ind)])
                    r_list.append(design)
                c_list.append(r_list)
            self.designs[model_type] = c_list
            
    def __call__(self, dataset):
        
        # calculate data RDM within each included category
        y_list = []
        n_cat = len(self.ucat)
        data_problem = False
        for i in range(n_cat):
            data_rdm = pdist(dataset.samples[self.cat==self.ucat[i],:],
                             'correlation')
            if np.any(np.isnan(data_rdm)):
                data_problem = True
                break
            y_list.append(data_rdm)
        
        # if any nans encountered in the data rdm, just output zeros
        if data_problem:
            if self.output == 'full':
                return tuple(np.zeros(self.n_perm+1))
            else:
                return np.float64(0)

        # calculate correlation between the data and the best-fitting
        # model-based RDM, separately for each model type
        rho = {}
        for model_type in self.designs.keys():
            n_model = len(self.designs[model_type])
            rho_type = np.zeros((n_cat, self.n_perm+1))
            for i in range(n_cat):
                # data to fit
                y = y_list[i]
                for j in range(self.n_perm+1):
                    # fit model
                    design = self.designs[model_type][i][j]
                    if self.fit == 'nnls':
                        x = optim.nnls(design, y)[0]
                    else:
                        x = linalg.lstsq(design, y)[0]

                    # assess fit
                    yhat = np.dot(design, x)
                    r = stats.spearmanr(yhat, y)[0]
                    if not np.isnan(r):
                        rho_type[i,j] = r
            rho[model_type] = rho_type

        if '-' in self.contrast:
            # this is a contrast of two model types
            mlist = self.contrast.split('-')
            stat_perm = np.mean(rho[mlist[0]] - rho[mlist[1]], 0)
        else:
            # just testing if one model type is significantly
            # correlated
            stat_perm = np.mean(rho[self.contrast], 0)
            
        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)
        
class SimFitCVContrastRDM(Measure):

    def __init__(self, cat, models, n_perm, type1, type2, output='full',
                 fit='nnls', n_xval=10):
        Measure.__init__(self)
        self.n_perm = n_perm
        self.output = output
        self.cat = cat
        self.ucat = unique_ordered(cat)
        self.type1 = type1
        self.type2 = type2
        self.fit = fit
        self.n_xval = n_xval

        # because we're building a model of the full RDM in a given
        # sphere, will replace missing values instead of excluding
        # them, so models that do define those values can contribute
        self.designs = {}

        # scramble within each category, keeping same scrambles across
        # all models
        rand_ind = []
        self.xval_ind = []
        for i, c in enumerate(self.ucat):
            # indices to scramble the within-category sub-matrix
            c_ind = []
            n_item = np.count_nonzero(cat == c)
            for j in range(n_perm):
                c_ind.append(random.sample(range(n_item), n_item))
            rand_ind.append(c_ind)

            x_ind = []
            for x in range(n_xval):
                fold = {}
                # choose a random half of the data
                fold['test'] = np.array(random.sample(range(n_item), n_item/2))
                fold['train'] = np.setdiff1d(np.arange(n_item), fold['test'])
                x_ind.append(fold)
            self.xval_ind.append(x_ind)
            
        # create model and random model dicts
        for model_type in models.keys():
            # make a [cat][xval][perm] nested list of train and test matrices
            n_model = len(models[model_type])
            
            # for each model type, category, and permutation, want a
            # design matrix for the nnls fit that is [features x models]
            c_list = []
            for i, c in enumerate(self.ucat):
                # prepare model RDMs for this category
                c_mats = []
                for mat in models[model_type]:
                    # get model similarity for this category
                    c_mat = mat[np.ix_(cat==c, cat==c)]

                    # replace any missing data with the median within
                    # this category
                    vec_ind = np.triu_indices(c_mat.shape[0], 1)
                    c_mat[np.isnan(c_mat)] = np.nanmedian(c_mat[vec_ind])
                    c_mats.append(c_mat)

                # prepare design matrices for all xval folds
                f_list = []
                for fold in self.xval_ind[i]:
                    r_list = []
                    
                    # prep design matrices for train and test blocks
                    f_dict = {}
                    for block, x_ind in fold.iteritems():
                        n_item = len(x_ind)
                        n_feat = (n_item ** 2 - n_item) / 2
                        design = np.ones((n_feat, n_model + 1))
                        for j, mat in enumerate(c_mats):
                            design[:,j+1] = squareform(mat[np.ix_(x_ind,x_ind)])
                        f_dict[block] = design
                    r_list.append(f_dict)

                    # same for randomized models
                    for r_ind in rand_ind[i]:
                        f_dict = {}
                        for block, x_ind in fold.iteritems():
                            n_item = len(x_ind)
                            n_feat = (n_item ** 2 - n_item) / 2
                            design = np.ones((n_feat, n_model + 1))
                            for j, mat in enumerate(c_mats):
                                # scramble, and then take the block
                                # (block indices stay still relative
                                # to scrambled items)
                                r_mat = mat[np.ix_(r_ind,r_ind)]
                                b_mat = r_mat[np.ix_(x_ind,x_ind)]
                                design[:,j+1] = squareform(b_mat)
                            f_dict[block] = design
                        r_list.append(f_dict)
                    f_list.append(r_list)
                c_list.append(f_list)
            self.designs[model_type] = c_list
            
    def __call__(self, dataset):
        
        # calculate data RDM within each included category
        y_list = []
        n_cat = len(self.ucat)
        data_problem = False
        for i in range(n_cat):
            # slower to deal with squareformed data, but much easier
            # to handle xval with item indexing instead of matrix
            # element indexing
            data_rdm = squareform(pdist(dataset.samples[self.cat==self.ucat[i],:],
                                        'correlation'))
            if np.any(np.isnan(data_rdm)):
                data_problem = True
                break
            y_list.append(data_rdm)
        
        # if any nans encountered in the data rdm, just output zeros
        if data_problem:
            if self.output == 'full':
                return tuple(np.zeros(self.n_perm+1))
            else:
                return np.float64(0)

        # calculate correlation between the data and the best-fitting
        # model-based RDM, separately for each model type
        rho = {}
        for model_type in self.designs.keys():
            n_model = len(self.designs[model_type])
            rho_type = np.zeros((n_cat, self.n_xval, self.n_perm+1))
            for i in range(n_cat):
                # data to fit
                y = y_list[i]
                for j in range(self.n_xval):
                    # get train and test data for this category and xval fold
                    data = {}
                    for block, x_ind in self.xval_ind[i][j].iteritems():
                        data[block] = squareform(y[np.ix_(x_ind,x_ind)])

                    # train and test the model for each permutation
                    for k in range(self.n_perm+1):
                        # train
                        model = self.designs[model_type][i][j][k]
                        if self.fit == 'nnls':
                            x = optim.nnls(model['train'], data['train'])[0]
                        else:
                            x = linalg.lstsq(model['train'], data['train'])[0]

                        # test
                        yhat = np.dot(model['test'], x);
                        r = stats.spearmanr(yhat,data['test'])[0]
                        if not np.isnan(r):
                            rho_type[i,j,k] = r
            # average over xval fold
            rho[model_type] = np.mean(rho_type, 1)

        # statistic is a difference of model types, averaged over
        # categories
        stat_perm = np.mean(rho[self.type1] - rho[self.type2], 0)
        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)
