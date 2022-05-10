"""Train and test a classifier on different parts of a dataset"""

import random
import numpy as np
import scipy.stats as stats
from mvpa2.measures.base import Measure
from sklearn.metrics import roc_curve, auc
from sim_react import perm_within, perm_z

class TrainTest(Measure):

    def __init__(self, clf, split_attr, run, n_perm=1000, output='full'):

        Measure.__init__(self)
        self.clf = clf
        self.split_attr = split_attr
        self.n_perm = n_perm
        self.output = output
        
        # scramble the test set within run
        rand_ind = perm_within(run, n_perm)
        n_item = len(run)
        rand_ind.insert(0, np.arange(n_item))
        self.rand_ind = rand_ind

    def __call__(self, dataset):
        
        # split into the chunks we are comparing
        split = dataset.sa[self.split_attr].value
        usplit = np.unique(split)
        if len(usplit) != 2:
            raise ValueError('Split attribute must have two unique values.')
        ds1 = dataset[split == usplit[0],:]
        ds2 = dataset[split == usplit[1],:]

        # train the classifier
        clf = self.clf
        clf.train(ds1)

        # test and get evidence for each category
        pred = clf.predict(ds2.samples)
        evid = np.array([p[1][0]-p[1][1] for p in clf.ca.probabilities])

        stat_perm = np.zeros(self.n_perm + 1)
        for i, ind in enumerate(self.rand_ind):
            # calculate AUC given the (actual or scrambled) category labels
            fpr, tpr, thresh = roc_curve(ds2.targets[ind], evid,
                                         pos_label=clf.ca.trained_targets[0])
            stat_perm[i] = auc(fpr, tpr)

        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)

class TrainTestSME(Measure):

    def __init__(self, clf, split, cat, correct, n_perm=100, output='full'):
        Measure.__init__(self)
        self.clf = clf
        self.n_perm = n_perm
        self.split_ind = tuple([np.nonzero(split==u)[0]
                                for u in np.unique(split)])
        self.output = output
        
        # scramble items within category
        rand_ind = perm_within(cat, n_perm)
        n_item = len(cat)
        rand_ind.insert(0, np.arange(n_item))

        # indices to get recalled and forgotten items
        self.r_ind = []
        self.f_ind = []
        for ind in rand_ind:
            r_correct = correct[ind]
            self.r_ind.append(np.nonzero(r_correct==1)[0])
            self.f_ind.append(np.nonzero(r_correct==0)[0])

    def __call__(self, dataset):
        ds1 = dataset[self.split_ind[0]]
        ds2 = dataset[self.split_ind[1]]

        clf = self.clf
        clf.train(ds1)
        pred = clf.predict(ds2.samples)
        evid = np.array([p[1][0]-p[1][1] for p in clf.ca.probabilities])

        stat_perm = np.zeros((self.n_perm + 1,))
        for i in range(self.n_perm + 1):
            # null hypothesis: there is no difference in category
            # discriminability between recalled and forgotten
            # items. This permutation test preserves the category
            # labels but scrambles the subsequent memory labels within
            # category
            fpr, tpr, thresh = roc_curve(ds2.targets[self.r_ind[i]],
                                         evid[self.r_ind[i]],
                                         pos_label=clf.ca.trained_targets[0])
            r_auc = auc(fpr, tpr)

            fpr, tpr, thresh = roc_curve(ds2.targets[self.f_ind[i]],
                                         evid[self.f_ind[i]],
                                         pos_label=clf.ca.trained_targets[0])
            f_auc = auc(fpr, tpr)
            stat_perm[i] = r_auc - f_auc

        if self.output == 'full':
            return tuple(stat_perm)
        else:
            return perm_z(stat_perm)
