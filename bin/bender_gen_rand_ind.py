#!/usr/bin/env python
#
# Generate random indices for group volume threshold analysis.

import sys
import random

n_subj = int(sys.argv[1])
n_subj_perm = int(sys.argv[2])
n_perm = int(sys.argv[3])
subj_ind = range(n_subj)
for i in range(n_perm):
    line = ""
    for j in range(n_subj):
        # randomly sample subject and permutation within subject with
        # replacement
        subj = random.randint(0, n_subj - 1)
        perm = random.randint(0, n_subj_perm)

        # get index within the tiled dataset. e.g., for (subject,
        # permutation) tuples with two subjects and three
        # permutations/subject: [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)]
        ind = subj * (n_subj_perm + 1) + perm
        if j == (n_subj - 1):
            line += "{}".format(ind)
        else:
            line += "{} ".format(ind)
    print(line)
