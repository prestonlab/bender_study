#!/usr/bin/env python
#
# Export a study-phase RDM for an anatomical mask.

from bender_study.subjutil import *

parser = SubjParser(include_log=False)
parser.add_argument("mask", help="name of mask file")
parser.add_argument("--suffix", "-s", help="suffix for beta images", default="_stim")
args = parser.parse_args()

import os

import numpy as np
from scipy.spatial.distance import pdist, squareform
from mvpa2.mappers.zscore import zscore

from bender_study import bender

bp = bender.BenderPath(args.subject, args.study_dir)
mask_file = bp.image_path("anatomy", "bbreg", "data", args.mask)

print("Loading BC study beta series with {} mask...".format(args.mask))
ds = bp.beta_dataset("bcxy_study", mask_file, suffix=args.suffix)
zscore(ds, chunks_attr="chunks")
ind = np.argsort(ds.sa.group)
ds = ds[ind, :]

rdm = squareform(pdist(ds.samples, "correlation"))

# directory for this ROI
model = "study" + args.suffix
rdm_dir = os.path.join(args.study_dir, "batch", "glm", model, "rdm", args.mask)
if not os.path.exists(rdm_dir):
    os.makedirs(rdm_dir)
rdm_file = os.path.join(rdm_dir, "{}_rdm.txt".format(args.subject))

print("Saving similarity values to: {}".format(rdm_file))
np.savetxt(rdm_file, rdm)
