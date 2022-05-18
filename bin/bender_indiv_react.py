#!/usr/bin/env python
#
# Export an RDM comparing pre-exposure and study patterns.

from bender_study.subjutil import *

parser = SubjParser(include_log=False)
parser.add_argument("mask", help="name of mask file")
parser.add_argument("--suffix", "-s", help="suffix for beta images", default="_stim")
args = parser.parse_args()

import os
import numpy as np
from scipy.spatial.distance import cdist
from mvpa2.mappers.fx import mean_group_sample
from mvpa2.mappers.zscore import zscore

from bender_study import bender

bp = bender.BenderPath(args.subject, args.study_dir)
mask_file = bp.image_path("anatomy", "bbreg", "data", args.mask)

print("Loading pre-exposure data with %s mask..." % args.mask)
ds = bp.beta_dataset("a_prex", mask_file, suffix=args.suffix)
zscore(ds, chunks_attr="chunks")
m = mean_group_sample(["group"])
dsm = ds.get_mapped(m)
ds = dsm

print("Loading BC study beta series with %s mask..." % args.mask)
ds2 = bp.beta_dataset("bcxy_study", mask_file, suffix=args.suffix)
zscore(ds2, chunks_attr="chunks")
ind = np.argsort(ds2.sa.group)
ds2 = ds2[ind, :]
ds2 = ds2[ds2.sa.cond < 5, :]

rdm = cdist(ds.samples, ds2.samples, "correlation")

# directory for this ROI
model = "study" + args.suffix
rdm_dir = os.path.join(args.study_dir, "batch", "glm", model, "react", args.mask)
if not os.path.exists(rdm_dir):
    os.makedirs(rdm_dir)
rdm_file = os.path.join(rdm_dir, "{}_rdm.txt".format(args.subject))

print("Saving similarity values to: {}".format(rdm_file))
np.savetxt(rdm_file, rdm)
