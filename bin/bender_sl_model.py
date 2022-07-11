#!/usr/bin/env python
#
# Searchlight contrasting RDM correlations for study data.

from bender_study.subjutil import SubjParser

parser = SubjParser(include_log=False)
parser.add_argument("mask", help="name of mask for searchlight")
parser.add_argument("model_types", help="model contrast (e.g. bx-cy)")
parser.add_argument("resname", help="name of results directory")
parser.add_argument(
    "--output", "-o", default="full", help="output type ('full' or 'z')"
)
parser.add_argument(
    "--cat", "-c", default=None, help="category to include (face,scene,[both])"
)
parser.add_argument("--suffix", "-s", help="suffix for beta images", default="_stim")
parser.add_argument("--radius", "-r", type=int, default=3, help="searchlight radius")
parser.add_argument(
    "--n-perm", "-p", type=int, default=100, help="number of permutations to run"
)
parser.add_argument(
    "--n-proc", "-n", type=int, default=None, help="processes for searchlight"
)
args = parser.parse_args()

import os
import numpy as np
from scipy.io import loadmat

from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.zscore import zscore
from mvpa2.measures.searchlight import sphere_searchlight

from bender_study import bender
from bender_study import rsa

bp = bender.BenderPath(args.subject, args.study_dir)
mask_file = bp.image_path("anatomy", "bbreg", "data", args.mask)

if not os.path.exists(mask_file):
    raise IOError("Mask does not exist: {}".format(mask_file))

model_dir = os.path.join(bp.study_dir, "batch", "models3", "trial")
if not os.path.exists(model_dir):
    raise IOError("Model directory does not exist: {}".format(model_dir))

# load all models
models = []
model_types = args.model_types.split("-")
for model_type in model_types:
    fname = "{}_wiki_w2v_{}.mat".format(bp.subject, model_type)
    mat_file = os.path.join(model_dir, fname)
    if not os.path.exists(mat_file):
        raise IOError("Model file does not exist: {}".format(mat_file))
    models.append(loadmat(mat_file)["rdm"])

# load data, sort, z-score
ds = bp.beta_dataset("bcxy_study", mask_file, suffix=args.suffix)
ind = np.argsort(ds.sa.group)
ds = ds[ind, :]
zscore(ds, chunks_attr="chunks")

# get BC trials
inc = ds.sa.cond < 5
n_trial = ds.shape[0]
ds = ds[inc, :]
for i in range(len(models)):
    if models[i].shape[0] == n_trial:
        models[i] = models[i][np.ix_(inc, inc)]

# add AC test accuracy
test = bp.read_period("ac_test", 5.5)
correct = test.sort("group").array("correct")

# sorted condition values are used to set the interaction, so set
# correct items to the lower value in the conditions
# vector.

# (correct: model 1 - model 2) - (incorrect: model 1 - model 2)
cond = np.zeros(ds.shape[0])
cond[correct == 1] = 1
cond[correct == 0] = 2

# re-sort the dataset and models to be in condition order
ind = np.argsort(cond)
ds = ds[ind, :]
cond = cond[ind]
for i in range(len(models)):
    models[i] = models[i][np.ix_(ind, ind)]

m = rsa.SimModelCond(cond, models, args.n_perm, output=args.output)
sl = sphere_searchlight(m, radius=args.radius, nproc=args.n_proc)
slmap = sl(ds)

print("Saving maps...")
res_dir = bp.path("rsa", args.resname)
if not os.path.exists(res_dir):
    os.makedirs(res_dir)

if args.output == "z":
    filepath = os.path.join(res_dir, "zstat.nii.gz")
    nifti = map2nifti(ds, slmap)
    nifti.to_filename(filepath)
else:
    for i in range(len(slmap)):
        filepath = os.path.join(res_dir, "perm%03d.nii.gz" % i)
        nifti = map2nifti(ds, slmap[i])
        nifti.to_filename(filepath)
