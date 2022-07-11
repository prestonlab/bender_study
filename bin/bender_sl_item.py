#!/usr/bin/env python
#
# Searchlight of category reactivation.

from bender_study.subjutil import *

parser = SubjParser(include_log=False)
parser.add_argument("mask", help="name of mask file")
parser.add_argument("resname", help="name of results directory")
parser.add_argument(
    "measure",
    choices=["item_react", "item_suppress", "item_react_sme", "item_suppress_sme"],
    help="measure of reactivation",
)
parser.add_argument(
    "--output", "-o", default="full", help="output type ('full' or 'z')"
)
parser.add_argument(
    "--cat", "-c", default=None, help="category to include (face,scene,[both])"
)
parser.add_argument(
    "--division",
    "-d",
    default="cat",
    choices=["cat", "subcat"],
    help="category labels to use",
)
parser.add_argument(
    "--rec",
    "-e",
    default=None,
    type=int,
    help="subsequent recall condition to include (0,1)",
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
from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.fx import mean_group_sample
from mvpa2.mappers.zscore import zscore
from mvpa2.measures.searchlight import sphere_searchlight

from bender_study import patio
from bender_study import bender
from bender_study import rsa


bp = bender.BenderPath(args.subject, args.study_dir)
mask_file = bp.image_path("anatomy", "bbreg", "data", args.mask)

# load pre-exposure data, with one pattern per A item
print("Loading pre-exposure data with %s mask..." % args.mask)
ds = bp.beta_dataset("a_prex", mask_file, suffix=args.suffix)

# filter to get just included trials
include = np.ones(ds.shape[0], dtype=bool)
if args.cat is not None:
    include = np.logical_and(include, ds.sa.category == args.cat)
if args.rec is not None:
    test = bp.read_period("ac_test", 5.5).sort("group")
    inc_groups = test.filter(correct=args.rec).array("group")
    rec_include = np.array([g in inc_groups for g in ds.sa.group])
    include = np.logical_and(include, rec_include)
ds = ds[include, :]

# normalize within run
zscore(ds, chunks_attr="chunks")

# get average pattern for each item
m = mean_group_sample(["group"])
dsm = ds.get_mapped(m)
ds = dsm

# load study data, with one pattern per BC presentation, sorted by
# group number
print("Loading BC study beta series with %s mask..." % args.mask)
ds2 = bp.beta_dataset("bcxy_study", mask_file, suffix=args.suffix)
if args.division == "cat":
    ds2.sa["targets"] = ds2.sa.category
elif args.division == "subcat":
    ds2.sa["targets"] = ds2.sa.subcategory
else:
    raise ValueError("Unknown division: {}".format(args.division))
ind = np.argsort(ds2.sa.group)
ds2 = ds2[ind, :]
ds2 = ds2[ds2.sa.cond < 5, :]

# filter to get just included trials
include = np.ones(ds2.shape[0], dtype=bool)
if args.cat is not None:
    # include only items in the specified category
    include = np.logical_and(include, ds2.sa.category == args.cat)
if args.rec is not None:
    # filter to get only subsequently recalled/forgotten items
    test = bp.read_period("ac_test", 5.5).sort("group")
    include = np.logical_and(include, test.array("correct") == args.rec)
ds2 = ds2[include, :]
zscore(ds2, chunks_attr="chunks")

# get category/subcategory labels (some analyses can use either)
if args.division == "cat":
    category = ds2.sa.category
elif args.division == "subcat":
    category = ds2.sa.cond
else:
    raise ValueError("Division must be 'cat' or 'subcat'.")

# combine datasets with special handling of attributes
ds_full = patio.combine_datasets((ds, ds2))
del ds

print("Initializing searchlight...")

if args.measure in ["item_react", "item_suppress"]:
    # within-item vs. within-(sub)category similarity
    if args.measure == "item_react":
        stat = "react"
    else:
        stat = "suppress"
    m = rsa.SimItemReact(
        ds_full.sa.dataset, category, args.n_perm, stat=stat, output=args.output
    )
elif args.measure in [
    "item_react_sme",
    "item_suppress_sme",
]:
    # label of whether each included item was correctly remembered
    test = bp.read_period("ac_test", 5.5).sort("group")
    correct = test.array("correct")[include]

    # item reactivation as a function of sequent memory or RT
    if "react" in args.measure:
        stat = "react_sme"
    else:
        stat = "suppress_sme"
    m = rsa.SimItemReactSME(
        ds_full.sa.dataset,
        category,
        correct,
        args.n_perm,
        stat=stat,
        output=args.output,
    )

print("Running searchlight...")
sl = sphere_searchlight(m, radius=args.radius, nproc=args.n_proc)
slmap = sl(ds_full)

print("Saving maps...")
res_dir = bp.path("mvpa", args.resname)
if not os.path.exists(res_dir):
    os.makedirs(res_dir)

if args.output == "z":
    filepath = os.path.join(res_dir, "zstat.nii.gz")
    nifti = map2nifti(ds2, slmap)
    nifti.to_filename(filepath)
else:
    for i in range(len(slmap)):
        filepath = os.path.join(res_dir, "perm%03d.nii.gz" % i)
        nifti = map2nifti(ds2, slmap[i])
        nifti.to_filename(filepath)
