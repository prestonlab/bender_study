#!/usr/bin/env python
#
# Searchlight comparing a model RDM to neural dissimilarity.

from bender_study.subjutil import SubjParser
from bender_study import sim_model

parser = SubjParser(include_log=False)
parser.add_argument("phase", help="experiment phase (prex or study)")
parser.add_argument("rdm", help="path to RDM to use")
parser.add_argument("mask", help="name of mask for searchlight")
parser.add_argument("resname", help="name of results directory")
parser.add_argument(
    "--feature-mask", default=None, help="name of mask of voxels to include as features"
)
parser.add_argument(
    "--output", "-o", default="full", help="output type (['full'],'z','roi')"
)
parser.add_argument(
    "--cat", "-c", default=None, help="category to include (face,scene,[both])"
)
parser.add_argument(
    "--division",
    "-d",
    default="cat",
    help="individual part of RDM to fit ([cat],subcat,none)",
)
parser.add_argument(
    "--rec",
    "-e",
    default=None,
    type=int,
    help="subsequent recall condition to include (0,1)",
)
parser.add_argument(
    "--measure", "-m", default="corr", help="measure to calculate ([corr],sme)"
)
parser.add_argument(
    "--partial", default=None, type=str, help="models to partial out (models3)"
)
parser.add_argument(
    "--fit",
    "-f",
    default="ls",
    type=str,
    help="regression type for partial correlation (ls,nnls)",
)
parser.add_argument("--suffix", "-s", help="suffix for beta images", default="_stim")
parser.add_argument("--radius", "-r", type=int, default=3, help="searchlight radius")
parser.add_argument(
    "--n-perm", "-p", type=int, default=100, help="number of permutations to run"
)
parser.add_argument(
    "--ind-file", "-i", default=None, help="MAT-file with random indices"
)
parser.add_argument(
    "--n-proc", "-n", type=int, default=None, help="processes for searchlight"
)
args = parser.parse_args()

import os
import numpy as np
import scipy.stats as stats
from scipy.io import loadmat
from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.fx import mean_group_sample
from mvpa2.mappers.zscore import zscore
from mvpa2.measures.searchlight import sphere_searchlight

from bender_study import bender

bp = bender.BenderPath(args.subject, args.study_dir)
mask_file = bp.image_path("anatomy", "bbreg", "data", args.mask)

if not os.path.exists(mask_file):
    raise IOError("Mask does not exist: {}".format(mask_file))

if not os.path.exists(args.rdm):
    raise IOError("RDM file does not exist: {}".format(args.rdm))

# load model
filename = os.path.splitext(os.path.basename(args.rdm))[0]

if os.path.splitext(args.rdm)[1] == ".mat":
    s = loadmat(args.rdm)
    if "rdm" in s:
        rdm = s["rdm"]
    else:
        rdm = s[filename]
else:
    rdm = np.loadtxt(args.rdm)

# load control models
if args.partial is not None:
    if args.partial == "models3":
        # standard models, version 3
        model_dir = os.path.join(bp.study_dir, "batch", "models3")
        model_name = filename.split("mat_")[1]
        control_models = ["hmax", "subcat", "wiki_w2v"]
    elif args.partial == "models3g":
        # standard models + geo model
        model_dir = os.path.join(bp.study_dir, "batch", "models3")
        model_name = filename.split("mat_")[1]
        control_models = ["hmax", "subcat", "wiki_w2v", "geo"]
    elif args.partial == "models3e":
        # standard models + geo and ego models
        model_dir = os.path.join(bp.study_dir, "batch", "models3")
        model_name = filename.split("mat_")[1]
        control_models = ["hmax", "subcat", "wiki_w2v", "geo", "ego"]
    elif args.partial == "models3c":
        model_dir = os.path.join(bp.study_dir, "batch", "models3")
        model_name = filename.split("mat_")[1]
        control_models = ["hmax", "cat", "subcat", "wiki_w2v"]
    elif args.partial == "models3ns":
        # standard models without subcategory
        model_dir = os.path.join(bp.study_dir, "batch", "models3")
        model_name = filename.split("mat_")[1]
        control_models = ["hmax", "wiki_w2v"]

    # load models other than the current one
    include = [m for m in control_models if m != model_name]
    control_rdms = []
    for model in include:
        cname = "mat_" + model
        s = loadmat(os.path.join(model_dir, cname + ".mat"))
        if "rdm" in s:
            crdm = s["rdm"]
        else:
            crdm = s[cname]
        control_rdms.append(crdm)
else:
    control_rdms = None

# load data
if args.phase == "prex":
    if args.feature_mask is not None:
        feature_file = bp.image_path("anatomy", "bbreg", "data", args.feature_mask)
        print("Using features within: {}".format(feature_file))
        ds = bp.beta_dataset(
            "a_prex", mask_file, suffix=args.suffix, include=feature_file
        )
        ds.fa.include = ds.fa.include.astype(bool)
    else:
        ds = bp.beta_dataset("a_prex", mask_file, suffix=args.suffix)
        ds.fa["include"] = np.ones(ds.shape[1], dtype=bool)
    zscore(ds, chunks_attr="chunks")
    m = mean_group_sample(["group"])
    dsm = ds.get_mapped(m)
    ds = dsm

elif args.phase == "study":
    ds = bp.beta_dataset("bcxy_study", mask_file, suffix=args.suffix)
    ind = np.argsort(ds.sa.group)
    ds = ds[ind, :]
    zscore(ds, chunks_attr="chunks")
    if args.measure in ["sme", "sme_cat", "corr"]:
        # filter dataset to get just BC trials
        inc = ds.sa.cond < 5
        n_trial = ds.shape[0]
        ds = ds[inc, :]

        # also filter rdm if necessary
        if rdm.shape[0] == n_trial:
            rdm = rdm[np.ix_(inc, inc)]

# filter by category
include = np.ones(ds.shape[0], dtype=bool)
category = ds.sa.category
if args.cat is not None:
    include = np.logical_and(include, category == args.cat)

# remove items with missing values in model
notmissing = np.sum(np.isnan(rdm), 0) < (len(rdm) - 1)
include = np.logical_and(include, notmissing)

# filter out any excluded items
ds = ds[include, :]
rdm = rdm[np.ix_(include, include)]
if control_rdms is not None:
    for i in range(len(control_rdms)):
        control_rdms[i] = control_rdms[i][np.ix_(include, include)]

if args.ind_file is not None:
    rand_ind = loadmat(args.ind_file)["rand_ind"]
    rand_ind = rand_ind[: args.n_perm, :]
    if args.cat is None:
        raise IOError(
            "If loading random indices from file, cannot use both categories."
        )

    cat_inc = include[category == args.cat]
    rand_filt = []
    for ind in rand_ind:
        # reorder include vector to match random indices, then filter
        ind_filt = ind[cat_inc[ind]]

        # get indices that match the filtered matrix
        subind = np.asarray(stats.rankdata(ind_filt) - 1, dtype=int)
        rand_filt.append(subind)
    rand_ind = rand_filt
else:
    rand_ind = None

if args.division == "cat":
    division = ds.sa.category
elif args.division == "subcat":
    division = ds.sa.cond
elif args.division == "none":
    division = np.ones(ds.sa.cond.shape)
else:
    raise ValueError("Unknown division type: {}".format(args.division))

if args.measure == "corr":
    if args.partial is not None:
        m = sim_model.SimRDMPartial2(
            division, rdm, control_rdms, args.n_perm, fit=args.fit, output=args.output
        )
    else:
        m = sim_model.SimRDM(
            division, rdm, args.n_perm, output=args.output, rand_ind=rand_ind
        )
elif args.measure == "sme":
    test = bp.read_period("ac_test", 5.5)
    correct = test.sort("group").array("correct")
    m = sim_model.SimRDMSME(division, correct[include], rdm, args.n_perm, args.output)
elif args.measure == "sme_cat":
    test = bp.read_period("ac_test", 5.5)
    correct = test.sort("group").array("correct")
    m = sim_model.SimRDMCatSME(
        division, correct[include], rdm, args.n_perm, args.output
    )
elif args.measure == "bcxy":
    cond = np.ones(ds.shape[0])
    cond[np.logical_or(ds.sa.cond == 1, ds.sa.cond == 2)] = 1
    cond[np.logical_or(ds.sa.cond == 3, ds.sa.cond == 4)] = 2
    cond[ds.sa.cond == 5] = 3
    contrast = np.array([1, 1, -2])
    m = sim_model.SimRDMContrast(cond, contrast, rdm, args.n_perm, output=args.output)

if args.output == "roi":
    print("Calculating statistic...")
    stat = m(ds)[0]

    res_dir = os.path.join(bp.study_dir, "batch", "rsa", args.resname, "roi", args.mask)
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)
    filepath = os.path.join(res_dir, "stat_{}.txt".format(bp.subject))

    print("Saving statistic in: {}".format(filepath))
    np.savetxt(filepath, stat)
else:
    print("Running searchlight...")
    sl = sphere_searchlight(m, radius=args.radius, nproc=args.n_proc)
    slmap, incmap = sl(ds)

    res_dir = bp.path("rsa", args.resname)
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

    print("Saving maps in: {}".format(res_dir))
    incfile = os.path.join(res_dir, "included.nii.gz")
    nifti = map2nifti(ds, incmap)
    nifti.to_filename(incfile)

    if args.output == "z":
        filepath = os.path.join(res_dir, "zstat.nii.gz")
        nifti = map2nifti(ds, slmap)
        nifti.to_filename(filepath)
    else:
        for i in range(len(slmap)):
            filepath = os.path.join(res_dir, "perm%03d.nii.gz" % i)
            nifti = map2nifti(ds, slmap[i])
            nifti.to_filename(filepath)
