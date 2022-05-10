#!/usr/bin/env python

from bender_study.subjutil import *
parser = SubjParser(include_log=False)
parser.add_argument('mask', help="name of mask file")
parser.add_argument('resname', help="name of results directory")
parser.add_argument('measure', help="measure of reactivation")
parser.add_argument('--output', '-o', default='full',
                    help="output type ('full' or 'z')")
parser.add_argument('--cat', '-c', default=None,
                    help="category to include (face,scene,[both])")
parser.add_argument('--division', '-d', default='cat',
                    help="category labels to use ([cat],subcat)")
parser.add_argument('--rec', '-e', default=None, type=int,
                    help="subsequent recall condition to include (0,1)")
parser.add_argument('--suffix', '-s', help="suffix for beta images",
                    default='_stim')
parser.add_argument('--radius', '-r', type=int, default=3,
                    help="searchlight radius")
parser.add_argument('--n-perm', '-p', type=int, default=100,
                    help="number of permutations to run")
parser.add_argument('--n-proc', '-n', type=int, default=None,
                    help="processes for searchlight")
args = parser.parse_args()

import os
import numpy as np
from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.fx import mean_group_sample
from mvpa2.mappers.zscore import zscore
from mvpa2.measures.searchlight import sphere_searchlight
from mvpa2.clfs.svm import LinearCSVMC

from bender_study import patio
from bender_study import bender
from bender_study import traintest
from bender_study import sim_react
from bender_study import sim_rdm_react


bp = bender.BenderPath(args.subject, args.study_dir)
mask_file = bp.image_path('anatomy', 'bbreg', 'data', args.mask)

if 'svm' in args.measure:
    # load localizer data for training a classifier
    print("Loading localizer data with %s mask..." % args.mask)
    ds = bp.block_dataset('c_loc', 2., mask=mask_file,
                          shift=4., suffix='_hpfsm')
    
    if args.division == 'cat':
        ds.sa['targets'] = ds.sa.category
    elif args.division == 'subcat':
        ds.sa['targets'] = ds.sa.subcategory
    else:
        raise ValueError('Unknown division: {}'.format(args.division))
        
    if args.cat is not None:
        # include only the specified categories
        include = ds.sa.category == args.cat
    else:
        # include only faces and scenes
        include = np.array([x[0] in 'fs' for x in ds.sa.category])
    ds = ds[include,:]
    zscore(ds, chunks_attr='chunks')
else:
    # load pre-exposure data, with one pattern per A item
    print("Loading pre-exposure data with %s mask..." % args.mask)
    ds = bp.beta_dataset('a_prex', mask_file, suffix=args.suffix)

    # filter to get just included trials
    include = np.ones(ds.shape[0], dtype=bool)
    if args.cat is not None:
        include = np.logical_and(include, ds.sa.category == args.cat)
    if args.rec is not None:
        test = bp.read_period('ac_test', 5.5).sort('group')
        inc_groups = test.filter(correct=args.rec).array('group')
        rec_include = np.array([g in inc_groups for g in ds.sa.group])
        include = np.logical_and(include, rec_include)
    ds = ds[include,:]
    
    zscore(ds, chunks_attr='chunks')
    m = mean_group_sample(['group'])
    dsm = ds.get_mapped(m)
    ds = dsm
    ds.sa['targets'] = ds.sa.condition

# load study data, with one pattern per BC presentation, sorted by
# group number
print("Loading BC study beta series with %s mask..." % args.mask)
ds2 = bp.beta_dataset('bcxy_study', mask_file, suffix=args.suffix)
if args.division == 'cat':
    ds2.sa['targets'] = ds2.sa.category
elif args.division == 'subcat':
    ds2.sa['targets'] = ds2.sa.subcategory
else:
    raise ValueError('Unknown division: {}'.format(args.division))
ind = np.argsort(ds2.sa.group)
ds2 = ds2[ind,:]
ds2 = ds2[ds2.sa.cond < 5,:]

# filter to get just included trials
include = np.ones(ds2.shape[0], dtype=bool)
if args.cat is not None:
    # include only items in the specified category
    include = np.logical_and(include, ds2.sa.category == args.cat)
if args.rec is not None:
    # filter to get only subquently recalled/forgotten items
    test = bp.read_period('ac_test', 5.5).sort('group')
    include = np.logical_and(include, test.array('correct')==args.rec)
ds2 = ds2[include,:]
zscore(ds2, chunks_attr='chunks')

# get category/subcategory labels (some analyses can use either)
if args.division == 'cat':
    category = ds2.sa.category
elif args.division == 'subcat':
    category = ds2.sa.cond

if 'sme' in args.measure:
    # label of whether each included item was correctly remembered
    test = bp.read_period('ac_test', 5.5).sort('group')
    correct = test.array('correct')[include]
elif 'rt' in args.measure:
    # information about the test period, sorted by item number
    test = bp.read_period('ac_test', 5.5).sort('group')
    correct = test.array('correct')[include]
    rt = test.array('rt')[include]

    # get just subsequently correct trials
    ds = ds[correct==1,:]
    ds2 = ds2[correct==1,:]
    category = category[correct==1]
    rt_correct = rt[correct==1]
    
    # for correct trials in each (sub)category, sort into fast and slow
    # responses
    fast = np.zeros((len(rt_correct),))
    for c in np.unique(category):
        m = np.median(rt_correct[category==c])
        fast[np.logical_and(category==c, rt_correct<=m)] = 1

    # for convenience, use same variable name for sme and rt measures
    correct = fast
    
# combine datasets with special handling of attributes
ds_full = patio.combine_datasets((ds, ds2))
del ds

print("Initializing searchlight...")

if args.measure == 'svm':
    # AUC for category decoding
    clf = LinearCSVMC(probability=1, enable_ca=['probabilities'])
    m = traintest.TrainTest(clf, 'dataset', ds2.sa.chunks, n_perm=args.n_perm,
                  output=args.output)
elif args.measure == 'svm_sme':
    # AUC for category coding, by subsequent memory
    clf = LinearCSVMC(probability=1, enable_ca=['probabilities'])
    m = traintest.TrainTestSME(clf, ds_full.sa.dataset, ds2.category,
                     correct, n_perm=args.n_perm, output=args.output)

elif args.measure == 'cat':
    # within vs. between-category similarity
    m = sim_react.SimCatReact(ds_full.sa.dataset, ds2.sa.category,
                    n_perm=args.n_perm, output=args.output)
elif args.measure in ['cat_sme', 'cat_rt']:
    # category separation as a function of subsequent memory or RT
    m = sim_react.SimCatReactSME(ds_full.sa.dataset, ds2.sa.category,
                       correct, n_perm=args.n_perm, output=args.output)

elif args.measure == 'subcat':
    # within vs. between subcategory similarity (within category)
    m = sim_react.SimSubcatReact(ds_full.sa.dataset, ds2.sa.category,
                       ds2.sa.cond, n_perm=args.n_perm, output=args.output)
elif args.measure in ['subcat_sme', 'subcat_rt']:
    # subcategory separation as a function of subsequent memory or RT
    m = sim_react.SimSubcatReactSME(ds_full.sa.dataset, ds2.sa.category,
                          ds2.sa.cond, correct, n_perm=args.n_perm,
                          output=args.output)

elif args.measure in ['item_react', 'item_suppress']:
    # within-item vs. within-(sub)category similarity
    if args.measure == 'item_react':
        stat = 'react'
    else:
        stat = 'suppress'
    m = sim_react.SimItemReact(ds_full.sa.dataset, category, args.n_perm,
                     stat=stat, output=args.output)
elif args.measure in ['item_react_sme', 'item_suppress_sme',
                      'item_react_rt', 'item_suppress_rt']:
    # item reactivation as a function of sequent memory or RT
    if 'react' in args.measure:
        stat = 'react_sme'
    else:
        stat = 'suppress_sme'
    m = sim_react.SimItemReactSME(ds_full.sa.dataset, category, correct, args.n_perm,
                        stat=stat, output=args.output)

elif args.measure == 'rdm':
    # correlation between prex and study RDMs
    m = sim_rdm_react.SimRDMReact('dataset', ds2.shape[0], args.n_perm)
    
print("Running searchlight...")
sl = sphere_searchlight(m, radius=args.radius, nproc=args.n_proc)
slmap = sl(ds_full)

print("Saving maps...")
res_dir = bp.path('mvpa', args.resname)
if not os.path.exists(res_dir):
    os.makedirs(res_dir)

if args.output == 'z':
    filepath = os.path.join(res_dir, 'zstat.nii.gz')
    nifti = map2nifti(ds2, slmap)
    nifti.to_filename(filepath)
else:
    for i in range(len(slmap)):
        filepath = os.path.join(res_dir, 'perm%03d.nii.gz' % i)
        nifti = map2nifti(ds2, slmap[i])
        nifti.to_filename(filepath)
