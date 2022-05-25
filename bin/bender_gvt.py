#!/usr/bin/env python
#
# Run group volume threshold analysis.


from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("inputds", help="input dataset hdf5 file")
parser.add_argument("stat", help="output file for actual stat")
parser.add_argument("thresh", help="output file for voxelwise threshold image")
parser.add_argument(
    "-a",
    "--alpha",
    dest="feature_thresh_prob",
    help="featurewise threshold probability",
    type=float,
    default=0.001,
)
parser.add_argument(
    "-b", "--n-blocks", default=1, type=int, help="number of compute blocks"
)
parser.add_argument("-p", "--n-proc", default=1, type=int, help="number of processes")
parser.add_argument(
    "-n", "--n-bootstrap", default=100000, type=int, help="number of bootstraps"
)
parser.add_argument(
    "-m", "--mmap-file", default=None, help="path to memory-mapped file"
)
parser.add_argument(
    "-r", "--rand-ind", default=None, help="path to file with random indices"
)
args = parser.parse_args()

import numpy as np
import nibabel as nib

from mvpa2.base import hdf5
from mvpa2.datasets.mri import map2nifti
from bender_study import gvt

print("Loading data and permutations...")
ds = hdf5.h5load(args.inputds)

if ds.shape[1] == 0:
    # mask is empty for this image; save out zeros
    print("Mask is empty. Saving zeros...")
    nib.save(map2nifti(ds, np.zeros(ds.shape[1])), args.stat)
    nib.save(map2nifti(ds, np.zeros(ds.shape[1])), args.thresh)
else:
    if args.rand_ind is not None:
        print("Loading random indices: {}...".format(args.rand_ind))
        bcombos = np.loadtxt(args.rand_ind, dtype=int)

    data_ind = [np.nonzero(ds.sa.chunks == c)[0][0] for c in np.unique(ds.sa.chunks)]
    stat = np.mean(ds.samples[data_ind, :], 0)
    print("Saving mean stat image...")
    nib.save(map2nifti(ds, stat), args.stat)

    print("Calculating threshold...")
    thrmap = gvt.thresh_map(
        ds,
        feature_thresh_prob=args.feature_thresh_prob,
        n_blocks=args.n_blocks,
        n_proc=args.n_proc,
        n_bootstrap=args.n_bootstrap,
        mmap_file=args.mmap_file,
    )

    print("Saving threshold map...")
    nib.save(map2nifti(ds, thrmap), args.thresh)
