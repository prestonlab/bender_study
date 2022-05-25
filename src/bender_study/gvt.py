"""Functions for running group volume thresholding."""

import os
import random
import numpy as np
from mvpa2.base.dataset import vstack
from mvpa2.datasets.mri import fmri_dataset


def loadcat(srcfiles, maskfile):
    """Load files as datasets and concatenate them."""
    for i, file in enumerate(srcfiles):
        print("Loading %s" % file)
        if i == 0:
            ds = fmri_dataset(file, mask=maskfile)
            ds.sa['chunks'] = [i]
            a = ds.a
        else:
            newds = fmri_dataset(file, mask=maskfile)
            newds.sa['chunks'] = [i]
            ds = vstack((ds, newds))
    ds.a = a
    return ds


def get_thresholding_map(data, p=0.001):
    """Return array of thresholds corresponding to a probability of such value in the input

    Thresholds are returned as an array with one value per column in the input
    data.

    Parameters
    ----------
    data : 2D-array
      Array with data on which the cumulative distribution is based.
      Values in each column are sorted and the value corresponding to the
      desired probability is returned.
    p : float [0,1]
      Value greater or equal than the returned threshold have a probability `p` or less.
    """
    # we need NumPy indexing logic, even if a dataset comes in
    data = np.asanyarray(data)
    p_index = int(len(data) * p)
    if p_index < 1:
        raise ValueError(
            "requested probability is too low for the given number of samples"
        )
    # threshold indices are all in one row of the argsorted inputs
    thridx = np.argsort(data, axis=0, kind="quicksort")[-p_index]
    return data[thridx, np.arange(data.shape[1])]


def thresh_map(
    ds,
    chunk_attr="chunks",
    n_bootstrap=100000,
    n_blocks=1,
    n_proc=1,
    feature_thresh_prob=0.001,
    mmap_file=None,
    bcombos=None,
):

    print("Preparing bootstrap...")
    chunk_samples = dict(
        [
            (c, np.where(ds.sa[chunk_attr].value == c)[0])
            for c in ds.sa[chunk_attr].unique
        ]
    )
    # pre-built the bootstrap combinations
    if bcombos is None:
        bcombos = [
            [random.sample(v, 1)[0] for v in chunk_samples.values()]
            for i in xrange(n_bootstrap)
        ]
        bcombos = np.array(bcombos, dtype=int)
    segwidth = ds.nfeatures / n_blocks
    # speed things up by operating on an array not a dataset
    ds_samples = ds.samples

    if mmap_file is not None:
        from joblib import load, dump

        if not os.path.exists(mmap_file):
            print("Writing memory mapped file: {}".format(mmap_file))
            _ = dump(ds_samples, mmap_file)
        else:
            print("Loading memory mapped file: {}".format(mmap_file))
        ds_samples = load(mmap_file, mmap_mode="r+")

    def featuresegment_producer(ncols):
        for segstart in xrange(0, ds.nfeatures, ncols):
            # one average map for every stored bcombo
            # this also slices the input data into feature subsets
            # for the compute blocks
            yield [
                np.mean(
                    # get a view to a subset of the features
                    # -- should be somewhat efficient as feature axis is
                    # sliced
                    ds_samples[sidx, segstart : segstart + ncols],
                    axis=0,
                )
                for sidx in bcombos
            ]

    if n_proc == 1:
        # Serial execution
        thrmap = np.hstack(  # merge across compute blocks
            [
                get_thresholding_map(d, feature_thresh_prob)
                # compute a partial threshold map for as many features
                # as fit into a compute block
                for d in featuresegment_producer(segwidth)
            ]
        )
    else:
        # Parallel execution
        verbose_level_parallel = 0
        # local import as only parallel execution needs this
        from joblib import Parallel, delayed
        from joblib.pool import has_shareable_memory

        if mmap_file is not None:
            print("Starting bootstrap with memory mapping...")
            thrmap = np.hstack(
                Parallel(
                    n_jobs=n_proc,
                    max_nbytes=None,
                    pre_dispatch=n_proc,
                    verbose=verbose_level_parallel,
                )(
                    delayed(get_thresholding_map)(d, feature_thresh_prob)
                    for d in featuresegment_producer(segwidth)
                )
            )

        else:
            # same code as above, just in parallel with joblib's Parallel
            print("Starting bootstrap...")
            thrmap = np.hstack(
                Parallel(
                    n_jobs=n_proc, pre_dispatch=n_proc, verbose=verbose_level_parallel
                )(
                    delayed(get_thresholding_map)(d, feature_thresh_prob)
                    for d in featuresegment_producer(segwidth)
                )
            )
    return thrmap
