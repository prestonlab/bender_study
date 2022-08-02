"""Analysis of task data."""

import os
import glob
import re
import json
from pkg_resources import resource_filename
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_runs(task):
    run_dict = {
        'localizer': ('func', 4),
        'exposure': ('func', 6),
        'ABstudy': ('beh', 4),
        'ABfeedback': ('beh', 4),
        'BCXYstudy': ('func', 6),
        'ACtest': ('func', 5),
        'BCXYtest': ('beh', 1),
        'ABtest': ('beh', 1),
    }
    if task not in run_dict:
        raise ValueError(f'Invalid task: {task}')
    data_type, n_run = run_dict[task]
    return data_type, n_run


def get_subjects(subset):
    """Get a list of included participant identifiers."""
    json_file = resource_filename('bender_study', 'data/subjects.json')
    with open(json_file, 'r') as f:
        numbers = json.load(f)
    subjects = [f'{n:02d}' for n in numbers[subset]]
    return subjects


def load_run_events(bids_dir, task, subject, run=None):
    """Load events for one run from a BIDS directory."""
    # get the run to load
    data_type, n_run = get_runs(task)
    if n_run == 1:
        if run is not None:
            if run != 1:
                raise ValueError(f'Invalid run number: {run}.')
            else:
                run = None
    elif run is None:
        raise ValueError('Must specify run number.')
    elif run > n_run:
        raise ValueError(f'Invalid run number: {run}.')

    # get the events file directory
    subj_dir = os.path.join(bids_dir, f'sub-{subject}', data_type)
    if not os.path.exists(subj_dir):
        raise IOError(f'Subject directory does not exist: {subj_dir}')

    # load the events
    if run is not None:
        file = os.path.join(
            subj_dir, f'sub-{subject}_task-{task}_run-{run}_events.tsv'
        )
    else:
        file = os.path.join(subj_dir, f'sub-{subject}_task-{task}_events.tsv')
    data = pd.read_table(file)

    # add a run column
    if run is not None:
        data['run'] = run
    else:
        data['run'] = 1
    return data


def load_subject_events(bids_dir, task, subject, runs=None):
    """Load events for one subject from a BIDS directory."""
    # directory with the events file
    data_type, n_run = get_runs(task)

    # number of runs for this task
    if runs is None:
        if n_run != 1:
            runs = range(1, n_run + 1)
        else:
            runs = None

    # load all task runs
    if runs is not None:
        data = pd.concat(
            [load_run_events(bids_dir, task, subject, run) for run in runs],
            axis=0,
            ignore_index=True,
        )
    else:
        data = load_run_events(bids_dir, task, subject)
    data['subject'] = subject
    return data


def load_merged_study_events(bids_dir, subject, runs=None):
    """Load BCXYstudy events with related test information."""
    study = load_subject_events(bids_dir, 'BCXYstudy', subject, runs)
    test = load_test_events(bids_dir, subjects=[subject])
    perf = test.groupby(['group', 'trial_type'])['correct'].mean()
    merged = pd.concat([study.set_index('group'), perf.unstack('trial_type')], axis=1)
    return merged


def load_test_events(bids_dir, tasks=None, subjects=None, **kwargs):
    """Read all test events."""
    if tasks is None:
        tasks = ['ABfeedback', 'BCXYtest', 'ACtest', 'ABtest']
    if subjects is None:
        subjects = get_subjects('react')
    df_all = []
    for subject in subjects:
        ds = pd.concat(
            [load_subject_events(bids_dir, task, subject, **kwargs) for task in tasks],
            axis=0, ignore_index=True
        )
        df_all.append(ds)
    df = pd.concat(df_all, axis=0, ignore_index=True)
    df['trial_type'] = df['trial_type'].astype('category')
    order = ['AB1', 'AB2', 'AB3', 'AB4', 'AC', 'BC', 'XY', 'AB']
    inc_order = [o for o in order if o in df['trial_type'].unique()]
    df['trial_type'] = df['trial_type'].cat.reorder_categories(inc_order, ordered=True)

    df['pair_type'] = df['pair_type'].astype('category')
    order = ['XY', 'AB', 'BC', 'AC']
    inc_order = [o for o in order if o in df['pair_type'].unique()]
    df['pair_type'] = df['pair_type'].cat.reorder_categories(inc_order, ordered=True)

    df['correct'] = df['response'] == df['target']
    df['correct'] = df['correct'].mask(df['response'].isna())
    return df


def read_images(pool, im_dir):
    """Read image files for a set of items."""
    images = {}
    for i, item in pool.iterrows():
        # file names should have underscores instead of spaces
        file_base = item.stim.replace(' ', '_')
        im_base = os.path.join(im_dir, item.subcategory, file_base)

        # search for a matching image file for this item
        res = [f for f in glob.glob(im_base + '.*') if re.search('\w+\.(png|jpg)', f)]
        if not res:
            raise IOError(f'No file found matching: {im_base}')
        elif len(res) > 1:
            raise IOError(f'Multiple matches for: {im_base}')
        im_file = res[0]

        # read the image
        image = plt.imread(im_file)

        # if grayscale, convert to RGB
        if image.ndim == 2:
            image = np.tile(image[:, :, np.newaxis], [1, 1, 3])
        images[item.stim] = image
    return images
