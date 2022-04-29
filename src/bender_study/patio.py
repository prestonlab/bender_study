"""Input/output utilities for use with PAT data."""

from collections import OrderedDict
import numpy as np
import os
import re
from events import Events


def read_log(filepath):
    """Read a standard PAT log file for one run."""

    # read the header
    with open(filepath) as f:
        l = f.readline()
    names = l[:-1].split("\t")

    t = np.loadtxt(filepath, skiprows=1)
    events = Events()
    for i in range(t.shape[0]):
        d = OrderedDict()
        for j in range(t.shape[1]):
            d[names[j]] = t[i, j]
        events.append(d)
    return events


def event_array(events, key):
    """Get a field of events as an array."""
    return np.array([e[key] for e in events])


def find_run_log(log_dir, period, run):
    """Find the log for a given run."""

    log_test = re.compile("^log_\d+-\d+-\d+_\D+.txt$")
    state_test = re.compile("\d+-\d+-\d+")

    # start with a list of all files and directories
    d = os.listdir(log_dir)
    found = False
    for name in d:
        full = os.path.join(log_dir, name)
        if log_test.match(name) and os.path.isfile(full):
            # this is a file that looks like a log
            match = state_test.search(name)
            if not match:
                continue
            else:
                # get the state string
                state = match.group()

            # determine run and period from the filename
            log_run = int(state.split("-")[1])
            log_period = name[match.end() + 1 : -4]

            # check if this is the requested run
            if log_run == run and log_period == period:
                found = True
                return full
    if not found:
        raise IOError("log not found.")


def event_match(events, evdef):
    """Find events that match a set of conditions."""
    inc = np.array([True for e in events])
    for cond in evdef:
        if cond[1] == "nan":
            # must use special test for NaNs
            cond_match = np.array([np.isnan(e[cond[0]]) for e in events])
        elif cond[1] == "!nan":
            cond_match = np.array([not np.isnan(e[cond[0]]) for e in events])
        else:
            cond_match = np.array([e[cond[0]] == cond[1] for e in events])
        inc = np.logical_and(inc, cond_match)
    return inc


def write_ev_onsets(
    events,
    ev_def,
    ev_name,
    out_dir,
    runid,
    onset_key="onset",
    weight_key=None,
    duration=None,
):
    """Write onset files compatible with FSL."""

    # get the events matching this EV definition
    inc = events.match(**ev_def)

    # get the information we need from the events
    onsets = events.array(onset_key)[inc]
    if duration is None:
        durations = events.array("duration")[inc]
    if weight_key is not None:
        weights = events.array(weight_key)[inc]

    # write out a standard EV onsets file
    filename = "%s_%s.txt" % (ev_name, runid)
    filepath = os.path.join(out_dir, filename)
    with open(filepath, "w") as f:
        for i, time in enumerate(onsets):
            if weight_key is not None:
                w = weights[i]
            else:
                w = 1
            if duration is None:
                trial_duration = durations[i]
            else:
                trial_duration = duration
            f.write("%.8f\t%.8f\t%.8f\n" % (time, trial_duration, w))


def write_mult_ev_onsets(events, ev_dict, out_dir, runid, **kwargs):
    """Write FSL onset files for multiple explanatory variables."""
    for key, val in ev_dict.items():
        write_ev_onsets(events, val, key, out_dir, runid, **kwargs)


def read_period_events(log_dir, period, n_run):
    """Read all events for a period."""

    events = Events()
    for i in range(n_run):
        log_file = find_run_log(log_dir, period, i + 1)
        run_events = read_log(log_file)
        events.extend(run_events)
    return events


def add_event_info(ds, events, shift=0.0):
    """Add information about events to a dataset."""

    from mvpa2.datasets import eventrelated as evr

    # translate the events into attribute format
    for key in events[0].keys():
        val = events[0][key]
        if isinstance(val, str):
            none_val = "rest"
        else:
            none_val = None
        ds.sa[key] = evr.events2sample_attr(
            events,
            ds.sa.time_coords,
            noinfolabel=none_val,
            condition_attr=key,
            onset_shift=shift,
        )

    onsets = np.array([e["onset"] for e in events]) + shift
    eindex = []
    etime = []
    etimeindex = []
    prev_index = None
    for i in range(ds.shape[0]):
        # find events starting before or on this sample
        match = np.nonzero(ds.sa.time_coords[i] >= onsets)[0]
        if len(match) != 0:
            # latest event starting before or at this time
            event_index = match[-1]
            eindex.append(event_index)

            # onset of this event
            event_onset = onsets[event_index]
            if prev_index is None or event_index != prev_index:
                # reset the event time index
                event_time_index = 0
            else:
                event_time_index += 1

            # time relative to event onset
            etime.append(ds.sa.time_coords[i] - event_onset)
            etimeindex.append(event_time_index)
            prev_index = event_index
        else:
            # this sample corresponds to no event
            eindex.append(None)
            etime.append(None)
            etimeindex.append(None)

    # add the new attributes to the dataset
    ds.sa["event_indices"] = eindex
    ds.sa["event_time_indices"] = etimeindex
    ds.sa["event_time_coords"] = etime


def combine_datasets(ds_tuple):
    """Combine datasets while dealing with attribute differences."""

    from mvpa2.base.dataset import vstack

    # find the intersection of all sample attributes
    s1 = set(ds_tuple[0].sa.keys())
    key_intersect = s1
    for i in range(1, len(ds_tuple)):
        s = set(ds_tuple[i].sa.keys())
        key_intersect = s1.intersection(s)

    # remove sample attributes not in the intersection; add a dataset
    # attribute
    for i, d in enumerate(ds_tuple):
        for k in d.sa.keys():
            if k not in key_intersect:
                del d.sa[k]
        d.sa["dataset"] = [i]

    # combine the datasets
    ds_full = vstack(ds_tuple)
    return ds_full
