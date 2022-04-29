"""Tools for analysis of Bender data."""

import pickle
import warnings

from bender_study.subjutil import *
from bender_study import patio
from bender_study.events import Events


def res_list(
    filepat,
    exclude=[
        None,
    ],
):
    res = []
    subjects = os.environ["SUBJIDS"].strip().split(":")
    inc_subj = []
    for subj in subjects:
        if subj in exclude:
            continue

        res_file = filepat % subj
        if not os.path.exists(res_file):
            continue

        with open(res_file, "r") as f:
            res_subj = pickle.load(f)
        res.append(res_subj)
        inc_subj.append(subj)
    return res, inc_subj


class BenderPath(SubjPath):
    def __init__(self, subject, study_dir=None):
        SubjPath.__init__(self, subject, study_dir)

    def get_category(self, cond):
        cond_dict = {1: "face", 2: "face", 3: "scene", 4: "scene", 5: "object"}
        subcond_dict = {1: "female", 2: "male", 3: "manmade", 4: "natural", 5: "object"}
        return cond_dict[cond], subcond_dict[cond]

    def read_log(self, period, run, duration=None, verbose=False):
        """Read a Bender log file for one run."""

        # find the log for this run
        log_dir = self.path("behav", "log")
        log = patio.find_run_log(log_dir, period, run)

        # read in raw events
        if verbose:
            print("Reading events from %s" % log)
        events = patio.read_log(log)

        # add duration field if specified
        if duration is not None:
            for i in range(len(events)):
                events[i]["duration"] = duration

        # translate condition codes to string equivalents
        for i in range(len(events)):
            cat, subcat = self.get_category(events[i]["cond"])
            events[i]["condition"] = cat
            events[i]["category"] = cat
            events[i]["subcategory"] = subcat
        return events

    def read_loc_blocks(self, run):
        """Read information about blocks in a localizer run."""

        import numpy as np

        n_run = self.get_n_runs("c_loc")

        # load information about picture stimuli
        stim_events = self.read_log("c_loc", run, verbose=True)

        # add block start and duration info
        block_events = Events()
        block_events.append(
            {
                "onset": 0,
                "block": 1,
                "category": "rest",
                "subcategory": "rest",
                "duration": 18,
            }
        )
        blocks = stim_events.array("block")
        ublocks = np.unique(blocks)
        for i, b in enumerate(ublocks):
            first = np.where(blocks == b)[0][0]
            e = stim_events[first]
            cat, subcat = self.get_category(e["cond"])
            block_events.append(
                {
                    "onset": e["onset"],
                    "block": i + 2,
                    "category": cat,
                    "subcategory": subcat,
                    "duration": 20,
                }
            )

        onset = block_events[-1]["onset"] + 20
        block_events.append(
            {
                "onset": onset,
                "block": len(ublocks) + 2,
                "category": "rest",
                "subcategory": "rest",
                "duration": 18,
            }
        )
        return block_events

    def read_period(self, period, duration=None, verbose=False):
        """Read all events for a period of Bender."""

        n_run = self.get_n_runs(period)
        events = Events()
        for i in range(n_run):
            run_events = self.read_log(period, i + 1, duration, verbose)
            for j in range(len(run_events)):
                run_events[j]["chunks"] = i + 1
            events.extend(run_events)
        return events

    def get_test_perf(self, field="correct"):
        """Read performance for all tests."""

        res = {}

        # AB training
        res["ab1"] = self.read_log("ab_feedback", 1).sort("group").array(field)
        res["ab2"] = self.read_log("ab_feedback", 2).sort("group").array(field)
        res["ab3"] = self.read_log("ab_feedback", 3).sort("group").array(field)
        res["ab4"] = self.read_log("ab_feedback", 4).sort("group").array(field)

        # BC/XY
        events = self.read_period("bcxy_test").sort("group")
        res["xy"] = events.array(field)[events.array("cond") == 5]
        res["bc"] = events.array(field)[events.array("cond") < 5]

        # AC
        res["ac"] = self.read_period("ac_test").sort("group").array(field)

        # AB final
        res["ab"] = self.read_period("ab_test").sort("group").array(field)

        return res

    def write_period_onsets(self, period, duration, ev_def, ev_name):
        """Write onset files for one condition for all runs in a period."""
        n_run = self.get_n_runs(period)
        run_period = period.split("_")[1]
        out_dir = self.path("behav", run_period)
        for i in range(1, n_run + 1):
            events = self.read_log(period, i, duration)
            patio.write_ev_onsets(events, ev_def, ev_name, out_dir, i)

    def write_mult_period_onsets(self, period, duration, ev_dict):
        """Write multiple onset files for each run in a period."""
        for key, val in ev_dict.items():
            self.write_period_onsets(period, duration, val, key)

    def get_n_runs(self, period):
        """Number of runs in a given period."""

        # period code may be the full behavioral code or the shortened
        # fMRI code
        d = {
            "p_loc": 1,
            "p_prex": 1,
            "c_loc": 4,
            "a_prex": 6,
            "p_study": 1,
            "p_feedback": 1,
            "ab_study": 4,
            "ab_feedback": 4,
            "bcxy_study": 6,
            "ac_test": 5,
            "bcxy_test": 1,
            "ab_test": 1,
            "loc": 4,
            "prex": 6,
            "study": 6,
            "test": 5,
        }
        return d[period]

    def get_all_runs(self, filename=None, check=False, subdir=None):
        """Generate all standard BOLD run files."""

        phases = ["loc", "prex", "study", "test"]
        n_runs = [4, 6, 6, 5]
        paths = dict()
        for i in range(len(phases)):
            phase_paths = []
            for j in range(n_runs[i]):
                run_name = "%s_%d" % (phases[i], j + 1)
                if subdir is not None:
                    run_dir = self.path("bold", subdir)
                    filename = run_name + ".nii.gz"
                else:
                    run_dir = self.path("bold", run_name)

                if not os.path.exists(run_dir):
                    warnings.warn("directory does not exist: %s" % run_dir)
                    continue

                if filename is not None:
                    # return the path to a specific file
                    run_file = impath(run_dir, filename)
                    if check and not os.path.exists(run_file):
                        raise IOError("File does not exist: %s" % run_file)
                    phase_paths.append(run_file)
                else:
                    # return just the directory
                    phase_paths.append(run_dir)
                paths[phases[i]] = phase_paths
        return paths

    def get_flat_runs(self, filename=None, check=False, subdir=None):
        """Get a flat list of standard BOLD run files."""

        run_dict = self.get_all_runs(filename, check, subdir)
        keys = run_dict.keys()
        runs = []
        for phase in keys:
            phase_runs = run_dict[phase]
            for run in phase_runs:
                runs.append(run)
        return runs

    def run_files(
        self, period, run=None, suffix="", subdir="antsreg/data", sepdirs=False
    ):
        """Get a list of run files for a Bender period."""

        # search for a specific run or all runs in this period
        if run:
            file_pattern = "%s_%d" % (period, run)
        else:
            file_pattern = "%s_\d" % period
        files = self.bold_files(
            sepdirs=sepdirs, file_pattern=file_pattern, subdir=subdir, suffix=suffix
        )
        return files

    def block_dataset(self, period, duration, mask=None, shift=0.0, suffix=""):
        """Read in BOLD data and label volumes with conditions."""

        from mvpa2.base.dataset import vstack
        from mvpa2.datasets.mri import fmri_dataset

        # get the files for BOLD dataa
        run_period = period.split("_")[1]
        files = self.run_files(run_period, suffix=suffix)
        ds = []
        for i in range(len(files)):
            # read the log and data
            events = self.read_log(period, i + 1, duration, verbose=True)
            print("Reading %s-%d data from %s" % (period, i + 1, files[i]))
            print("Masking with %s" % mask)
            run_ds = fmri_dataset(files[i], mask=mask)

            # add events information to the dataset
            patio.add_event_info(run_ds, events, shift)

            run_ds.sa["chunks"] = [i + 1]
            ds.append(run_ds)
        # concatenate all runs for this period
        fds = vstack(ds, a=0)
        return fds

    def event_dataset(self, period, duration, mask=None):
        """Read BOLD data with time points x voxels as features."""

        from mvpa2.base.dataset import vstack
        from mvpa2.datasets.mri import fmri_dataset
        from mvpa2.datasets import eventrelated as ev

        run_period = period.split("_")[1]
        files = self.run_files(run_period)
        ds = []
        for i, file in enumerate(files):
            events = self.read_log(period, i + 1, duration)
            runds = fmri_dataset(files[i], mask=mask)
            evds = ev.extract_boxcar_event_samples(
                runds, events=events, time_attr="time_coords"
            )
            ds.append(evds)
        # concatenate all runs
        fds = vstack(ds, a=0)
        return fds

    def beta_dataset(self, period, mask=None, suffix="_stim", include=None):
        """Read a betaseries with multiple runs."""

        from mvpa2.datasets.mri import fmri_dataset

        # load each run of behavioral data
        n_run = self.get_n_runs(period)
        events = Events()
        for i in range(n_run):
            if period == "c_loc":
                ev_time = self.read_loc_blocks(i + 1)
                ev_group = ev_time.reduce("block", rm_varying=True)
            else:
                ev_time = self.read_log(period, i + 1, verbose=True)
                ev_group = ev_time.reduce("group", rm_varying=True)

            # set the run
            ev_group.setfield("chunks", i + 1)

            # add to the full list of events
            events.extend(ev_group)

        # load functional data
        run_period = period.split("_")[1]
        beta_dir = os.path.join(
            self.study_dir, "batch", "glm", run_period + suffix, "beta"
        )
        beta_file = impath(beta_dir, self.subject + "_beta")
        print("Reading %s betaseries data from %s" % (period, beta_file))
        print("Masking with %s" % mask)
        if include is not None:
            ds = fmri_dataset(beta_file, mask=mask, add_fa={"include": include})
        else:
            ds = fmri_dataset(beta_file, mask=mask)

        # add attributes from events
        for key in events.keys():
            try:
                ds.sa[key] = events.array(key)
            except ValueError:
                print("Could not convert %s field." % key)
        return ds

    def get_fieldmaps(self, postfix="mag", subdir=None):
        """Get a list of all fieldmaps for a subject."""

        if subdir is not None:
            basepath = self.path("fieldmap", subdir)
        else:
            basepath = self.path("fieldmap")

        max_fm = 6
        all_fm = []
        for i in range(max_fm):
            mag_file = impath(basepath, "fieldmap_mag%d" % (i + 1))
            filepath = impath(basepath, "fieldmap_%s%d" % (postfix, (i + 1)))
            if os.path.exists(mag_file):
                all_fm.append(filepath)
        return all_fm

    def get_fieldmap_dict(self, subdir=None):
        """Get a dictionary with all types of fieldmap images."""

        fm_types = ["mag", "mag_brain", "phase", "rads_brain"]
        fm = dict()
        for i in range(len(fm_types)):
            fm[fm_types[i]] = self.get_fieldmaps(fm_types[i], subdir)

        fm2 = []
        for i in range(len(fm["mag"])):
            d = dict()
            for key in fm.keys():
                d[key] = fm[key][i]
            fm2.append(d)
        return fm2

    def phase_fieldmap(self):
        """Get a dictionary mapping experiment phase to fieldmap scan number."""

        # defaults are specified here; exceptions are coded in
        # run_fieldmap
        n = len(self.get_fieldmaps())
        phase_fm = dict()
        if n == 2:
            # for some early participants, just ran two fieldmaps
            phase_fm["loc"] = 1
            phase_fm["prex"] = 1
            phase_fm["study"] = 2
            phase_fm["test"] = 2
        elif n == 3:
            # for most, ran a fieldmap after study
            phase_fm["loc"] = 1
            phase_fm["prex"] = 1
            phase_fm["study"] = 2
            phase_fm["test"] = 3
        else:
            # some participants needed another fieldmap during prex
            phase_fm["loc"] = 1
            phase_fm["prex"] = 2
            phase_fm["study"] = 3
            phase_fm["test"] = 4
        return phase_fm

    def run_fieldmap(self):
        """Get fieldmap scan number to use for each functional scan."""

        run_fm = dict()
        phase_fm = self.phase_fieldmap()
        for (phase, fm) in phase_fm.iteritems():
            n_runs = self.get_n_runs(phase)
            run_fm[phase] = [fm] * n_runs

        if self.subject == "bender_08":
            # localizer failed on day 1; reran on day 2. Must use
            # fieldmap from day 2 for unwarping
            run_fm["loc"] = [2, 2, 2, 2]
        elif self.subject == "bender_15":
            # day 2 fieldmap 2 accidentally run with bad prescription,
            # but immediately redone. So same as standard
            # three-fieldmap run, but using fieldmaps 1, 2, and 4
            run_fm["prex"] = [1, 1, 1, 1, 1, 1]
            run_fm["study"] = [2, 2, 2, 2, 2, 2]
            run_fm["test"] = [4, 4, 4, 4, 4]
        elif self.subject == "bender_16":
            # had to take out of scanner after prex run 2; ran new
            # fieldmap for the last 4
            run_fm["prex"] = [1, 1, 2, 2, 2, 2]
        elif self.subject == "bender_30":
            # took out of scanner after prex 1
            run_fm["prex"] = [1, 2, 2, 2, 2, 2]
            # problems with motion on day 2
            run_fm["study"] = [3, 4, 3, 4, 4, 4]
            run_fm["test"] = [5, 6, 6, 6, 6]
        elif self.subject == "bender_33":
            # ghosting in prex 1, so reran shim and fieldmap
            run_fm["prex"] = [1, 2, 2, 2, 2, 2]
        return run_fm

    def fieldmap_anat(self):
        """Get the anatomical scan for each fieldmap."""

        # look up the fieldmap used for each functional scan
        run_fm = self.run_fieldmap()

        # corresponding anatomical for each fieldmap is 1 for loc and
        # prex, and 2 for study and test
        day1 = [run_fm["loc"], run_fm["prex"]]
        day1_fm = set([run for phase in day1 for run in phase])
        day2 = [run_fm["study"], run_fm["test"]]
        day2_fm = set([run for phase in day2 for run in phase])
        if self.subject == "bender_08":
            # used localizer is actually on day 2
            day1_fm = set(
                [
                    1,
                ]
            )
            day2_fm = set(
                [
                    2,
                ]
            )

        fm_anat = dict()
        n = len(self.get_fieldmaps())
        for i in range(1, n + 1):
            if i in day1_fm:
                fm_anat[i] = 1
            elif i in day2_fm:
                fm_anat[i] = 2
            else:
                print("Warning: could not determine anatomical for fieldmap %d" % i)

        return fm_anat
