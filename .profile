# Profile to be sourced before running bender_study scripts.

case $USER in
    morton)
        environ=local
        ;;
    mortonne)
        environ=remote
        ;;
    *)
        echo "Error: unknown user $USER. Please edit .profile for your environment."
        ;;
esac

if [[ $environ = remote ]]; then
    # activate LMod modules for Lonestar 6
    module load python3/3.9.7
    module load launcher/3.10
    module use /work2/03206/mortonne/software/modules
    module load fsl/5.0.11
    module load ants/2.1.0
    module load gsl
    module load afni/16.1.20

    # activate the Python virtual environment
    source $STOCKYARD/software/venv/bender_study/bin/activate

    # subjects in main analysis (good images and behavior)
    export SUBJIDFORMAT=bender_%02d
    export SUBJNOS=2:4:5:6:7:8:10:11:12:13:14:15:16:17:18:20:21:22:23:24:25:26:28:29:31:32:34:35:37:38
    export SUBJIDS=$(subjids $SUBJNOS)

    # lists of functional runs
    export LOCRUNS=loc_1:loc_2:loc_3:loc_4
    export PREXRUNS=prex_1:prex_2:prex_3:prex_4:prex_5:prex_6
    export STUDYRUNS=study_1:study_2:study_3:study_4:study_5:study_6
    export TESTRUNS=test_1:test_2:test_3:test_4:test_5
    export ALLRUNS=$LOCRUNS:$PREXRUNS:$STUDYRUNS:$TESTRUNS

    # base study directory
    export STUDY=bender
    export STUDYDIR=$STOCKYARD/lonestar/bender

    # directory to save job output
    export BATCHDIR=$STOCKYARD/lonestar/bender/batch/launchscripts
fi

# variables for notebooks
export BENDER_BIDS=$HOME/Dropbox/data/bender
export BENDER_FIGURES=$HOME/Dropbox/work/bender/figs/react_r1
export BENDER_RESULTS=$HOME/work/bender
