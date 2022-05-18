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
export STUDYDIR=$STOCKYARD/lonestar/bender

# directory to save job output
export BATCHDIR=$STOCKYARD/lonestar/bender/batch/launchscripts
