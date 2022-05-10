#!/usr/bin/env python

import argparse
import os
from glob import glob

from bender_study.subjutil import *

parser = SubjParser()
args = parser.parse_args()

# get information about all directories associated with this subject
sp1 = SubjPath(args.subject, args.study_dir)
sp2 = SubjPath(args.subject + 'a', args.study_dir)
if os.path.exists(os.path.join(args.study_dir, args.subject + 'b')):
    d2split = True
    sp3 = SubjPath(args.subject + 'b', args.study_dir)
    sp = [sp1, sp2, sp3]
    day = [1, 2, 2]
else:
    d2split = False
    sp = [sp1, sp2]
    day = [1, 2]

# do logging in the first subject directory, which will contain
# everything
log = sp1.init_log('clean', 'preproc', args)
log.start()

# remove DTI directory
for i in range(len(sp)):
    dti_dir = sp[i].path('dti')
    if os.path.exists(dti_dir):
        log.run('rm -r ' + dti_dir)
        
# backup anatomical and fieldmap directories
other_anat = []
other_fm = []
for i in range(len(sp)):
    dirpath = sp[0].path('anatomy', 'other_day%d' % day[i])
    log.run('mkdir -p %s' % dirpath)
    other_anat.append(dirpath)

    dirpath = sp[0].path('fieldmap', 'other_day%d' % day[i])
    log.run('mkdir -p %s' % dirpath)
    other_fm.append(dirpath)

# move highres
for i in range(len(sp)):
    highres_file = sp[i].path('anatomy', 'highres001.nii.gz')
    new_file = sp[0].path('anatomy', 'highres%d.nii.gz' % day[i])
    if os.path.exists(highres_file):
        log.run('mv %s %s' % (highres_file, new_file))

# reorient coronals
coronal_ind = 1
for i in range(len(sp)):
    coronals = sp[i].glob('anatomy', 'other', 'T2coronal*')
    coronals.sort()
    for j in range(len(coronals)):
        new_file = sp[0].path('anatomy', 
                              'coronal%d.nii.gz' % coronal_ind)
        log.run('fslreorient2std %s %s' % (coronals[j], new_file))
        coronal_ind += 1

# prepare fieldmaps
fm_ind = 1
for i in range(len(sp)):
    # change to standard filenames. Different versions of dcm2nii
    # create nifti files differently. Assuming the 2011-08-12
    # version, which produces:
    # 1) magnitude 1 - fieldmapsSSSa1001
    # 2) magnitude 2 - fieldmapsSSSa2001
    # 3) phase       - fieldmapsSSSa2001
    # The series number (SSS) varies, but phase is always one more
    # than the magnitude images. Magnitude 1 seems to have less
    # dropout, so using that for registration (unlikely to matter
    # much though).
    fm_files = sp[i].glob('fieldmap', '*001.nii.gz')
    fm_files.sort()
    for j in range(0, len(fm_files), 3):
        mag_file = sp[0].path('fieldmap',
                              'fieldmap_mag%d.nii.gz' % fm_ind)
        phase_file = sp[0].path('fieldmap',
                                'fieldmap_phase%d.nii.gz' % fm_ind)
        log.run('mv %s %s' % (fm_files[j], mag_file))
        log.run('mv %s %s' % (fm_files[j+2], phase_file))
        fm_ind += 1

    # move remaining scans
    scans = sp[i].glob('fieldmap', '*001.nii.gz')
    for f in scans:
        log.run('mv %s %s' % (f, other_fm[i]))
    
# rename functional scans

# remove incomplete runs
sp[0].rm_partial_bold('loc', log)
sp[0].rm_partial_bold('prex', log)
sp[1].rm_partial_bold('study', log)
if d2split:
    sp[2].rm_partial_bold('test', log)
else:
    sp[1].rm_partial_bold('test', log)

bold1 = sp[0].path('bold')
bold2 = sp[1].path('bold')

log.run('rename_func.sh %s %s %s %s' % (
    bold1, bold1, 'functional_loc_', 'loc_'))
log.run('rename_func.sh %s %s %s %s' % (
    bold1, bold1, 'functional_prex_', 'prex_'))
log.run('rename_func.sh %s %s %s %s' % (
    bold2, bold1, 'functional_study_', 'study_'))
if d2split:
    bold3 = sp[2].path('bold')
    log.run('rename_func.sh %s %s %s %s' % (
        bold3, bold1, 'functional_test_', 'test_'))
else:
    log.run('rename_func.sh %s %s %s %s' % (
        bold2, bold1, 'functional_test_', 'test_'))

# merge remaining files across days

# move other anatomical scans to day directories
for i in range(len(sp)):
    log.run('mv %s %s' %
            (sp[i].path('anatomy', 'other', '*'), other_anat[i]))
    log.run('rmdir %s' % sp[i].path('anatomy', 'other'))

# move raw files to new separate subject code directories
log.run('mv %s %s' % (sp[1].path('raw', args.subject + 'a'),
                      sp[0].path('raw')))
if d2split:
    log.run('mv %s %s' % (sp[2].path('raw', args.subject + 'b'),
                          sp[0].path('raw')))

# move existing log files to day directories
for i in range(1, len(sp)):
    log_dir = sp[0].path('logs', 'day%d' % day[i])
    log.run('mkdir -p %s' % log_dir)
    log.run('mv %s %s' % (sp[i].path('logs', '*'), log_dir))

# attempt to remove remaining dirs
day2_dirs = ['anatomy', 'behav', 'bold', 'fieldmap',
             'logs', 'model', 'raw']
for i in range(1, len(sp)):
    for d in day2_dirs:
        if os.path.exists(sp[i].path(d)):
            log.run('rmdir %s' % sp[i].path(d))
    log.run('rmdir %s' % sp[i].subj_dir)

log.finish()
