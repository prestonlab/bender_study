#!/usr/bin/env python

from bender_study.subjutil import *
import os

parser = SubjParser()
parser.add_argument('runid', help="run identifier")
parser.add_argument('-k', '--keep', help="keep intermediate files",
                    action='store_true')
args = parser.parse_args()

from bender_study.bender import BenderPath

bp = BenderPath(args.subject, args.study_dir)
log = bp.init_log('epireg_%s' % args.runid, 'preproc', args)

log.start()

number = int(args.subject.split('_')[1])
if number >= 22:
    ees = 0.000385 # echo spacing 0.77 ms / GRAPPA 2
else:
    ees = 0.00047 # echo spacing 0.94 ms / GRAPPA 2
pedir = 'y-'

run_fm = bp.run_fieldmap()
map_dir = bp.path('fieldmap')
phase = args.runid.split('_')[0]
run = int(args.runid.split('_')[1])
fm_scan = run_fm[phase][run-1]
fm_anat = bp.fieldmap_anat()
anat_scan = fm_anat[fm_scan]

# structural files (use the one collected that day)
highres = bp.image_path('anatomy', 'orig%d' % anat_scan)
highres_brain = bp.image_path('anatomy', 'orig_brain%d' % anat_scan)
highres_mask = bp.image_path('anatomy', 'brainmask%d' % anat_scan)
wm_mask = bp.image_path('anatomy', 'wm%d' % anat_scan)

# epi files
epi_dir = bp.path('bold', args.runid)
epi_input = impath(epi_dir, 'bold_cor_mcf_avg')
epi_output = impath(epi_dir, 'bold_cor_mcf_avg_unwarp_brain')

# prepare output directory
fm_dir = bp.path('bold', args.runid, 'fm')
out_base = os.path.join(fm_dir, 'epireg')
log.run('mkdir -p %s' % fm_dir)

fmap = impath(map_dir, 'fieldmap_rads_brain%d' % fm_scan)
fmapmag = impath(map_dir, 'fieldmap_mag_cor%d' % fm_scan)
fmapmagbrain = impath(map_dir, 'fieldmap_mag_cor_brain%d' % fm_scan)

# run epi_reg
cmd = 'epi_reg_ants --fmap=%s --fmapmag=%s --fmapmagbrain=%s --wmseg=%s --echospacing=%.06f --pedir=%s -v --epi=%s --t1=%s --t1brain=%s --out=%s --noclean' % (
    fmap, fmapmag, fmapmagbrain, wm_mask, ees, pedir, epi_input, highres,
    highres_brain, out_base)
log.run(cmd)

# convert shift map to a warp
shift = impath(fm_dir, 'epireg_fieldmaprads2epi_shift')
warp = impath(fm_dir, 'epireg_epi_warp')
log.run('convertwarp -r %s -s %s -o %s --shiftdir=%s --relout' % (
    epi_input, shift, warp, pedir))

# unwarp the average run image for registration purposes
epi_unwarped = impath(fm_dir, 'epireg_epi_unwarped')
log.run('applywarp -i %s -r %s -o %s -w %s --interp=spline --rel' % (
    epi_input, epi_input, epi_unwarped, warp))

# transform the anatomical brain mask into functional space
mask_reg = impath(fm_dir, 'brainmask')
log.run('flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour' % (
    highres_mask, epi_unwarped, os.path.join(fm_dir, 'epireg_inv.mat'),
    mask_reg))

# dilate to make a tighter brain extraction than the liberal one
# originally used for the functionals
log.run('fslmaths %s -kernel sphere 3 -dilD %s' % (mask_reg, mask_reg))

# mask the unwarped epi
log.run('fslmaths %s -mas %s %s' % (epi_unwarped, mask_reg, epi_output))

if not args.keep:
    log.run('rm -f %s/{epireg,epireg_fieldmap*,epireg_1vol,epireg_fast*,epireg_epi_unwarped,epireg_warp}.nii.gz' % fm_dir)

log.finish()
