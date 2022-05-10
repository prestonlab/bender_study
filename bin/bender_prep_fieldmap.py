#!/usr/bin/env python
#
# Prepare fieldmaps for unwarping.

from bender_study.subjutil import *

parser = SubjParser()
parser.add_argument("--dte", help="delta TE", default=2.46)
args = parser.parse_args()

from bender_study.bender import BenderPath

bp = BenderPath(args.subject, args.study_dir)
log = bp.init_log("prepfm", "preproc", args)

reg_xfm = bp.path("fieldmap", "antsreg", "transforms")

log.start()
log.run("mkdir -p %s" % reg_xfm)

run_fm = bp.run_fieldmap()
fm_anat = bp.fieldmap_anat()
fm_runs = set([item for sublist in run_fm.values() for item in sublist])
fm_dir = bp.path("fieldmap")
for i in fm_runs:
    # use the first magnitude image
    mag = impath(fm_dir, "fieldmap_mag%d" % i)

    # correct bias of magnitude image
    mag_cor = impath(fm_dir, "fieldmap_mag_cor%d" % i)
    log.run("N4BiasFieldCorrection -d 3 -i %s -o %s" % (mag, mag_cor))

    # determine correct highres to use (assuming better registration
    # within days)
    highres_brain = bp.image_path("anatomy", "orig_brain%d" % fm_anat[i])
    highres_mask = bp.image_path("anatomy", "brainmask%d" % fm_anat[i])

    # register the corrected magnitude image to the corrected highres
    xfm_base = os.path.join(reg_xfm, "fieldmap%d-orig%d_" % (i, fm_anat[i]))
    xfm_file = xfm_base + "0GenericAffine.mat"
    cmd = "antsRegistration -d 3 -r [{ref},{mov},1] -t Rigid[0.1] -m MI[{ref},{mov},1,32,Regular,0.25] -c [1000x500x250x100,1e-6,10] -f 8x4x2x1 -s 3x2x1x0vox -n BSpline -w [0.005,0.995] -o {xfm}".format(
        ref=highres_brain, mov=mag_cor, xfm=xfm_base
    )
    log.run(cmd)

    # use the highres brain mask to mask the magnitude image
    mask_reg = impath(fm_dir, "brainmask%d" % i)
    mag_brain = impath(fm_dir, "fieldmap_mag_brain%d" % i)
    mag_cor_brain = impath(fm_dir, "fieldmap_mag_cor_brain%d" % i)
    log.run(
        "antsApplyTransforms -i %s -o %s -r %s -t [%s,1] -n NearestNeighbor"
        % (highres_mask, mask_reg, mag, xfm_file)
    )
    log.run("fslmaths %s -fillh26 %s" % (mask_reg, mask_reg))
    log.run("fslmaths %s -mas %s %s" % (mag, mask_reg, mag_brain))
    log.run("fslmaths %s -mas %s %s" % (mag_cor, mask_reg, mag_cor_brain))

    # convert the phase image to radians
    phase = impath(fm_dir, "fieldmap_phase%d" % i)
    rads = impath(fm_dir, "fieldmap_rads_brain%d" % i)
    cmd = "fsl_prepare_fieldmap SIEMENS %s %s %s %.2f" % (
        phase,
        mag_brain,
        rads,
        args.dte,
    )
    log.run(cmd)

log.finish()
