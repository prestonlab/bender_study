#!/bin/bash
#
# Prepare permuted statistics images for group volume threshold analysis.

if [[ $# -lt 1 ]]; then
    echo "Prepare template space slices for GVT analysis."
    echo
    echo "Usage: bender_prep_gvt.sh filepath"
    echo
    exit 1
fi

mask=$STUDYDIR/gptemplate/highres_brain_all/gp_template_mni_affine_mask.nii.gz
while getopts ":m:" opt; do
    case $opt in
    m)
        mask=$OPTARG
        ;;
    *)
        echo "Unknown option: ${opt}"
        exit 1
        ;;
    esac
done
shift $((OPTIND-1))

filepath=$1

outdir=$SCRATCH/$STUDY/batch/$filepath
mkdir -p "${outdir}"
template=$STUDYDIR/gptemplate/highres_brain_all/gp_template_mni_affine.nii.gz

# split the mask into slices to match the standard space images
cp "${template}" "${outdir}"/bg_image.nii.gz
cp "${mask}" "${outdir}"/mask.nii.gz
mkdir -p "${outdir}"/mask_z
fslsplit "${mask}" "${outdir}"/mask_z/mask_z -z
