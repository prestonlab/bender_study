#!/bin/bash
#
# Collect slice-level group volume threshold results.

if [[ $# -lt 2 ]]; then
    echo "Collect slice-level GVT results and backup results."
    echo
    echo "Usage:   bender_post_gvt.sh filepath subjects"
    echo 'Example: bender_post_gvt.sh mvpa/cat_react_sim $SUBJNOS'
    echo
    exit 1
fi

filepath=$1
subjects=$2

cd "${SCRATCH}/${STUDY}/batch/${filepath}" || exit

# check input files
if [[ ! -f mask.nii.gz ]]; then
    echo "Error: mask file missing."
    exit 1
fi

nslices=$(fslval mask.nii.gz dim3)
maxsliceno=$(printf '%04d' $((nslices-1)))
if [[ ! -f thresh_z${maxsliceno}.nii.gz ]]; then
    echo "Error: threshold slices missing."
    exit 1
fi

if [[ ! -f stat_z${maxsliceno}.nii.gz ]]; then
    echo "Error: stat slices missing."
    exit 1
fi

echo "merging slices..."
fslmerge -z thresh thresh_z*
if [[ -s thresh.nii.gz ]]; then
    # if threshold volume successfully assembled, then delete the
    # slice images to save space
    rm thresh_z*

    # remove the individual subject slice images. Keep the group slice
    # datasets for now, in case we want to rerun gvt with a different
    # voxel alpha or something
    for subject in $(subjids -s ' ' "${subjects}"); do
    	rm -f "${SCRATCH}/${STUDY}/${subject}/${filepath}"/allperm_std_z*.nii.gz
    done
else
    echo "Error: problem merging threshold slices."
    exit 1
fi

echo "merging slices..."
fslmerge -z stat stat_z*
if [[ -s stat.nii.gz ]]; then
    # if threshold volume successfully assembled, then delete the
    # slice images to save space
    rm stat_z*
else
    echo "Error: problem merging stat slices."
    exit 1
fi

echo "creating thresholded stat image..."
fslmaths stat -sub thresh -bin thresh_mask
fslmaths stat -mas thresh_mask stat_thresh

echo "backing up main results on WORK..."
wdir=${STUDYDIR}/batch/${filepath}
mkdir -p "${wdir}"
imcp bg_image mask thresh thresh_mask stat stat_thresh "${wdir}"
