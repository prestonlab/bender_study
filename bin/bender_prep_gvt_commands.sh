#!/bin/bash
#
# Generate commands to prepare slices for group voxel threshold analysis.

if [[ $# -lt 2 ]]; then
    echo "Generate commands to prepare slices for GVT analysis."
    echo
    echo "Usage: bender_prep_gvt_commands.sh filepath subjects"
    echo
    exit 1
fi

filepath=$1
nolist=$2

subjids=$(subjids "${nolist}" | tr ':' ' ')

# concatenate images for each slice and convert to dataset
outdir=${SCRATCH}/${STUDY}/batch/${filepath}
template=${STUDYDIR}/gptemplate/highres_brain_all/gp_template_mni_affine.nii.gz
nslices=$(fslval "${template}" dim3)

for i in $(seq 1 "${nslices}"); do
    permlist=""
    sliceno=$(printf '%04d' $((i-1)))
    for id in ${subjids}; do
        perm=${SCRATCH}/${STUDY}/${id}/${filepath}/allperm_std_z${sliceno}.nii.gz
        if [[ -z ${permlist} ]]; then
            permlist=${perm}
        else
            permlist="${permlist} ${perm}"
        fi
    done
    echo "bender_dsmerge.py ${outdir}/mask_z/mask_z${sliceno}.nii.gz ${outdir}/perm_z${sliceno}.hdf5 ${permlist}"
done
