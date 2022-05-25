#!/bin/bash
#
# Print commands to run all slices in a group volume threshold analysis.

if [[ $# -lt 1 ]]; then
    echo "Print commands to run all slices in a GVT analysis."
    echo
    echo "Usage:   bender_gvt_commands.sh [-o] permdir [gvt options]"
    echo 'Example: bender_gvt_commands $SCRATCH/bender/batch/mvpa/cat_react -a 0.01'
    echo
    echo "-o"
    echo "    Overwrite existing threshold files."
    echo
    exit 1
fi

overwrite=false
while getopts ":o" opt; do
    case $opt in
    o)
        overwrite=true
        ;;
    *)
        echo "Unknown option: ${opt}"
        exit 1
        ;;
    esac
done
shift $((OPTIND-1))

permdir=$1
shift

template=${STUDYDIR}/gptemplate/highres_brain_all/gp_template_mni_affine.nii.gz
nslices=$(fslval "${template}" dim3)

for i in $(seq 1 "${nslices}"); do
    sliceno=$(printf '%04d' $((i-1)))
    permds=${permdir}/perm_z${sliceno}.hdf5
    stat=${permdir}/stat_z${sliceno}.nii.gz
    thresh=${permdir}/thresh_z${sliceno}.nii.gz
    
    if [[ ${overwrite} = true || ! -s ${thresh} ]]; then
        echo "gvt.py $permds $stat $thresh -r $permdir/rand_ind.txt $*"
    fi
done
