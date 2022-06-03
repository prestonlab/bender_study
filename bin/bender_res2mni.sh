#!/bin/bash
#
# Transforms residual images to template space.

if [[ $# -lt 3 ]]; then
    echo "Transform residual images to template space."
    echo
    echo "Usage: bender_res2mni.sh subject runid model"
    echo
    exit 1
fi

subject=$1
runid=$2
model=$3
shift 3

# directory with level 1 model results
model_dir=$SCRATCH/$STUDY/$subject/model/$model/${runid}.feat
if [[ ! -d $model_dir ]]; then
    echo "Error: model directory not found: $model_dir" 1>&2
    exit 1
fi

# residuals image
input=$model_dir/stats/res4d.nii.gz
if [[ ! -s $input ]]; then
    echo "Error: residuals not found: $input" 1>&2
    exit 1
fi

# transform to standard space
output=$model_dir/stats/res4d_std.nii.gz
bender_func2mni.sh "$@" "${input}" "${output}" "${subject}"
