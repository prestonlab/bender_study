#!/bin/bash
#
# Transform an image from MNI space to the custom template space.

if [ $# -lt 1 ]; then
    echo "Usage: bender_mni2template.sh [-n interp] input output"
    exit 1
fi

interp=Linear
while getopts ":n:" opt; do
    case $opt in
        n)
            interp=$OPTARG
            ;;
    esac
done
shift $((OPTIND-1))

input=$1
output=$2

tdir=$STUDYDIR/gptemplate/highres_brain_all

antsApplyTransforms -d 3 -i $input -r $tdir/gp_template_mni_affine.nii.gz -o $output -n $interp -t $tdir/gp_template2mni_1InverseWarp.nii.gz
