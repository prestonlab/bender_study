#!/bin/bash

if [ $# -lt 4 ]; then
    echo "indiv_clust_mask.sh - Create individual masks from a group cluster"
    echo
    echo "Usage: indiv_clust_mask.sh [-a anat] [-r radius] -g subjects cluster_mask clustind maskname"
    echo
    echo "-a anat"
    echo "    reference anatomical image"
    echo
    echo "-r radius"
    echo "    radius for dilating individual-subject masks"
    echo
    echo "-g"
    echo "    if set, will intersect individual masks with a dilated"
    echo "    gray-matter mask"
    echo
    echo "subjects"
    echo "    colon-separated list of subject numbers to include"
    echo
    echo "cluster_mask"
    echo "    image with labeled clusters"
    echo
    echo "clustind"
    echo "    index of cluster to use for mask"
    echo
    echo "maskname"
    echo "    name of mask to create"
    echo
    exit 1
fi

radius=0
refanat=""
gray=false
while getopts ":r:a:g" opt; do
    case $opt in
        r)
            radius=$OPTARG
            ;;
        a)
            refanat=$OPTARG
            ;;
        g)
            gray=true
            ;;
    esac
done
shift $((OPTIND-1))

subjects=$1
cluster_mask=$2
clustind=$3
maskname=$4

if [ $(imtest $cluster_mask) = 0 ]; then
    echo "cluster mask does not exist: $cluster_mask"
    exit 1
fi

parent=$(dirname $cluster_mask)

# make group-level mask (saved to same directory as the cluster mask)
echo "Making group mask..."
group_mask=$parent/${maskname}.nii.gz
fslmaths $cluster_mask -thr $clustind -uthr $clustind -bin $group_mask
n_vox=$(fslstats $group_mask -V | cut -d ' ' -f 1)
echo "$maskname: $n_vox voxels"

# make individual masks
echo "Transforming to native spaces..."
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4
if [ -n "$refanat" ]; then
    flags="-a $refanat"
else
    flags=""
fi
parallel -P 6 bender_mni2func.sh "$flags" -n NearestNeighbor $group_mask $STUDYDIR/{}/anatomy/bbreg/data/${maskname}.nii.gz {} ::: $(subjids -s " " $subjects)

# expand the mask
if (( $(bc <<< "$radius > 0") )); then
    echo "Expanding native masks using a radius of $radius..."
    parallel -P 12 fslmaths $STUDYDIR/{}/anatomy/bbreg/data/${maskname}.nii.gz -kernel sphere $radius -dilD $STUDYDIR/{}/anatomy/bbreg/data/${maskname}.nii.gz ::: $(subjids -s " " $subjects)
fi

# intersect with a minimally dilated gray matter mask
# made using: fslmaths b_gray -kernel sphere 1.75 -dilD
if [ $gray = true ]; then
    echo "Intersecting with gray matter mask..."
    parallel -P 12 fslmaths $STUDYDIR/{}/anatomy/bbreg/data/${maskname}.nii.gz -mas $STUDYDIR/{}/anatomy/bbreg/data/b_gray_dil3mm.nii.gz $STUDYDIR/{}/anatomy/bbreg/data/${maskname}.nii.gz ::: $(subjids -s " " $subjects)
fi
