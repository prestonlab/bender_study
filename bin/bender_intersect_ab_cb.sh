#!/bin/bash
#
# Create intersection clusters for A>B and C>B contrasts.

if (( $# < 3 )); then
    echo "Usage: bender_intersect_ab_cb.sh ab_res_name cb_res_name intersect_res_name"
    exit 1
fi

ab_res_name=$1
cb_res_name=$2
intersect_res_name=$3

ab_dir=$STUDYDIR/batch/rsa/$ab_res_name
cb_dir=$STUDYDIR/batch/rsa/$cb_res_name
intersect_dir=$STUDYDIR/batch/rsa/$intersect_res_name

if [[ ! -d $ab_dir ]]; then
    echo "Error: AB results directory not found: $ab_dir"
fi

if [[ ! -d $cb_dir ]]; then
    echo "Error: CB results directory not found: $cb_dir"
fi

# hippocampus
for roi in b_hip b_prc; do
    mask=$STUDYDIR/batch/glm/study_stim2/$roi/mask.nii.gz
    if [[ ! -f $mask ]]; then
        echo "Error: $roi mask not found: $mask"
    fi

    res_dir=$intersect_dir/$roi
    mkdir -p "$res_dir"
    imcp "$mask" "$res_dir/mask"

    imcp "$ab_dir/stat_thresh" "$res_dir/ab"
    "$FSLDIR/bin/cluster" -i "$res_dir/ab" -t 0.001 --minextent=10 --oindex="$res_dir/ab_clusters10"
    if [[ $roi = b_hip ]]; then
        index=216
    else
        index=225
    fi
    fslmaths "$res_dir/ab_clusters10" -thr "$index" -uthr "$index" -bin "$res_dir/ab_cluster_mask"

    imcp "$cb_dir/stat_thresh" "$res_dir/cb"
    "$FSLDIR/bin/cluster" -i "$res_dir/cb" -t 0.001 --minextent=10 --oindex="$res_dir/cb_clusters10"
    if [[ $roi = b_hip ]]; then
        index=132
    else
        index=181
    fi
    fslmaths "$res_dir/cb_clusters10" -thr "$index" -uthr "$index" -bin "$res_dir/cb_cluster_mask"

    fslmaths "$res_dir/ab_cluster_mask" -mul "$res_dir/cb_cluster_mask" -bin "$res_dir/intersect_cluster_mask"
    fslmaths "$res_dir/ab_cluster_mask" -sub "$res_dir/cb_cluster_mask" -bin "$res_dir/ab_exclusive_cluster_mask"
    fslmaths "$res_dir/cb_cluster_mask" -sub "$res_dir/ab_cluster_mask" -bin "$res_dir/cb_exclusive_cluster_mask"

    fslmaths "$res_dir/ab_exclusive_cluster_mask" -mul 1 "$res_dir/temp1"
    fslmaths "$res_dir/cb_exclusive_cluster_mask" -mul 2 "$res_dir/temp2"
    fslmaths "$res_dir/intersect_cluster_mask" -mul 3 "$res_dir/temp3"
    fslmaths "$res_dir/temp1" -add "$res_dir/temp2" -add "$res_dir/temp3" "$res_dir/partitions"
    imrm "$res_dir"/temp{1,2,3}

    cp "$STUDYDIR"/gptemplate/highres_brain_all/gp_template_mni_affine.nii.gz "$res_dir"/bg_image.nii.gz
done
