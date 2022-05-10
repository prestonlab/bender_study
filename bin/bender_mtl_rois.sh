#!/bin/bash
#
# Transform drawn medial temporal masks to the template space.

mnidir=$STUDYDIR/gptemplate/tofu_mni
tdir=$STUDYDIR/gptemplate/highres_brain_all

mkdir -p $tdir/roi/mtl
if [ ! -e $tdir/drawn_mtl.nii.gz ]; then
    # project labels into template space, based on template-MNI
    # nonlinear registration
    bender_mni2template.sh -n MultiLabel $mnidir/drawn_mtl.nii.gz $tdir/drawn_mtl.nii.gz
fi

imcp $tdir/drawn_mtl $tdir/roi/mtl
imcp $tdir/b_gray $tdir/roi/mtl

# take HPC ROI from the FS segmentation data. Included here just to
# make sure the other ROIs don't run into HPC
imcp $tdir/b_hip_fshs_dil1 $tdir/roi/mtl/hip
cd $tdir/roi/mtl

# get bilateral mask for each ROI
mergelabels prc drawn_mtl 1 2
mergelabels erc drawn_mtl 3 4
#mergelabels hip drawn_mtl 5 6
mergelabels phc drawn_mtl 7 8

# smooth to prep for max labels calculation
sigma=$(python -c 'print(2/2.355)')
fslmaths prc -uthr 0 -bin empty
files="empty"
for roi in hip phc prc erc; do
    echo "expanding ${roi}..."
    fslmaths $roi -s $sigma ${roi}_sm

    files="$files ${roi}_sm"
done

# get the best label for each voxel
echo "calculating max labels..."
fslmerge -t mtl_all $files
fslmaths mtl_all -Tmaxn mtl_max_ind

# within closest voxels for ROI, expand mask and intersect with grey
# matter. Don't include hip here, as it's already finalized
i=1
for roi in phc prc erc; do
    echo "finalizing ${roi}..."
    i=$((i+1))
    fslmaths mtl_max_ind -thr $i -uthr $i -bin ${roi}_max
    fslmaths ${roi} -kernel sphere 3.1 -dilD ${roi}_dil
    fslmaths ${roi}_dil -mas ${roi}_max ${roi}_max_dil
    fslmaths ${roi}_max_dil -mas b_gray ${roi}_max_dil_mask
    imcp ${roi}_max_dil_mask $tdir/b_${roi}2
done
