#!/bin/bash
#
# Transform searchlight results to template space.

if [[ $# -lt 2 ]]; then
    echo "Transform subject permuted searchlight results to template space."
    echo
    echo "Usage:   bender_perm2mni.sh [-m mask] [-i interp] subject filepath"
    echo "Example: bender_perm2mni.sh bender_02 mvpa/cat_react_sim"
    echo
    echo "If no mask specified, will use the bender template mask."
    echo
    exit 1
fi

interp=BSpline
mask=$STUDYDIR/gptemplate/highres_brain_all/gp_template_mni_affine_mask.nii.gz
while getopts ":i:m:" opt; do
    case $opt in
    i)
        interp=$OPTARG
        ;;
    m)
        mask=$OPTARG
        ;;
    *)
        echo "Error: unknown option $opt."
        exit 1
        ;;
    esac
done
shift $((OPTIND - 1))

subject=$1
filepath=$2
sdir=${STUDYDIR}/${subject}

# reference image
template=$STUDYDIR/gptemplate/highres_brain_all/gp_template_mni_affine.nii.gz

# transforms
orig2template_warp=$sdir/anatomy/antsreg/transforms/orig-template_Warp.nii.gz
orig2template=$sdir/anatomy/antsreg/transforms/orig-template_Affine.txt
orig2refvol=$sdir/anatomy/bbreg/transforms/highres-refvol

# save transformed data to scratch to save space on work
indir=$STUDYDIR/$subject/$filepath
outdir=$SCRATCH/$STUDY/$subject/$filepath
mkdir -p "${outdir}"

# concatenate first; runs faster transforming as a timeseries rather
# than as individual volumes
echo "Concatenating permutations..."
fslmerge -t "${outdir}"/allperm "${indir}"/perm*.nii.gz

echo "Transforming to template space..."
if [[ ! -f ${orig2refvol}.txt ]]; then
    if [[ ! -f ${orig2refvol}.mat ]]; then
        echo "Error: anatomical to functional registration missing."
        exit 1
    fi

    refvol=$sdir/BOLD/antsreg/data/refvol.nii.gz
    if [[ -f $sdir/anatomy/orig2.nii.gz ]]; then
        src=$sdir/anatomy/orig2.nii.gz
    else
        src=$sdir/anatomy/orig.nii.gz
    fi
    c3d_affine_tool -ref "${refvol}" -src "${src}" "${orig2refvol}.mat" -fsl2ras -oitk "${orig2refvol}.txt"
fi

antsApplyTransforms -d 3 -e 3 -i "${outdir}/allperm.nii.gz" -o "${outdir}/allperm_std.nii.gz" -r "${template}" -n "${interp}" -t "${orig2template_warp}" -t "${orig2template}" -t ["${orig2refvol}.txt",1]

# mask to remove any small values outside the mask caused by interpolation
fslmaths "${outdir}/allperm_std" -mas "${mask}" "${outdir}/allperm_std"

# create slice files. This will be important when running GVT, since
# it takes a huge amount of memory
echo "Splitting template space file..."
fslsplit "${outdir}"/allperm_std "${outdir}"/allperm_std_z -z
