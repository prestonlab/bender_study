#!/bin/bash
#
# Estimate smoothness based on standard-space residuals.

if [[ $# -lt 4 ]]; then
    echo "Estimate smoothness based on standard-space residuals."
    echo
    echo "Usage:   bender_res_smoothness.sh [-n nthreads] subject run model mask"
    echo "Example: bender_res_smoothness.sh bender_02 study_1 study_stim2 $WORK/bender/gptemplate/highres_brain_all/mpfc.nii.gz 24"
    echo
    exit 1
fi

nthreads=1
while getopts ":n:" opt; do
    case $opt in
    n)
        nthreads=$OPTARG
        ;;
    *)
        echo "Unknown option: ${opt}"
        exit 1
        ;;
    esac
done
shift $((OPTIND-1))

subject=$1
run=$2
model=$3
mask=$4

model_dir=$SCRATCH/$STUDY/$subject/model/$model/${run}.feat/stats

if [[ ! -e $model_dir/res4d_std.nii.gz ]]; then
    echo "Error: res4d_std image not found."
    exit 1
fi

# make a directory for calculations within this mask. Put on work,
# since this will take up very little space
maskname=$(basename "${mask}" .nii.gz)
outdir=$STUDYDIR/$subject/model/$model/${run}.feat/stats/$maskname

mkdir -p "$outdir"
rm -f "$outdir"/acf*
cp "$mask" "$outdir"/mask.nii.gz

# run smoothness estimation based on the residuals in standard space
cd "$outdir" || exit
export OMP_NUM_THREADS=$nthreads
3dFWHMx -mask mask.nii.gz -acf acf -input "$model_dir"/res4d_std.nii.gz -out acf_vol -arith > acf_smoothness
