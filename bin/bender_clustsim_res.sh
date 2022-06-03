#!/bin/bash
#
# Calculate average smoothness and estimate null cluster distribution.

if [[ $# -lt 4 ]]; then
    echo "Calculate average smoothness and estimate null cluster distribution."
    echo
    echo "Usage:   bender_clustsim_res.sh [-n nthreads] subjects runs model mask"
    echo 'Example: bender_clustsim_res.sh $SUBJNOS $STUDYRUNS study_stim2 $WORK/bender/gptemplate/highres_brain_all/mpfc.nii.gz'
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

subjects=$1
runs=$2
model=$3
mask=$4

if [[ ! -e $mask ]]; then
    echo "Error: mask does not exist: $mask"
    exit 1
fi

maskname=$(basename "$mask" .nii.gz)
outdir=$STUDYDIR/batch/glm/$model/$maskname
mkdir -p "$outdir"

# get list of ACF files to concatenate
reslist=""
for subject in $(subjids -s ' ' "$subjects"); do
    for run in $(echo "$runs" | tr ':' ' '); do
        # original file
        mask_dir=$STUDYDIR/$subject/model/$model/${run}.feat/stats/$maskname
        acf_file=$mask_dir/acf_smoothness
        if [[ ! -e $acf_file ]]; then
            echo "Error: acf file missing: $acf_file"
            exit 1
        fi

        # create trimmed file with just the ACF parameters (no
        # Gaussian parameters)
        acf_trim=$mask_dir/acf_smoothness_trim
        tail -n 1 < "$acf_file" > "$acf_trim"
        if [[ -z $reslist ]]; then
            reslist="$acf_trim"
        else
            reslist="$reslist $acf_trim"
        fi
    done
done

# concatenate ACF estimates from all runs
cat $reslist > "$outdir"/acf_smoothness

# extract parameters in correct format for 3dClustSim
acfpar=$(bender_res_mean.py "$outdir"/acf_smoothness | awk '{print $1,$2,$3}')
echo "ACF parameters: $acfpar"
echo "$acfpar" > "$outdir"/acf_par

# estimate null distribution of cluster sizes within this mask
export OMP_NUM_THREADS=$nthreads
cp "$mask" "$outdir"/mask.nii.gz
cd "$outdir" || exit
3dClustSim -mask mask.nii.gz -acf "$acfpar" -iter 2000 -nodec -prefix clustsim
