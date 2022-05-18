#!/bin/bash
#
# Apply smoothing using FSL SUSAN as in FEAT 6.0.

if [[ $# -lt 4 ]]; then
    cat <<EOF
Usage: bender_smooth_susan [-f hpf_sigma] [-vk] input mask fwhm output
Example: bender_smooth_susan -f 32 loc_1 mask 4.0 loc_1_sm

Options:
-f
    Highpass filter (sigma in volumes)

-v
    Print detailed information

-k
    Keep intermediate files
EOF
    exit 1
fi

verbose=0
keep=0
highpass=""
while getopts ":vkf:" opt; do
    case $opt in
	v)
	    verbose=1
	    ;;
	k)
	    keep=1
	    ;;
	f)
	    highpass=$OPTARG
	    ;;
	*)
	    echo "Invalid option: $opt"
	    exit 1
    esac
done
shift $((OPTIND-1))

func=$1
mask=$2
smooth=$3
output=$4

# create temporary directory
name=$(basename "$func" .nii.gz)
parent=$(dirname "$func")
dirname=${name}.susan
while [[ -d $parent/$dirname ]]; do
    dirname=${dirname}+
done
sdir=$parent/$dirname
mkdir -p "$sdir"

# copy mask
if [[ $(imtest "$mask") = 0 ]]; then
    echo "Error: mask does not exist: $mask"
    exit 1
fi

fslmaths "$func" "$sdir/prefiltered_func_data" -odt float
imcp "$mask" "$sdir/mask"

pd=$(pwd)
cd "$sdir" || exit

# prepare functional data
funcdata=prefiltered_func_data
median_intensity=$(fslstats $funcdata -k mask -p 50)
fslmaths $funcdata -mas mask prefiltered_func_data_thresh
funcdata=prefiltered_func_data_thresh

# spatial filtering
smoothsigma=$(python3 -c "print(${smooth} / 2.355)")

# the 0.1 percentile within mask is similar to using the 2nd
# percentile over the whole image, but better controlled
int1=$(fslstats $funcdata -k mask -p 0.1)
susan_int=$(python3 -c "print((${median_intensity} - ${int1}) * 0.75)")
if [[ $verbose == 1 ]]; then
    echo "int1: $int1"
    echo "median: $median_intensity"
    echo "susan int: $susan_int"
    echo "smooth sigma: $smoothsigma"
fi
fslmaths "$funcdata" -Tmean mean_func
susan "$funcdata" "$susan_int" "$smoothsigma" 3 1 1 mean_func "$susan_int" \
      prefiltered_func_data_smooth

# optional temporal filtering
if [[ -n $highpass ]]; then
    fslmaths prefiltered_func_data_smooth -Tmean tempMean
    fslmaths prefiltered_func_data_smooth -bptf "$highpass" -1 -add tempMean \
        prefiltered_func_data_tempfilt
    cd "$pd" || exit
    fslmaths "$sdir/prefiltered_func_data_tempfilt" -mas "$sdir/mask" "$output"
else
    cd "$pd" || exit
    fslmaths "$sdir/prefiltered_func_data_smooth" -mas "$sdir/mask" "$output"
fi

# delete temporary files
if [[ $keep == 0 ]]; then
    cd "$parent" && rm -rf "$dirname"
    cd "$pd" || exit
fi
