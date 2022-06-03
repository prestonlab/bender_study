#!/bin/bash
#
# Run small-volume correction on GVT results.

if [[ $# -lt 2 ]]; then
    echo "Run small-volume correction on GVT results."
    echo
    echo "bender_svc.sh [-o] filepath model mask1 mask2 ..."
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

filepath=$1
model=$2

shift 2

recalc=$overwrite
for mask; do
    echo "$mask"
    table=$STUDYDIR/batch/$filepath/$mask/cluster_cope1_std.txt
    if [[ ! -e $table || $overwrite = true ]]; then
        if [ $recalc = true ]; then
            # must redo stat image thresholding
            bender_apply_clustsim_res.sh -o "$filepath" "$model" "$mask"
            recalc=false
        else
            # just redo the mask-specific stats
            bender_apply_clustsim_res.sh "$filepath" "$model" "$mask"
        fi
    else
        cat "$table"
    fi
done
