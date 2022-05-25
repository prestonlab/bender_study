#!/bin/bash
#
# Run group volume threshold analysis on permutation test results.

if [[ $# -lt 2 ]]; then
    echo "Usage: bender_run_gvt.sh [-m mask] [-i interp] [-a alpha] [-p nperm] [-A project] filepath subjects"
    exit 1
fi

alpha=0.01
interp=BSpline
mask=${STUDYDIR}/gptemplate/highres_brain_all/gp_template_mni_affine_mask.nii.gz
n_subj_perm=100
dry_run=false
project=ANTS
while getopts ":a:i:m:p:A:d" opt; do
    case $opt in
    a)
        alpha=$OPTARG
        ;;
    i)
        interp=$OPTARG
        ;;
    m)
        mask=$OPTARG
        ;;
    p)
        n_subj_perm=$OPTARG
        ;;
    A)
        project=$OPTARG
        ;;
    d)
        dry_run=true
        ;;
    *)
        echo "Unknown option: ${opt}"
        exit 1
        ;;
    esac
done
shift $((OPTIND-1))

filepath=$1
subjects=$2
analysis=$(basename "${filepath}")
log_file=${BATCHDIR}/gvt_${analysis}.log

echo "Options:"
echo "name:   ${analysis}"
echo "alpha:  ${alpha}"
echo "nperm:  ${n_subj_perm}"
echo "interp: ${interp}"
echo "mask:   ${mask}"
echo "log:    ${log_file}"

maxperm=$(printf '%03d' "${n_subj_perm}")
for subject in $(subjids -s ' ' "${subjects}"); do
    if [[ ! -f ${STUDYDIR}/${subject}/${filepath}/perm${maxperm}.nii.gz ]]; then
        echo "Error: missing permutation file for ${subject}."
        exit 1
    fi
done

if [[ -f ${log_file} ]]; then
    rm "${log_file}"
fi

# merge permutations, transform images to template space, split into
# slices
if [[ $dry_run = true ]]; then
    echo "bender_perm2mni.sh -m ${mask} -i ${interp} {} ${filepath}"
else
    jobid1=$(slaunch -J gvt_perm2mni "bender_perm2mni.sh -m $mask -i $interp {} $filepath" "$subjects" -N 13 -n 26 -a 12 -r 01:00:00 -A "$project" | tee -a "$log_file" | getjid)
fi

# prep gvt directory, split mask for creating slice file datasets
if [[ $dry_run = true ]]; then
    echo "bender_prep_gvt.sh -m $mask $filepath $subjects"
else
    jobid2=$(ezlaunch "bender_prep_gvt.sh -m $mask $filepath $subjects" -N 1 -n 1 -r 00:05:00 -A "$project" -d "$jobid1" | tee -a "$log_file" | getjid)
fi

# merge slice files into datasets with all subjects and all permutations
prep_file=$BATCHDIR/prep_gvt_${analysis}.sh
bender_prep_gvt_commands.sh "$filepath" "$subjects" > "$prep_file"
if [[ $dry_run = true ]]; then
    echo "ezlaunch -s $prep_file"
else
    jobid3=$(ezlaunch -s "$prep_file" -N 1 -n 24 -r 00:30:00 -A "$project" -d "$jobid2" | tee -a "$log_file" | getjid)
fi

# run gvt on each slice
run_file=$BATCHDIR/run_gvt_${analysis}.sh
perm_dir=$SCRATCH/$STUDY/batch/$filepath
mkdir -p "$perm_dir"
n_subj=$(echo "$subjects" | tr ':' ' ' | wc -w)
n_perm=100000
bender_gen_rand_ind.py "$n_subj" "$n_subj_perm" "$n_perm" > "$perm_dir"/rand_ind.txt
bender_gvt_commands.sh "$perm_dir" -a "$alpha" > "$run_file"
if [[ $dry_run = true ]]; then
    echo "ezlaunch -s $run_file"
else
    jobid4=$(ezlaunch -s "$run_file" -N 10 -n 10 -r 03:00:00 -A "$project" -d "$jobid3" | tee -a "$log_file" | getjid)
fi

# merge slices, calculate thresholded image, backup main results
if [[ $dry_run = true ]]; then
    echo "bender_post_gvt.sh $filepath"
else
    ezlaunch "bender_post_gvt.sh $filepath" -N 1 -n 1 -r 00:10:00 -A "$project" -d "$jobid4" | tee -a "$log_file"
fi
