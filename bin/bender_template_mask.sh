#!/bin/bash

maskname=$1
subjects=$2

datalist=""
for subject in $(subjids -s ' ' $subjects); do
    file=$WORK/bender/$subject/anatomy/antsreg/data/$maskname.nii.gz
    if [ -z "$datalist" ]; then
	    datalist="$file"
    else
	    datalist="$datalist $file"
    fi
done

cd $WORK/bender/gptemplate/highres_brain_all

fslmerge -t ${maskname}_all $datalist
fslmaths ${maskname}_all -Tmean ${maskname}_mean
fslmaths ${maskname}_mean -thr 0.9 -kernel sphere 3 -dilD -bin -mas b_gray $maskname
