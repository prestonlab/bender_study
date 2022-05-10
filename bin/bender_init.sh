#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage:   bender_init.sh subject"
    echo "Example: bender_init.sh bender_01"
    exit 1
fi

subject=$1

echo "Converting dicom files:"
echo "Day 1..."
convert_dicom.py $subject
echo "Day 2..."
convert_dicom.py ${subject}a
if [ -d ${STUDYDIR}/${subject}b ]; then
    echo "Day 2 part 2..."
    convert_dicom.py ${subject}b
fi

echo "Renaming nifti files:"
echo "Day 1..."
rename_nifti.py $subject
echo "Day 2..."
rename_nifti.py ${subject}a
if [ -d ${STUDYDIR}/${subject}b ]; then
    echo "Day 2 part 2..."
    rename_nifti.py ${subject}b
fi
