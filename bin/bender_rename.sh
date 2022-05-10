#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage:   bender_rename.sh subject"
    echo "Example: bender_rename.sh 1 2"
    exit 1
fi

jobfile=$BATCHDIR/batch_rename.sh
if [ -e $jobfile ]; then
    rm $jobfile
fi

subjs=`subjids -s ' ' $@`
for id in $subjs; do
    for suffix in '' a b; do
	subject=${id}${suffix}
	if [ -d ${STUDYDIR}/${subject}/anatomy ]; then
	    echo $subject
	    echo "rename_nifti.py $subject" >> $jobfile
	fi
    done
done

cd $BATCHDIR
bash -x $jobfile
