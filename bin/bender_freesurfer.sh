#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Usage: run_freesurfer.sh subject nthreads"
    exit 1
fi

if [ -u $STUDYDIR ]; then
    echo "STUDYDIR unset; quitting."
    exit 1
fi

if [ ! -d $STUDYDIR ]; then
    echo "STUDYDIR does not exist; quitting."
    exit 1
fi

subject=$1
nthreads=$2
subjdir=$STUDYDIR/$subject

if [ ! -f $subjdir/anatomy/highres.nii.gz ]; then
    echo "ERROR: Highres file not found."
    exit 1
fi

if [ ! -f $subjdir/anatomy/coronal.nii.gz ]; then
    echo "ERROR: Coronal file not found."
    exit 1
fi

source $FREESURFER_HOME/SetUpFreeSurfer.sh
# parallel flag indicates that multiple threads should be used when
# possible. openmp and itkthreads flags set the number of threads for
# different applications (probably openmp for most FS tools,
# itkthreads for a few things borrowed from ITK or ANTS). Using 12
# threads instead of 24 because hemispheres are also sometimes run in
# parallel
recon-all -s ${subject} -sd $subjdir/anatomy/ \
	  -i $subjdir/anatomy/highres.nii.gz -all \
	  -hippocampal-subfields-T2 $subjdir/anatomy/coronal.nii.gz T2 \
	  -parallel -openmp $nthreads -itkthreads $nthreads
