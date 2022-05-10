#!/bin/bash
#
# Merge inferior frontal gyrus and insula masks.

subject=$1

anatdir=$WORK/bender/$subject/anatomy/antsreg/data
fslmaths $anatdir/b_ifg -add $anatdir/b_insula -bin $anatdir/b_ifg_insula

anatdir=$WORK/bender/$subject/anatomy/bbreg/data
fslmaths $anatdir/b_ifg -add $anatdir/b_insula -bin $anatdir/b_ifg_insula
