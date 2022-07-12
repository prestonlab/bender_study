#!/bin/bash
#
# Pull betaseries files and masks.

src=$1
dest=$2
shift 2

rsync -azvu "$src" "$dest" \
    --include="batch/" \
    --include="batch/glm/" \
    --include="batch/glm/study_stim2/" \
    --include="batch/glm/study_stim2/rdm/" \
    --include="batch/glm/study_stim2/rdm/study_wiki_w2v_fix_cont_ac_bx_sme_*_dil1c/" \
    --include="bender_??_rdm.txt" \
    --prune-empty-dirs \
    --exclude="*" \
    "$@"
