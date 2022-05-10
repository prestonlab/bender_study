#!/bin/bash
#
# Create masks based on significant clusters.

# roi name, relative path to analysis directory, cluster index in the
# standard cluster100 image
roi_spec=$WORK/bender/batch/sl_roi_spec.txt
cat <<EOF > $roi_spec
lprc,mvpa/cat_react_item2,43
rphc,mvpa/cat_react_item2,52
rifg,mvpa/cat_react_item2,53
phpc,mvpa/item_suppress_gvt,32
ampfc,mvpa/cat_react_item_sme2,58
rifg,mvpa/cat_react_item_sme2,61
mmpfc,rsa/study_wiki_w2v_fix_cont_a_bc_sme,56
rhpc,rsa/study_wiki_w2v_fix_cont_ac_bx_sme,60
rprc,rsa/study_wiki_w2v_fix_cont_ac_bx_sme,55
EOF

while read entry; do
    echo "${entry}"

    # unpack spec
    roi=$(echo "${entry}" | cut -d , -f 1)
    filepath=$(echo "${entry}" | cut -d , -f 2)
    cluster=$(echo "${entry}" | cut -d , -f 3)

    statdir=$STUDYDIR/batch/$filepath
    analysis=$(basename "${filepath}")
    maskname=${analysis}_${roi}_dil1nn
    cluster -i "${statdir}/stat_thresh" -t 0.001 --minextent=100 -o "${statdir}/cluster_mask100" > "${statdir}/cluster100.txt"
    bender_clust_mask.sh -a 2 -r 1.75 "${SUBJNOS}" "${statdir}/cluster_mask100.nii.gz" "${cluster}" "${maskname}"
done < $roi_spec
