#!/usr/bin/env python

from bender_study.subjutil import *

parser = SubjParser()
args = parser.parse_args()

sp = SubjPath(args.subject, args.study_dir)
log = sp.init_log('regdays', 'preproc', args)
log.start()

img_dir = sp.path('anatomy')
xfm_dir = sp.path('anatomy', 'antsreg', 'transforms')

log.run('mkdir -p %s' % xfm_dir)

img_names = ['highres','coronal']
    
for img_name in img_names:
    # get the base images
    names = []
    for i in range(2):
        name = '%s%d' % (img_name, i + 1)
        names.append(name)
    image1 = sp.image_path('anatomy', names[0])
    image2 = sp.image_path('anatomy', names[1])
    merged = sp.image_path('anatomy', img_name)
    mask = sp.image_path('anatomy', img_name + '_brain')
    if not os.path.exists(image2) and os.path.exists(image1):
        merge = False
    else:
        merge = True
        
    if img_name == 'coronal':
        if merge:
            log.run('merge_anat.sh -c %s %s %s' % (
                image1, image2, merged))
        else:
            log.run('N4BiasFieldCorrection -i {} -o {}'.format(
                image1, merged))
        log.run('bet {} {} -f 0.01'.format(merged, mask))
    else:
        if merge:
            log.run('merge_anat.sh %s %s %s' % (
                image1, image2, merged))
        else:
            log.run('N4BiasFieldCorrection -i {} -o {}'.format(
                image1, merged))
        log.run('bet {} {}'.format(merged, mask))

# register coronal to highres
highres_merge = sp.image_path('anatomy', 'highres_brain')
coronal_merge = sp.image_path('anatomy', 'coronal_brain')
xfm_base = os.path.join(xfm_dir, 'coronal-highres_')
log.run('antsRegistration -d 3 -r [{highres},{coronal},1] -t Rigid[0.1] -m MI[{highres},{coronal},1,32,Regular,0.25] -c [1000x500x250x100,1e-6,10] -f 8x4x2x1 -s 3x2x1x0vox -n BSpline -w [0.005,0.995] -o {xfm}'.format(
    highres=highres_merge, coronal=coronal_merge, xfm=xfm_base))
c2h = os.path.join(xfm_dir, xfm_base + '0GenericAffine.mat')

# create registered coronal image
coronal_reg = sp.image_path('anatomy', 'coronal-highres')
if not merge:
    # only one coronal; just carry out the transformation
    log.run('antsApplyTransforms -i {input} -o {output} -r {ref} -t {c2h} -n BSpline'.format(
        input=coronal_merge, output=coronal_reg, ref=highres_merge, c2h=c2h))
else:
    # transform coronal 2 (coronal 2 to highres)
    cmerge_dir = sp.path('anatomy', 'coronal1_coronal2')
    coronal2 = impath(cmerge_dir, 'fix_cor')
    coronal2_reg = sp.image_path('anatomy', 'coronal2-highres')
    log.run('antsApplyTransforms -i {input} -o {output} -r {ref} -t {c2h} -n BSpline'.format(
        input=coronal2, output=coronal2_reg, ref=highres_merge, c2h=c2h))

    # transform coronal 1 (coronal 1 to coronal 2, coronal 2 to highres)
    coronal1 = impath(cmerge_dir, 'mov_cor')
    coronal1_reg = sp.image_path('anatomy', 'coronal1-highres')
    c2c = os.path.join(cmerge_dir, 'mov-fix_0GenericAffine.mat')
    log.run('antsApplyTransforms -i {input} -o {output} -r {ref} -t {c2h} -t {c2c} -n BSpline'.format(
        input=coronal1, output=coronal1_reg, ref=highres_merge,
        c2h=c2h, c2c=c2c))

    # merge registered coronals
    log.run('merge_anat.sh -m {} {} {}'.format(
        coronal1_reg, coronal2_reg, coronal_reg))

log.finish()
