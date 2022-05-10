function rdms = bender_load_rdms(beta_dir, analysis, roi, rdm_type, subjects)
%BENDER_LOAD_RDMS   Load individual RDMs for an ROI.
%
%  rdms = bender_load_rdms(beta_dir, analysis, roi, rdm_type, subjects)
%
%  INPUTS
%  beta_dir : string
%      Main directory with RDMs (e.g., 'study_stim2', 'models3')
%
%  analysis : string
%      Analysis name (e.g., 'cat_react_item2')
%
%  roi : string
%      Name of ROI mask (e.g., 'lprc')
%
%  rdm_type : string
%     RDM type ('rdm' 'react' 'model').
%
%  subjects : numeric vector
%     Subject numbers in behavioral format.
%
%  OUTPUTS
%  rdms : [items x items x subjects] numeric array
%      RDM for each subject.

if strcmp(rdm_type, 'model')
    rdm_dir = fullfile('~/work/bender/batch', beta_dir, 'trial');
else
    rdm_dir = fullfile('~/work/bender/batch/glm', beta_dir, rdm_type);
end

rdms = [];
for i = 1:length(subjects)
    subjid = sprintf('bender_%02d', subjects(i) - 100);
    if strcmp(rdm_type, 'model')
        rdm_file = fullfile(rdm_dir, ...
                            sprintf('%s_%s_%s.mat', ...
                                    subjid, analysis, roi));
        disp(rdm_file)
        rdm = getfield(load(rdm_file, 'rdm'), 'rdm');
    else
        mask_name = sprintf('%s_%s', analysis, roi);
        rdm_file = fullfile(rdm_dir, mask_name, sprintf('%s_rdm.txt', subjid));
        disp(rdm_file)
        rdm = load(rdm_file);
    end
    rdms = cat(3, rdms, rdm);
end
