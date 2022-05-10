function roi = bender_roi_rdms(rdm_type, roi_names, roi_types, fnames, ...
                               beta_name, subjects, suffix)
%BENDER_ROI_RDMS   Load reactivation or within-phase RDMs.
%
%  roi = bender_roi_rdms(rdm_type, roi_names, roi_types, fnames, 
%                        beta_name, subjects)
%
%  INPUTS
%  rdm_type - char
%      Type of RDM to load ('react' or 'rdm').
%
%  roi_names - cell array of strings
%      Name of each ROI to analyze.
%
%  roi_types - cell array of strings
%      Name of the analysis that generated each ROI.
%
%  fnames - cell array of strings
%      Name of the field of each ROI in the output struct.
%
%  beta_name - string
%      Name of the betaseries used for analysis.
%
%  subjects - numeric array
%      All subject numbers.
%
%  OUTPUTS
%  roi - struct
%      Struct with a field for each ROI. Each will contain an
%      [items x items x subjects] matrix of dissimilarity values.

roi = struct();
for i = 1:length(roi_names)
    rdms = bender_load_rdms(beta_name, roi_types{i}, ...
                            [roi_names{i} suffix], rdm_type, subjects);
    roi.(fnames{i}) = rdms(1:120,1:120,:);
end
