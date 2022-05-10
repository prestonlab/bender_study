function rho = roi_model_corr_trial(rdm, model, varargin)
%ROI_MODEL_CORR_TRIAL   Calculate trial-level model correlations.
%
%  For each subject and trial, calculate correlations between ROI
%  RDM rows and model RDM rows. Can specify trial subsets to
%  include in the row correlations.
%
%  rho = roi_model_corr_trial(model_type, model_name, roi_type,
%      roi_name, subjects, ...)
%
%  OUTPUTS
%  rho: [subjects x items] matrix
%      correlations between ROI RDM rows and model RDM rows (not
%      including each trial's self-similarity)

n_subj = size(rdm, 3);
n_item = size(rdm, 1);

def.subj_ind = true(n_subj, n_item);
opt = propval(varargin, def);

rho = NaN(n_subj, n_item);
for i = 1:n_subj
    for j = 1:n_item
        % items to include in calculating correlation
        include = opt.subj_ind(i,:);
        
        % remove the self-similarity, which is non-informative
        include(j) = false;
        
        % get correlations with other items
        model_vec = model(include,j,i);
        roi_vec = rdm(include,j,i);
        
        rho(i,j) = corr(model_vec, roi_vec, 'Rows', 'Complete', ...
                        'Type', 'Spearman');
    end
end
