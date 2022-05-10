function [trial, stats] = bender_roi_react_stats(react_roi, correct, rois)
%BENDER_ROI_REACT_STATS    Calculate reactivation stats in ROIs.
%
%  [trial, stats] = bender_roi_react_stats(react_roi, correct)
%
%  INPUTS
%  react_roi - struct
%      Struct with a field for each ROI. Each must contain an
%      [items x items x subjects] matrix of dissimilarity values.
%
%  correct - [subjects x trials] numeric array
%      For each trial, correct (1) or incorrect (0).
%
%  OUTPUTS
%  trial - struct
%      Trial-level self, within, and between-category similarity.
%
%  stats - struct
%      Average similarities over all trials and within just correct
%      or incorrect trials.

if nargin < 3
    rois = fieldnames(react_roi);
end
    
n_subject = size(correct, 1);
n_item = 120;
category = ones(1, 120);
category(61:end) = 2;

self_mat = eye(n_item) == 1;
within_mat = category == category';

trial = struct;
stats = struct('mean', struct, 'corr', struct, 'incorr', struct);
for i = 1:length(rois)
    f = rois{i};
    if isfield(trial, f)
        f = [f '2'];
    end
    
    % reactivation stats for different subsets of trials
    rdms = react_roi.(rois{i});
    z = atanh(1 - rdms);
    stats.mean.(f) = react_rdm_stats(z, category);
    stats.corr.(f) = react_rdm_stats(z, category, ...
                                     'subj_ind', correct==1);
    stats.incorr.(f) = react_rdm_stats(z, category, ...
                                       'subj_ind', correct==0);

    % trial-level self-similarity, within- and between-category similarities
    within = NaN(n_subject, n_item);
    between = NaN(n_subject, n_item);
    for j = 1:n_subject
        for k = 1:n_item
            % similarity between the study pattern for this trial
            % and each item's pre-exposure pattern
            trial_z = z(:,k,j);
            trial_self = self_mat(:,k);
            trial_within = within_mat(:,k);
            
            self(j,k) = trial_z(trial_self);
            within(j,k) = mean(trial_z(~trial_self & trial_within));
            between(j,k) = mean(trial_z(~trial_within));
        end
    end
    trial.self.(f) = self;
    trial.within.(f) = within;
    trial.between.(f) = between;
    trial.item.(f) = self - within;
    trial.category.(f) = within - between;
end
