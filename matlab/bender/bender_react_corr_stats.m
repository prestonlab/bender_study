function con = bender_react_corr_stats(trial, correct, react_types, rois)
%BENDER_REACT_CORR_STATS   Test for correlation of reactivation between ROIs.
%
%  con = bender_react_corr_stats(trial, correct, react_types, rois)

n_subj = size(trial.(react_types{1}).(rois{1}), 1);

% results are sorted by group number, so category is in consistent
% order
category = ones(1, 120);
category(61:end) = 2;
cat_mat = repmat(category, [n_subj 1]);

con = struct;
for i = 1:length(react_types)
    f = react_types{i};
    % use robust regression to estimate slope of reactivation
    % between each pair of ROIs, separately for correct and
    % incorrect trials and different categories
    m_con = bender_roi_react_corr(trial.(f), rois, ...
                                  true(size(correct)), cat_mat);
    c_con = bender_roi_react_corr(trial.(f), rois, ...
                                  correct==1, cat_mat);
    i_con = bender_roi_react_corr(trial.(f), rois, ...
                                  correct==0, cat_mat);

    % test for significant connectivity/significant subsequent
    % memory effects
    con.(f).mean = roi_corr_stats(m_con, rois);
    con.(f).corr = roi_corr_stats(c_con, rois);
    con.(f).incorr = roi_corr_stats(i_con, rois);
    con.(f).sme = roi_corr_stats(c_con - i_con, rois);    
end
