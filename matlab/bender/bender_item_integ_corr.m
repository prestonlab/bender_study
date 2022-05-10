function b = bender_item_integ_corr(react_trial, rdm_roi, model_rdm, ...
                                    item_rois, integ_rois, correct)
%BENDER_ITEM_INTEG_CORR   Correlation between item reactivation and integration.
%
%  b = bender_item_integ_corr(react_trial, rdm_roi, model_rdm, item_rois,
%                             integ_rois, correct)

n_subj = size(correct, 1);
acc = [1 0];
b = NaN(n_subj, length(acc), length(item_rois), length(integ_rois));
for i = 1:length(acc)
    for j = 1:length(item_rois)
        for k = 1:length(integ_rois)
            % for each trial, what was the self correlation?
            r_self = react_trial.self.(item_rois{j});
            
            % what was the correlation with the AC model, for each
            % trial with the included accuracy code?
            r_ac = roi_model_corr_trial(rdm_roi.(integ_rois{k}), model_rdm, ...
                                        'subj_ind', correct==acc(i));

            % relationship between the two across included trials
            p = robust_reg_rows(r_self, r_ac, correct==acc(i));
            b(:,i,j,k) = p(:,2);
        end
    end
end
