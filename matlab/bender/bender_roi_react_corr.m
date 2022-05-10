function b = bender_roi_react_corr(item, rois, include, category)
%BENDER_ROI_REACT_CORR   Test whether ROI connectivity predicts performance.
%
%  b = bender_roi_react_corr(item, rois, include, category)

ucat = unique(category);
n_cat = length(ucat);
n_subj = size(include, 1);
b = NaN(n_subj, length(rois), length(rois));
for i = 1:length(rois)
    item1 = item.(rois{i});
    for j = i+1:length(rois)
        item2 = item.(rois{j});
        p = NaN(n_subj, n_cat);
        for k = 1:n_cat
            % use robust regression to calculate relationship
            % between regions, for each subject
            pc = robust_reg_rows(item1, item2, ...
                                 include & category == ucat(k));

            % slope is second coefficient
            p(:,k) = pc(:,2);
        end
        subj_b = mean(p, 2);
        b(:,i,j) = subj_b;
    end
end
