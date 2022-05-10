function print_item_cat_react_sme(stats, rois, fig_file)
%PRINT_ITEM_CAT_REACT_SME   Plot item and category reactivation SMEs by ROI.
%
%  print_item_cat_react_sme(stats, rois, fig_file)

n_subj = length(stats.mean.(rois{1}).item);
mat = NaN(length(rois), 4, n_subj);
for i = 1:length(rois)
    s_corr = stats.corr.(rois{i});
    s_incorr = stats.incorr.(rois{i});
    
    % everyone has at least one incorrect item, but not always one
    % incorrect item per category
    corr_item = nanmean(s_corr.item_cat - s_corr.within_cat, 2);
    incorr_item = nanmean(s_incorr.item_cat - s_incorr.within_cat, 2);
    corr_cat = nanmean(s_corr.within_cat, 2) - s_corr.between;
    incorr_cat = nanmean(s_incorr.within_cat, 2) - s_incorr.between;
    
    mat(i,:,:) = [corr_item incorr_item corr_cat incorr_cat]';
end

disp('ROIs:')
disp(rois)
disp('Stats:')
disp({'item corr' 'item incorr' 'cat corr' 'cat incorr'})
disp('reactivation/suppression:')
ttest_tab(mat, 0, 'dim', 3);
disp('item reactivation SME:')
ttest_tab(mat(:,1,:), mat(:,2,:), 'dim', 3);
disp('category reactivation SME:')
ttest_tab(mat(:,3,:), mat(:,4,:), 'dim', 3);

clf
y = nanmean(mat, 3);
err = nanstd(mat, [], 3) ./ sqrt(sum(~isnan(mat), 3) - 1);
colors = [.3255 .5608 .4980   % item reactivation, correct
          .5569 .7529 .7255   % item reactivation, incorrect
          .3843 .2863 .4784   % category reactivation, correct
          .6431 .5373 .6902]; % category reactivation, incorrect
[hbar, herr] = ebar([], y, err);
for i = 1:length(hbar)
    hbar(i).FaceColor = colors(i,:);
    hbar(i).LineStyle = 'none';
    herr(i).CapSize = 10;
    herr(i).LineWidth = 1;
end
set(gca, 'XTick', 1:length(rois), 'XLim', [0.5 length(rois)+.5], ...
         'XTickLabel', upper(rois))
%l = legend({'corr. item' 'incorr. item' 'corr. category' 'incorr. category'});
%l.Location = 'NorthWest';
ylabel('reactivation')
set_fig_style(gcf, 'preston')
legend off
set(gca, 'YLim', [-.04 .04], 'YTick', -.04:.02:.04)

if ~isempty(fig_file)
    print(gcf, '-depsc', fig_file);
end
