function print_item_react_sme(stats, rois, fig_file)
%PRINT_ITEM_REACT_SME   Plot item and category reactivation SMEs by ROI.
%
%  print_item_react_sme(stats, rois, fig_file)

n_subj = length(stats.mean.(rois{1}).item);
mat = NaN(length(rois), 2, n_subj);
for i = 1:length(rois)
    s_corr = stats.corr.(rois{i});
    s_incorr = stats.incorr.(rois{i});
    
    % everyone has at least one incorrect item, but not always one
    % incorrect item per category
    corr_item = nanmean(s_corr.item_cat - s_corr.within_cat, 2);
    incorr_item = nanmean(s_incorr.item_cat - s_incorr.within_cat, 2);
    
    mat(i,:,:) = [corr_item incorr_item]';
end

disp('Within Region')
disp('ROIs:')
disp(rois)
disp('correct:')
ttest_tab(mat(:,1,:), 0, 'dim', 3);
disp('incorrect:')
ttest_tab(mat(:,2,:), 0, 'dim', 3);
disp('item reactivation SME:')
ttest_tab(mat(:,1,:), mat(:,2,:), 'dim', 3);

disp('Region Differences')
disp('region correct:')
ttest_tab(mat(1,1,:), mat(2,1,:), 'dim', 3);
disp('region incorrect:')
ttest_tab(mat(1,2,:), mat(2,2,:), 'dim', 3);
disp('region x memory interaction:')
ttest_tab(mat(1,1,:)-mat(1,2,:), mat(2,1,:)-mat(2,2,:), 'dim', 3);

clf
y = nanmean(mat, 3);
err = nanstd(mat, [], 3) ./ sqrt(sum(~isnan(mat), 3) - 1);
colors = [.3 .3 .3
          .7 .7 .7];
%colors = [.3255 .5608 .4980   % item reactivation, correct
%          .5569 .7529 .7255]; % category reactivation, incorrect
[hbar, herr] = ebar([], y, err);
for i = 1:length(hbar)
    hbar(i).FaceColor = colors(i,:);
    hbar(i).LineStyle = 'none';
    herr(i).CapSize = 15;
    herr(i).LineWidth = 1;
end
set(gca, 'XTick', 1:length(rois), 'XLim', [0.5 length(rois)+.5], ...
         'XTickLabel', rois)
%ylabel('reactivation')
set_fig_style(gcf, 'preston')
legend off
set(gca, 'YLim', [-.12 .12], 'YTick', -.12:.04:.12)

aspect = [3.2 4];
scale = 1.5;
set(gcf, 'PaperPosition', [0 0 aspect(1)*scale aspect(2)*scale]);

if ~isempty(fig_file)
    print(gcf, '-depsc', fig_file);
end
