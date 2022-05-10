function print_item_cat_react(stats, rois, fig_file)
%PRINT_ITEM_CAT_REACT   Plot item and category reactivation by ROI.
%
%  print_item_cat_react(stats, rois, fig_file)

n_subj = length(stats.mean.(rois{1}).item);
mat = NaN(length(rois), 2, n_subj);
for i = 1:length(rois)
    m_item = nanmean(stats.mean.(rois{i}).item_cat - ...
                     stats.mean.(rois{i}).within_cat, 2);
    m_cat = nanmean(stats.mean.(rois{i}).within_cat, 2) - ...
                    stats.mean.(rois{i}).between;
    mat(i,:,:) = [m_item m_cat]';
end

disp('ROIs')
disp(rois)
disp('item reactivation')
ttest_tab(mat(:,1,:), 0, 'dim', 3);
disp('category reactivation')
ttest_tab(mat(:,2,:), 0, 'dim', 3);

clf
y = nanmean(mat, 3);
err = nanstd(mat, [], 3) ./ sqrt(sum(~isnan(mat), 3) - 1);
colors = [.3255 .5608 .4980; % item reactivation
          .3843 .2863 .4784]; % category reactivation          
[hbar, herr] = ebar([], y, err);
for i = 1:length(hbar)
    hbar(i).FaceColor = colors(i,:);
    hbar(i).LineStyle = 'none';
    herr(i).CapSize = 10;
    herr(i).LineWidth = 1;
end
set(gca, 'XTick', 1:length(rois), 'XLim', [0.5 length(rois)+.5], ...
         'XTickLabel', rois)
%l = legend({'item' 'category'});
%l.Location = 'NorthEast';
%ylabel('reactivation')
set_fig_style(gcf, 'preston')
legend off
set(gca, 'YLim', [-.03 .03], 'YTick', -.03:.01:.03)

if ~isempty(fig_file)
    print(gcf, '-depsc', fig_file);
end
