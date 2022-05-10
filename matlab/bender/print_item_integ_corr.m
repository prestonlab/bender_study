function print_item_integ_corr(b, item_rois, integ_rois, fig_file)
%PRINT_ITEM_INTEG_CORR   Print reactivation-integration correlation stats.
%
%  print_item_integ_corr(b, item_rois, integ_rois, fig_file)

disp('item ROIs:')
disp(item_rois)
disp('integration ROIs:')
disp(integ_rois)
disp('Reactivation connectivity for correct trials:')
ttest_tab(b(:,1,:,:), 0, 'dim', 1);
disp('Reactivation connectivity for incorrect trials:')
ttest_tab(b(:,2,:,:), 0, 'dim', 1);
disp('Subsequent memory: item ROI x integration ROI')
ttest_tab(b(:,1,:,:), b(:,2,:,:), 'dim', 1);

% pHPC correlation with item reactivation regions
close all
%colors = [1.0000 0.2784 0.5647
%          1.0000 0.5784 0.8647];
colors = [.3 .3 .3
          .7 .7 .7];
cmat = squeeze(b(:,1,end,:));
imat = squeeze(b(:,2,end,:));

mat = permute(cat(3, cmat, imat), [2 3 1]);
y = nanmean(mat, 3);
err = nanstd(mat, [], 3) ./ sqrt(sum(~isnan(mat), 3) - 1);
[hbar, herr] = ebar([], y, err);
for i = 1:length(hbar)
    herr(i).CapSize = 20;
    hbar(i).LineStyle = 'none';
    hbar(i).FaceColor = colors(i,:);
end
ylabel('r_{self} - r_{AC} coupling')
set(gca, 'XTick', 1:2, 'XLim', [.5 2.5], ...
         'XTickLabel', {'aHPC' 'RPRC'})
set(gca, 'YLim', [-.1 .25], 'YTick', -.1:.05:.25)
%l = legend({'correct' 'incorrect'});
%l.Location = 'NorthEast';
set_fig_style(gcf, 'preston')
legend('off')
aspect = [4 4];
scale = 1.5;
set(gcf, 'PaperPosition', [0 0 aspect(1)*scale aspect(2)*scale]);

if ~isempty(fig_file)
    print(gcf, '-depsc', fig_file)
end
