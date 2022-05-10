function print_item_react_corr(con, rois, fig_file)
%PRINT_ITEM_REACT_CORR   Plot item reactivation correlation stats.
%
%  print_item_react_corr(con, rois, fig_file)

disp('main effect:')
disp(con.self.mean.sigtab)
disp('correct:')
disp(con.self.corr.sigtab)
disp('incorrect:')
disp(con.self.incorr.sigtab)
disp('SME:')
disp(con.self.sme.sigtab)

% pHPC correlation with item reactivation regions
clf
%colors = [1.0000 0.2784 0.5647
%          1.0000 0.5784 0.8647];
colors = [.3 .3 .3
          .7 .7 .7];
cmat = con.self.corr.b(:,1:3,4);
imat = con.self.incorr.b(:,1:3,4);
mat = permute(cat(3, cmat, imat), [2 3 1]);
y = nanmean(mat, 3);
err = nanstd(mat, [], 3) ./ sqrt(sum(~isnan(mat), 3) - 1);
[hbar, herr] = ebar([], y, err);
for i = 1:length(hbar)
    herr(i).CapSize = 20;
    hbar(i).LineStyle = 'none';
    hbar(i).FaceColor = colors(i,:);
end
ylabel('coupling with pHPC')
set(gca, 'XTick', 1:3, 'XLim', [.5 3.5], ...
         'XTickLabel', {'RPHC' 'LPRC' 'RIFG'})
set(gca, 'YLim', [-.1 .8], 'YTick', 0:.2:.8)
%l = legend({'correct' 'incorrect'});
%l.Location = 'NorthWest';
set_fig_style(gcf, 'preston')
legend('off')
if ~isempty(fig_file)
    print(gcf, '-depsc', fig_file)
end

disp('Reactivation connectivity with pHPC:')
disp('ROIs:')
disp(rois(1:3))
disp('correct')
ttest_tab(mat(:,1,:), 0, 'dim', 3);
disp('incorrect')
ttest_tab(mat(:,2,:), 0, 'dim', 3);
