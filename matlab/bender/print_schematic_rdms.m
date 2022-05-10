function print_schematic_rdms(model, a_items, c_items, fig_dir)
%PRINT_SCHEMATIC_RDMS   Print plots of schematic RDMs.
%
%  print_schematic_rdms(model, a_items, c_items, fig_dir)

% construct the corresponding model vectors
a_vecs = [];
ac_vecs = [];
for i = 1:length(a_items)
    a_ind = find(strcmp(model.items, a_items{i}));
    c_ind = find(strcmp(model.items, c_items{i}));
    
    a_vecs = [a_vecs; model.vectors(a_ind,:)];
    ac_vecs = [ac_vecs; model.vectors(a_ind,:) + model.vectors(c_ind,:)];
end

% create rdms
a_rdm = squareform(pdist(a_vecs, 'correlation'));
ac_rdm = squareform(pdist(ac_vecs, 'correlation'));

% print plots
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end

% A
plot_rdm(a_rdm);
set(gca, 'CLim', [.1 .25])
colorbar off
set(gca, 'Position', [0 0 1 1])
set(gcf, 'PaperPosition', [0 0 1 1])
print(gcf, '-dpng', fullfile(fig_dir, 'model_a.png'));

% neural
shift = [ 0   0  .025  .035
          0   0  .05  .05
         .025  .05   0   0
         .035  .05   0   0];
plot_rdm(a_rdm + shift);
set(gca, 'CLim', [.1 .3])
colorbar off
set(gca, 'Position', [0 0 1 1])
set(gcf, 'PaperPosition', [0 0 1 1])
print(gcf, '-dpng', fullfile(fig_dir, 'model_a_noise.png'));

% AC
plot_rdm(ac_rdm);
set(gca, 'CLim', [.05 .2])
colorbar off
set(gca, 'Position', [0 0 1 1])
set(gcf, 'PaperPosition', [0 0 1 1])
print(gcf, '-dpng', fullfile(fig_dir, 'model_ac.png'));

% colorbar
imagesc(linspace(0, 1, 1000)');
axis xy
set(gca, 'Position', [0 0 1 1])
set(gcf, 'PaperPosition', [0 0 1 30])
set(gca, 'XTick', [], 'YTick', [], 'Visible', 'off')
print(gcf, '-dpng', fullfile(fig_dir, 'colorbar.png'));
