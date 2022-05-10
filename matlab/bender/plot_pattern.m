function mat = plot_pattern(fig_file, varargin)
%PLOT_PATTERN   Plot a randomly generated pattern or RDM.
%
%  Used to quickly create sample vectors and matrices for
%  schematics.
%
%  mat = plot_pattern(fig_file, ...)
%
%  OPTIONS
%  shape - [1 x dimensions] array
%      rows and columns of the matrix to be generated. ([5 5])
%
%  mat - [rows x columns] array
%      values for the matrix to plot. If not specified, a random matrix
%      will be generated. ([])
%
%  colormap - char
%      name of a colormap to use. See documentation for colormap. ('gray')
%
%  rdm - boolean
%      if true, will create a random symmetric dissimilarity
%      matrix. (false)

def.shape = [5 5];
def.mat = [];
def.colormap = 'gray';
def.rdm = false;
opt = propval(varargin, def);

if isempty(opt.mat)
    if opt.rdm
        if length(unique(opt.shape)) > 1
            error('For an RDM, shape must be symmetric.')
        end
        n = (prod(opt.shape) - opt.shape(1)) / 2;
        mat = squareform(rand(1, n));
    else
        mat = rand(opt.shape);
    end
else
    mat = opt.mat;
    opt.shape = size(mat);
end

clf
imagesc(mat, [0 1])
colormap(opt.colormap)
axis off image

ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
ax.Position = [0 0 1 1];

fig = gcf;
fig.PaperPosition = [0 0 (opt.shape(2)/opt.shape(1)) 1];
print(gcf, '-dpng', fig_file)
