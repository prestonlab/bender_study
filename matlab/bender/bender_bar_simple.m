function h = bender_bar_simple(mat, varargin)
%BENDER_BAR   Make a standard Bender bar plot.
%
%  bender_bar(mat, ...)
%
%  OPTIONS:
%   tests
%   plot_indiv (false)
%   gray (false)
%   oneplot (true)
%   label_type ('none')
%   outline (false)

def.tests = {'AB1' 'AB2' 'AB3' 'AB4' 'AC' 'BC' 'XY' 'AB'};
def.test_labels = {'1' '2' '3' '4' 'AC' 'BC' 'XY' 'AB'};
def.y_label = 'fraction correct';
def.x_label = '';
[def.col, def.lcol] = bender_colors();
def.fit = [];
def.plot_indiv = false;
def.indiv_color = true;
def.plot_err = false;
def.plot_chance = false;
def.gray = false;
def.oneplot = true;
def.label_type = 'below';
def.fig_style = 'preston';
def.outline = false;
def.xbar = [1 2 3 4 6 7 8 9];
def.sep_days = true;
opt = propval(varargin, def, 'strict', false);

if opt.oneplot
  clf
end

if isstruct(mat)
  c = struct2cell(mat);
  mat = cat(1, c{:})';
end

hold on
m = nanmean(mat, 1);
for i = 1:size(mat, 2)
  x = opt.xbar(i);
  
  % get colors for this test
  if opt.gray
    col = [.6 .6 .6];
    lcol = [.8 .8 .8];
  else
    if any(isstrprop(opt.tests{i}, 'digit')) && length(opt.tests{i}) > 2
      test_type = lower(opt.tests{i}(1:2));
    else
      test_type = lower(opt.tests{i});
    end
    col = opt.col.(test_type);
  end
  
  % plot bars
  h.bar(i) = bar(x, m(i), 'FaceColor', col, 'EdgeColor', 'none');
  
  % add individual subject performance
  if opt.plot_indiv
    xi = repmat(x, [size(mat, 1) 1]);
    if opt.indiv_color && ~isempty(opt.lcol) && isfield(opt.lcol, test_type)
      lcol = opt.lcol.(test_type);
      h.sub(i) = plot(xi, mat(:,i), 'ok', 'MarkerFaceColor', lcol, ...
                   'MarkerEdgeColor', 'none', 'MarkerSize', 12);
    else
      h.sub(i) = plot(xi, mat(:,i), 'ok', 'MarkerSize', 14, 'LineWidth', 1);
    end
  end
  
  if opt.plot_err
    % standard error for participants with defined data
    sem = nanstd(mat(:,i)) / sqrt(nnz(~isnan(mat(:,i)))-1);
    h.err(i) = errorbar(x, m(i), sem);
    set(h.err(i), 'Color', 'k', 'LineWidth', 1, 'CapSize', 20)
  end
  
  if opt.plot_chance
    plot(get(gca, 'XLim'), [1/3 1/3], '-k', 'LineWidth', 1);
  end
end

set(get(h.bar(1), 'Baseline'), 'LineWidth', 1)

switch opt.label_type
  case 'below'
    set(gca, 'XTick', opt.xbar, 'XTickLabel', opt.test_labels)
  case 'legend'
    l = legend(h.bar, opt.tests);
    set(gca, 'XTick', [], 'XTickLabel', {});
  case 'none'
    set(gca, 'XTick', [], 'XTickLabel', {});
    opt.x_label = '';
end

set(gca, 'XLim', [0 max(opt.xbar) + 1])
axis square

if strcmp(opt.label_type, 'below')
  set(gca, 'FontSize', 22)
end
% xl = xlabel('AB');
% pos = xl.Position;
% pos(1) = 2.5;
% xl.Position = pos;

if opt.outline
  for i = 1:length(hbar)
    set(h.bar(i), 'LineStyle', '-', 'EdgeColor', 'k', 'LineWidth', 1)
  end
end

if opt.sep_days
  plot([5 5], get(gca, 'YLim'), '--k', 'LineWidth', 1)
end

y_lim = get(gca, 'YLim');
set(gca, 'LineWidth', 1)
plot(get(gca, 'XLim'), [y_lim(1) y_lim(1)], '-k', 'LineWidth', 1)
