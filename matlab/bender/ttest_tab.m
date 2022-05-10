function ttest_tab(varargin)
%TTEST_TAB   Run a t-test and display the results in a table.
%
%  ttest_tab(...)
%
%  All arguments are passed to ttest. Outputs are run through
%  squeeze to remove any singleton dimensions.

[hy, pv, ci, st] = ttest(varargin{:});

if length(varargin) > 1 && isnumeric(varargin{2}) && ~isscalar(varargin{2})
    % comparing two paired samples
    mat = varargin{1} - varargin{2};
    mat1 = varargin{1};
    mat2 = varargin{2};
else
    % one population vs. 0 or some other baseline
    mat = varargin{1};
end

% parse the "dim" input if given
dim_ind = cellfun(@(x) ischar(x) && strcmp(x, 'dim'), varargin);
if any(dim_ind)
    dim = varargin{find(dim_ind)+1};
else
    dim = find(size(mat) > 1, 1);
end

% summary stats
m = squeeze(nanmean(mat, dim));
n = squeeze(sum(~isnan(mat), dim));
sem = squeeze(nanstd(mat, [], dim)) ./ sqrt(n-1);

% hypothesis testing stats
tstat = squeeze(st.tstat);
df = squeeze(st.df);
pv = squeeze(pv);

if exist('mat1', 'var')
    m1 = squeeze(nanmean(mat1, dim));
    n1 = squeeze(sum(~isnan(mat1), dim));
    sem1 = squeeze(nanstd(mat1, [], dim)) ./ sqrt(n1-1);
    
    m2 = squeeze(nanmean(mat2, dim));
    n2 = squeeze(sum(~isnan(mat2), dim));
    sem2 = squeeze(nanstd(mat2, [], dim)) ./ sqrt(n2-1);

    d = m;
    dsem = sem;
    t = table(m1, sem1, m2, sem2, d, dsem, tstat, df, pv);
else
    t = table(m, sem, tstat, df, pv);
end
disp(t);
