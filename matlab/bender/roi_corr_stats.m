function stats = roi_corr_stats(b, rois, ttype)
%ROI_CORR_STATS   Calculate statistics for correlations between ROIs.
%
%  stats = roi_corr_stats(b, rois)
%
%  INPUTS
%  b - [subjects x roi x roi] numeric array
%      Slope of a relationship between each pair of ROIs, for
%      example calculated using robust regression.
%
%  rois - cell array of strings
%      Name for each ROI.
%
%  OUTPUTS
%  stats - struct
%      Results of t-tests across subjects, and table of significant
%      ROI pairs (uncorrected).

if nargin < 3
    ttype = 'ttest';
end

if isrow(rois)
    rois = rois';
end

stats.b = b;

% test for ROI pairs with a slope different from zero, across subjects
switch ttype
  case 'ttest'
    [hy, pv, ci, st] = ttest(b);
    pv = squeeze(pv);
  case 'signrank'
    pv = NaN(length(rois), length(rois));
    for i = 1:length(rois)
        for j = 1:length(rois)
            try
                pv(i,j) = signrank(b(:,i,j), 0);
            catch
            end
        end
    end
end

mb = squeeze(nanmean(b));
nb = squeeze(sum(~isnan(b)));
semb = squeeze(nanstd(b) / sqrt(nnz(~isnan(b))-1));

[i, j] = find(pv < .05);
ind = sub2ind(size(pv), i, j);
roi1 = rois(i);
roi2 = rois(j);
pvalue = pv(ind);
slope = mb(ind);
sem = semb(ind);
tstat = st.tstat(ind);
df = nb(ind)-1;

stats.pv = pv;
if strcmp(ttype, 'ttest')
    stats.t = squeeze(st.tstat);
    stats.df = squeeze(st.df);
end
stats.sigtab = table(roi1, roi2, slope, sem, tstat, df, pvalue);
