function print_behav(perf, rtc, fig_dir)
%PRINT_BEHAV   Plot test accuracy and reaction time.
%
%  print_behav(perf, rtc, fig_dir)

% accuracy
h = bender_bar_simple(perf, 'plot_err', true, 'plot_chance', true);
ylabel('fraction correct')
if ~isempty(fig_dir)
    print(gcf, '-depsc', fullfile(fig_dir, 'test_perf.eps'));
end

% rt
h = bender_bar_simple(rtc, 'plot_err', true);
ylabel('reaction time (ms)')
if ~isempty(fig_dir)
    print(gcf, '-depsc', fullfile(fig_dir, 'test_rt_corr.eps'));
end

% per-subject threshold for above-chance performance
bino_thresh = binoinv(0.95, 120, 1/3) / 120;

% stats
disp('AC accuracy')
ttest_tab(perf.ac, 1/3);
disp('AC binomial test')
fprintf('threshold: %.2f\n', bino_thresh);
fprintf('n above:   %d/%d\n', nnz(perf.ac >= bino_thresh), length(perf.ac));
fprintf('\n');
disp('AC RT (correct trials)')
ttest_tab(rtc.ac);
disp('BC-XY accuracy')
ttest_tab(perf.bc, perf.xy);
disp('BC-XY RT (correct trials)')
ttest_tab(rtc.bc, rtc.xy);
disp('AB final accuracy')
ttest_tab(perf.ab);
