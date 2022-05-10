
%% path setup

addpath('~/analysis/aperture')
init_aperture()
addpath('~/analysis/bender/matlab')
addpath('~/analysis/wiki2vec/ana')
addpath('~/analysis/pat/ana')
addpath('~/analysis/aperture/plotting/figures')

%% colormaps

colors = [228,26,28
          55,126,184
          77,175,74
          152,78,163
          255,127,0]/256;
offset = 0.15;
cmapdir = '~/.fsleyes';
filebase = 'colormap5-Set1-%d.cmap';
for i = 1:size(colors, 1)
    c = [colors(i,:) - offset
         colors(i,:) + offset];
    c(c < 0) = 0;
    c(c > 1) = 1;
    filepath = fullfile(cmapdir, sprintf(filebase, i));
    save(filepath, 'c', '-ascii');
end

%% schematic neural patterns

fig_dir = '~/work/bender/figs3/schematic';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end
pat_names = {'a1' 'a2' 'a3' 'enc1' 'enc2'};
for i = 1:length(pat_names)
    fig_name = sprintf('pattern_%s.png', pat_names{i});
    fig_file = fullfile(fig_dir, fig_name);
    plot_pattern(fig_file, 'shape', [5 5]);
end
fig_file = fullfile(fig_dir, 'pattern_circle.png');
plot_pattern(fig_file, 'shape', [9 9]);

%% schematic model RDMs

% load pool
load ~/work/bender/batch/bender_pool.mat
pool_labels = {pool.name};

% get sample items with corresponding similarity in wiki2vec
model = load('~/work/bender/batch/models3/allstim/mat_wiki_w2v.mat');
a_items = {'Julia Roberts' 'Barack Obama' 'Yosemite' 'Eiffel Tower'};
c_items = {'scissors' 'cup' 'butterfly' 'picnic basket'};
fig_dir = '~/work/bender/figs3/schematic';
close all
print_schematic_rdms(model, a_items, c_items, fig_dir);

%% wiki2vec model MDS plots

cat_names = {'face' 'scene' 'object'};
cat_ind = {1:60 61:120 121:4:360};
for i = 1:length(cat_names)
    clf
    if strcmp(cat_names{i}, 'object')
        plot_mds(model.rdm, {pool.crop}, 'ind', cat_ind{i});
    else
        plot_mds(model.rdm, {pool.crop}, 'ind', cat_ind{i}, ...
                 'alpha', {pool.alpha});
    end
    fig_name = sprintf('model_mds_%s.png', cat_names{i});
    fig_file = fullfile(fig_dir, fig_name);
    print(gcf, '-dpng', '-r300', fig_file);
end

%% subject item model RDMs and stimuli

model_dir = '~/work/bender/batch/models3/trial';
fig_dir = '~/work/bender/figs3/schematic';
models = {'a' 'bx' 'cy'};
n_row = 12;
n_col = 10;
close all
for i = 1:length(models)
    model = load(fullfile(model_dir, sprintf('bender_02_wiki_w2v_%s.mat', ...
                                             models{i})));
    plot_rdm(model.rdm(1:120,1:120));
    colorbar off
    set(gca, 'Position', [0 0 1 1])
    set(gcf, 'PaperPosition', [0 0 1 1])
    fig_file = fullfile(fig_dir, sprintf('model_%s_bender_02.png', models{i}));
    print(gcf, '-dpng', fig_file);
    
    clf
    n = 0;
    mat = [];
    for j = 1:n_row
        row = [];
        for k = 1:n_col
            n = n + 1;
            item = pool(strcmp(model.items{n}, pool_labels));
            row = cat(2, row, item.picture);
        end
        mat = cat(1, mat, row);
    end
    imshow(mat);
    set(gca, 'Position', [0 0 1 1]);
    set(gcf, 'PaperPosition', [0 0 10 12]);
    print(gcf, '-djpeg90', ...
          fullfile(fig_dir, sprintf('stim_bender_02_%s.jpg', models{i})))
end

%% behavioral performance plots

% load and sort data
close all
data = bender_load_data();
s = sort_data(data, {'group' 'cue_group'});
subjects = s.bcxy_study.subj_number(:,1);

% calculate averages by test
test = bender_tests(data);
perf = testfun(@(x) mean(x.correct, 2)', test);
rtc = bender_test_rt(test, 1, 'mean');

fig_dir = '~/work/bender/figs3/plots';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end

% plots and stats
print_behav(perf, rtc, fig_dir);

%% load item reactivation and trial RDMs

% spec for all searchlight ROIs
rois = {'rphc' 'lprc' 'rifg' ... % item reactivation
        'phpc' 'ampfc' 'rifg' ... % item suppression
        'mmpfc' ... % semantic reactivation
        'rhpc' 'rprc'}; % semantic integration
roi_types = {'cat_react_item2' 'cat_react_item2' 'cat_react_item2' ...
             'item_suppress_gvt' 'cat_react_item_sme2' 'cat_react_item_sme2' ...
             'study_wiki_w2v_fix_cont_a_bc_sme' ...
             'study_wiki_w2v_fix_cont_ac_bx_sme' 'study_wiki_w2v_fix_cont_ac_bx_sme'};
beta_name = 'study_stim2';
roi_names = {'RPHC' 'LPRC' 'RIFG' ...
             'pHPC' 'amPFC' 'aIFG' ...
             'mmPFC' ...
             'RHPC' 'RPRC'};

% reactivation RDMs and stats
react_roi = bender_roi_rdms('react', rois, roi_types, roi_names, ...
                            beta_name, subjects, '_dil1c');
[react_trial, react_stats] = bender_roi_react_stats(react_roi, ...
                                                  s.ac_test.correct);

% study-phase RDMs
rdm_roi = bender_roi_rdms('rdm', rois, roi_types, roi_names, ...
                          beta_name, subjects, '_dil1c');

%% overall reactivation/suppression

% reactivation by region
close all
rois = {'RPHC' 'LPRC' 'RIFG' 'pHPC'};
print_item_cat_react(react_stats, rois, ...
                     '~/work/bender/figs3/plots/react_item_cat.eps');
print_item_cat_react_sme(react_stats, rois, ...
                         '~/work/bender/figs3/plots/react_item_cat_sme.eps');

% connectivity between regions
con = bender_react_corr_stats(react_trial, s.ac_test.correct, {'self'}, rois);
print_item_react_corr(con, rois, ...
                      '~/work/bender/figs3/plots/phpc_react_coupling.eps');

%% reactivation x subsequent inference

% reactivation by region and inference accuracy
close all
rois = {'aIFG' 'amPFC'};
print_item_cat_react_sme(react_stats, rois, ...
                         '~/work/bender/figs3/plots/react_sme_cat.eps');
print_item_react_sme(react_stats, rois, ...
                     '~/work/bender/figs3/plots/react_sme.eps');

% connectivity between regions
con = bender_react_corr_stats(react_trial, s.ac_test.correct, {'self'}, rois);
disp(con.self.mean.sigtab);

%% does item suppression in aMPFC predict reduced integration in
%% MTL?

% load semantic model RDMs
subjects = s.bcxy_study.subj_number(:,1);
model_names = {'a' 'bx' 'cy' 'ab' 'bcxy' 'ac' 'abc'};
fnames = {'a' 'b' 'c' 'ab' 'bc' 'ac' 'abc'};
model = bender_trial_models(model_names, fnames, subjects);

% calculate correlation between trial reactivation and trial integration
item_rois = {'aIFG' 'amPFC'};
integ_rois = {'RHPC' 'RPRC'};
b = bender_item_integ_corr(react_trial, rdm_roi, model.ac, ...
                           item_rois, integ_rois, s.ac_test.correct);

print_item_integ_corr(b, item_rois, integ_rois, ...
                      '~/work/bender/figs3/plots/react_integ_coupling.eps');
