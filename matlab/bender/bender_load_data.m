function data = bender_load_data(dataroot, include)
%BENDER_LOAD_DATA   Load behavioral data for Bender.
%
%  data = bender_load_data(dataroot, include)

if nargin < 1 || isempty(dataroot)
  system = getenv('TACC_SYSTEM');
  switch system
    case 'ls5'
      dataroot = '/work/03206/mortonne/lonestar/bender/behav';
    otherwise
      dataroot = '/Users/morton/experiments/nemesis/exp/bender/data';
  end
end

subj = get_subj_dirs(dataroot, 'subj1*');

if nargin < 2
  % same as $SUBJNOS defined in .profile
  includeno = {2 4 5 6 7 8 10 11 12 13 14 15 16 17 18 20 21 22 23 ...
               24 25 26 28 29 31 32 34 35 37 38};
  include = cellfun(@(x) sprintf('subj1%02d', x), includeno, ...
                    'UniformOutput', false);
end

% remove original versions of subjects that have been manually modified
c = cellfun(@(x) strfind(x, '_orig'), {subj.id}, 'UniformOutput', false);
inc = cellfun(@isempty, c);
subj = subj(inc);

% load and concatenate data for all subjects
data = struct;
for i = 1:length(subj)
  if ~isempty(include) && ~ismember(subj(i).id, include)
    continue
  end

  hdr_file = fullfile(dataroot, subj(i).id, 'header.mat');
  load(hdr_file)

  % fix paths
  ind = strfind(hdr.file, subj(i).id);
  dataroot_orig = hdr.file(1:ind-2);
  hdr = struct_strrep(hdr, dataroot_orig, dataroot);
  if strcmp(subj(i).id, 'subj108')
    % localizer was completely redone due to display issues on day
    % 1; completely replace with the second localizer (different
    % order, etc.)
    hdr.data.c_loc = struct_strrep(hdr.data.c_loc, 'subj108', 'subj808');
    subj_data = combine_run_data(hdr);
    subj_data.c_loc.subj_number(:) = 108;
  else
    subj_data = combine_run_data(hdr);
  end
  
  data = cat_data(data, subj_data);
end

% fix for when a field is missing for the first subject
f = fieldnames(data);
for i = 1:length(f)
  c = struct2cell(data.(f{i}));
  n = cellfun(@(x) size(x, 1), c);
  if ~isunique(n)
    f2 = fieldnames(data.(f{i}));
    for j = 1:length(f2)
      if n(j) < max(n)
        x = data.(f{i}).(f2{j});
        data.(f{i}).(f2{j}) = padcat(1, NaN, NaN, x);
      end
    end
  end
end
