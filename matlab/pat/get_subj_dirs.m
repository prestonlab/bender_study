function subj = get_subj_dirs(dataroot, subjstr)
%GET_SUBJ_DIRS   Get data directories for a set of subjects.
%
%  subj = get_subj_dirs(dataroot, subjstr)

d = dir(fullfile(dataroot, subjstr));
subjects = {d.name};
subj = [];
for i = 1:length(subjects)
  this_subj.id = subjects{i};
  subj_dir = fullfile(dataroot, subjects{i});
  this_subj.dir = subj_dir;
  this_subj.sess = struct('number', 1, 'dir', subj_dir);
  
  subj = addobj(subj, this_subj);
end

