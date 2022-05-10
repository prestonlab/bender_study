function res = testfun(func, test, varargin)
%TESTFUN   Apply a function to a struct with information about tests.
%
%  Given a standard test information structure (one field per test;
%  each field contains a structure; this structure contains fields
%  with [subjects X trials] arrays), run some function on each
%  test. Given a structure with fields that contain information
%  about one subject, the function must return a numeric
%  scalar. Results will be organized by test, with each containing
%  a [1 X subjects] vector.
%
%  res = testfun(func, test, ...)

def.index_field = 'subj_number';
def.UniformOutput = true;
opt = propval(varargin, def);

% get a list of all subject codes (or some other user-specified index)
% that appear in any of the tests
test_index = structfun(@(s) unique(s.(opt.index_field))', test, ...
                       'UniformOutput', false);
c = struct2cell(test_index);
index = nanunique([c{:}]);
n_index = length(index);

f = fieldnames(test);
for i = 1:length(f)
  if opt.UniformOutput
    x = NaN(1, n_index);
  else
    x = cell(1, n_index);
  end
  f_index = test.(f{i}).(opt.index_field);
  for j = 1:n_index
    % get the trials of interest (e.g. trials for one subject)
    t = trial_subset(f_index(:,1) == index(j), test.(f{i}));
    if size_data(t, 1) == 0
      continue
    end

    % run the function
    if opt.UniformOutput
      x(j) = func(t);
    else
      x{j} = func(t);
    end
  end
  res.(f{i}) = x;
end
