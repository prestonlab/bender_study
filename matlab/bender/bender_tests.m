function test = bender_tests(data)
%BENDER_TESTS   Put information about Bender tests in standard format.
%
%  test = bender_tests(data)

% sort data so there is one subject per row
s = sort_data(data, {'run' 'trial'});

% split the AB test by repetition (run)
test = struct;

if isfield(s, 'ab_feedback')
  n_run = length(unique(data.ab_feedback.run));
  for i = 1:n_run
    f = ['ab' num2str(i)];
    test.(f) = trial_subset(s.ab_feedback.run(1,:) == i, s.ab_feedback, 2);
  end
end

if isfield(s, 'ac_test')
  test.ac = s.ac_test;
end

if isfield(s, 'bcxy_test')
  % BC and XY are mixed within each run; separate them out
  test.bc = cond_subset(s.bcxy_test.cond ~= 5, s.bcxy_test);
  test.xy = cond_subset(s.bcxy_test.cond == 5, s.bcxy_test);
end
  
if isfield(s, 'ab_test')
  test.ab = s.ab_test;
end

% sort by group number
test = sort_data(test, {'group' 'cue_group'});
