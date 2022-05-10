function res = bender_test_rt(test, correct, stat)
%BENDER_TEST_RT   Reaction time for each bender test type.
%
%  res = bender_test_rt(data, correct, stat)
%
%  test - bender test structure.
%
%  correct - 1 or 0 to indicate correct trials or incorrect trials.
%
%  stat - 'all' to return all RTs or 'mean' for just the mean.

if nargin < 3
  stat = 'all';
end
  
get_rt = @(s, correct) s.rt(s.correct == correct);
switch stat
  case 'all'
    res = testfun(@(x) get_rt(x, correct), test, 'UniformOutput', false);
  case 'mean'
    res = testfun(@(x) nanmean(get_rt(x, correct)), test);
  otherwise
    error('Unknown stat: %s', stat);
end
