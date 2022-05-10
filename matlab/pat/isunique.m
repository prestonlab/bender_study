function tf = isunique(x)
%ISUNIQUE   Test whether an array contains only unique elements.
%
%  tf = isunique(x)

if length(unique(x)) < length(x)
  tf = false;
else
  tf = true;
end
