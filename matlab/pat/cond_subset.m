function y = cond_subset(mask, x)
%COND_SUBSET   Get a subset of each column of all matrices on a struct.
%
%  y = cond_subset(mask, x)

if length(unique(sum(mask, 2))) > 1
  error('Mask must include the same number of columns for each row.')
end

y = struct;
f = fieldnames(x);
for i = 1:length(f)
  fullmat = x.(f{i});
  if size(fullmat, 2) ~= size(mask, 2)
    continue
  end
  
  if isnumeric(fullmat)
    mat = [];
  elseif iscell(fullmat)
    mat = {};
  end
  
  for j = 1:size(fullmat, 1)
    mat = [mat; fullmat(j,mask(j,:))];
  end
  y.(f{i}) = mat;
end

