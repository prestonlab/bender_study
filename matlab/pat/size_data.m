function x = size_data(data, dim)
%SIZE_DATA   Size of a data structure with multiple fields.
%
%  Returns the maximum size of each requested dimension over
%  all fields in data.
%
%  x = size_data(data, dim)

n_dim = max(structfun(@ndims, data));
if nargin > 1
  x = ones(1, max(dim));
else
  x = ones(1, n_dim);
end

for i = 1:n_dim
  x(i) = max(structfun(@(x) size(x, i), data));
end

if nargin > 1
  x = x(dim);
end
