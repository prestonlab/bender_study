function sorted = sort_period_data(data, sort_fields, index, n_trial)
%SORT_PERIOD_DATA   Sort data from a period by some property of the trials.
%
%  Sort a data structure with a number of fields, including
%  vectors (one element per run) and matrices (runs X trials). Can
%  specify any field to do the sorting, for example the group number
%  of the stimulus presented on a given trial. This is not always
%  unique; for example, there may be multiple presentations of a
%  stimulus. You can specify additional fields to do multiple sorts,
%  in descending priority. Trials will be collapsed over run, so that
%  each field is vectorized. Vector fields will be expanded to
%  matrices so their information will be preserved in this new format;
%  scalar or irregular (width that doesn't match the number of
%  trials) fields will be omitted.
%
%  Fields used for sorting must contain numeric data in [runs x
%  trials] matrices. Every field must have the same number of rows.
%
%  sorted = sort_period_data(data, sort_fields, index, n_trial)
%
%  INPUTS:
%         data:  structure with one or more fields, with any of the
%                following types of arrays:
%                 vector - has shape of [runs X 1].
%                 matrix - has shape of [runs X trials].
%
%  sort_fields:  cell array of strings with the names of fields to
%                sort by. Sorting will prioritize the fields in
%                descending order.
%
%        index:  numeric index with one element per run. The order
%                of the rows in the sorted struct is determined by
%                the order of the values in index.
%
%      n_trial:  (optional) number of trials contained in each
%                standard data matrix. If not specified, this is
%                assumed to be the maximum number of columns
%                over all of the subfields.
%
%  OUTPUTS:
%   sorted:  sorted version of data.
%
%  EXAMPLE:
%  s = struct('trial', [1 2 3; 1 2 3], 'run', [1 1 1; 2 2 2], ...
%             'item', [1 2 3; 4 5 6]);
%  sort_period_data(s, {'run' 'trial'}, [1 1]')

if ischar(sort_fields)
  sort_fields = {sort_fields};
end

f = fieldnames(data);
n_field = length(f);
n_sort = length(sort_fields);
c = cellfun(@(x) size(x), struct2cell(data), 'UniformOutput', false);
fsize = padcat(1, 1, c{:});
max_dim = size(fsize, 2);
maxsize = max(fsize, [], 1);
n_row = fsize(:,1);
n_col = fsize(:,2);
if length(unique(n_row)) > 1
  error('Number of rows must be the same for each field.')
end
n_row = n_row(1);

if nargin < 4 || isempty(n_trial)
  n_trial = maxsize(2);
end
if nargin < 3 || isempty(index)
  index = ones(n_row, 1);
end

% expand vector fields
for i = 1:n_field
  if all(fsize(i,2:end) == 1)
    data.(f{i}) = repmat(data.(f{i}), 1, n_trial);
    n_col(i) = n_trial;
  end
end

uindex = unique(index);
n_index = length(uindex);
sorted = struct;
for i = 1:n_index
  d = trial_subset(index == uindex(i), data);
  s = struct;
  
  % determine sorting for this data subset
  sort_mat = [];
  for j = 1:n_sort
    if isfield(data, sort_fields{j})
      sort_mat = [sort_mat d.(sort_fields{j})(:)];
    end
  end
  
  if isempty(sort_mat)
    error('Found none of the sort fields on the data struct.')
  end
  
  [~, sort_ind] = sortrows(sort_mat);

  for j = 1:n_field
    if n_col(j) == n_trial
      if (max_dim > 2 && fsize(j,3) > 1) && ...
            (max_dim < 4 || all(fsize(j,4:end) == 1))
        % 3D matrix
        matsize = size(d.(f{j}));
        temp = reshape(d.(f{j}), [prod(matsize(1:2)), matsize(3)]);
        x = permute(temp(sort_ind,:), [3 1 2]);
      else
        % 2D matrix
        x = d.(f{j})(sort_ind);
        if iscolumn(x)
          x = x';
        end
      end
      s.(f{j}) = x;
    end
  end
  sorted = cat_data(sorted, s, 'padval', NaN);
end


