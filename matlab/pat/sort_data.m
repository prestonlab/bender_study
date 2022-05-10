function sorted = sort_data(data, sort_fields, varargin)
%SORT_DATA   Sort by fields in a data structure.
%
%  Arranges data matrices by an index (for example, subject), which
%  appears on the rows of each data matrix, and sorts by other
%  attributes of the data within each row (for example, run and
%  trial).
%
%  sorted = sort_data(data, sort_fields, ...)
%
%  INPUTS:
%         data:  struct with a subfield for each period in an
%                experiment. Each period subfield must contain a set of
%                matrices where different rows indicate different
%                runs. If a period subfield has a column for each trial,
%                or has only one column, it will be included in the
%                sorted output; otherwise, it will be omitted. See
%                sort_period_data for details.
%
%  sort_fields:  names of fields to use to sort the data from each
%                period.
%
%  OPTIONS:
%   index_fields - list of fields to use for defining the rows. If
%                  empty, all data in each period will be vectorized.
%                  ({'subj_number' 'subno' 'subNum'})

def.index_fields = {'subj_number' 'subno' 'subNum'};
opt = propval(varargin, def);

if ischar(opt.index_fields)
  opt.index_fields = {opt.index_fields};
end

sorted = struct;
f = fieldnames(data);
for i = 1:length(f)
  if isempty(fieldnames(data.(f{i})))
    % no participant has run this period yet
    continue
  end
  
  % determine how rows will be defined
  if ~isempty(opt.index_fields)
    % define the index based on a field or the conjunction of
    % multiple fields
    c = {};
    for j = 1:length(opt.index_fields)
      if isfield(data.(f{i}), opt.index_fields{j})
        x = data.(f{i}).(opt.index_fields{j});
        c = [c {x(:,1)}];
      end
    end
    index = make_index(c{:});
  else
    % put every data field into a row vector
    index = [];
  end

  % sort each subfield
  sorted.(f{i}) = sort_period_data(data.(f{i}), sort_fields, index);
end

