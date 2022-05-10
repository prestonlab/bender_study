function data = cat_data(data1, data2, varargin)
%CAT_DATA   Concatenate all fields of two data structures.
%
%   DATA = CAT_DATA(DATA1, DATA2) creates DATA by concatenating
%   each field in DATA1 with the corresponding field in DATA2.  If
%   a field in DATA1 is not in DATA2, an error is raised.
%
%   data = cat_data(data1, data2, ...)
%
%  INPUTS:
%      data1:  a data structure
%
%      data2:  a data structure
%
%  OPTIONS:
%   scalars - a cell array of strings, each string being the name
%             of a field to be made a scalar or not (as indicated
%             by include) in the resulting data structure ({})

def.scalars = {};
def.padval = 0;
opt = propval(varargin, def);

% check the input data structs
data = struct();
if ~(isstruct(data1) || isstruct(data2))
  error('You must supply at least one structure.');
end
if isempty(data2)
  data = data1;
  return
end
if isempty(data1)
  data = data2;
  return
end

% first use the fieldnames on data1, then any new fields on
% data2. This will preserve the order of fields
fnames = [fieldnames(data1); setdiff(fieldnames(data2), fieldnames(data1))];

% fields to scalarize
f_scalar = ismember(fnames, opt.scalars);

for f = 1:length(fnames)
  field = fnames{f};
  
  if isfield(data1, field) && ~isfield(data2, field)
    % this field only exists on the first struct
    data.(field) = data1.(field);
  elseif isfield(data2, field) && ~isfield(data1, field)
    % this field only exists on the second struct
    data.(field) = data2.(field);
  elseif isstruct(data1.(field)) || isstruct(data2.(field))
    % recursively concatenate structure field
    data.(field) = cat_data(data1.(field), data2.(field));
  else
    % if one array is numeric, and the other is a cell array,
    % convert the numeric to cell
    if iscell(data1.(field)) && isnumeric(data2.(field))
      data2.(field) = num2cell(data2.(field));
    elseif iscell(data2.(field)) && isnumeric(data1.(field))
      data1.(field) = num2cell(data1.(field));
    end
    
    if isnumeric(data1.(field))
      % pad with zeros
      data.(field) = padcat(1, opt.padval, data1.(field), data2.(field));
    else
      % use the default padding for this type of data
      data.(field) = padcat(1, [], data1.(field), data2.(field));
    end
  end
  
  if f_scalar(f)
    % make this field a scalar
    data.(field) = make_scalar(data.(field), field);
  end
end
%endfunction


function field = make_scalar(field, fname)
%MAKE_SCALAR  Takes a field from the data structure and attempts to
%             make it a scalar value.  If the values within field
%             are not unique, an error will be raised.

if ~isscalar(field)
  
  if iscell(field)
    % unique nonempty elemnts
    uniq_field = unique(field(~cellfun(@isempty, field)));
    if isempty(uniq_field)
      uniq_field = {[]};
    end
  else
    % unique non-NaN elements
    uniq_field = unique(field(~isnan(field)));
    if isempty(uniq_field)
      uniq_field = [NaN];
    end
  end
  
  %uniq_field = unique(field);
  if ~isscalar(uniq_field)
    errmsg = sprintf(['Field %s does not contain unique values; ' ...
                      'unable to make a scalar.'], fname);
    error(errmsg);
  end

  field = uniq_field;
end
%endfunction