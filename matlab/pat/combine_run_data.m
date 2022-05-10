function data = combine_run_data(hdr, periods)
%COMBINE_RUN_DATA   Combine data from all phases and runs.
%
%  data = combine_run_data(hdr)

if nargin < 2
  periods = hdr.par.periods;
end

data = struct;
for i = 1:length(periods)
  if ~isfield(hdr.data, periods{i})
    data.(periods{i}) = struct;
    continue
  end

  % get all run data structs
  data.(periods{i}) = struct;
  data_files = hdr.data.(periods{i}).mat;
  if isfield(hdr.data.(periods{i}), 'rec')
    rec_files = hdr.data.(periods{i}).rec;
  else
    rec_files = {};
  end
  c = cell(1, length(data_files));
  for j = 1:length(data_files)
    if exist(data_files{j}, 'file')
      run_data = getfield(load(data_files{j}, 'data'), 'data');
      run_data.run = j;
      
      if ~isempty(rec_files)
        score_file = [rec_files{j}(1:end-4) '.csv'];
        if exist(score_file, 'file')
          [correct, resp] = read_cr_score(score_file);
          if length(correct) < length(run_data.trial)
            error('Failed to read in all scores for %s.', ...
                  data_files{j});
          end

          run_data.correct = correct;
          run_data.resp = resp;
        end
      end
      c{j} = run_data;
    end
  end
  
  if all(cellfun(@isempty, c))
    data.(periods{i}) = struct;
    continue
  end
  c = c(~cellfun(@isempty, c));
  
  % get fieldnames for data from all runs
  f = cellfun(@fieldnames, c, 'UniformOutput', false);
  all_f = unique(cat(1, f{:}), 'stable');
  
  % get datatype for every field
  t = cellfun(@(x) structfun(@isnumeric, x), c, 'UniformOutput', false);
  if length(t) > 1
    isnum = padcat(2, NaN, t{:});
  else
    isnum = t{1};
  end
  
  % concatenate the runs
  for j = 1:length(c)
    for k = 1:length(all_f)
      if ~isfield(c{j}, all_f{k})
        % field doesn't exist for this run; add an appropriate placeholder
        if isnum(k,~isnan(isnum(k,:)))
          c{j}.(all_f{k}) = NaN;
        else
          c{j}.(all_f{k}) = {[]};
        end
      end
    end
    data.(periods{i}) = cat_data(data.(periods{i}), c{j});
  end
end

