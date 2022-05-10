function [correct, resp] = read_cr_score(csv_file);

fid = fopen(csv_file, 'r');
c = textscan(fid, '%d%d%d%d%.3f%.3f%.0f%s%s%d%s', 'Delimiter', ',', ...
             'HeaderLines', 1);

correct = c{end-1}';
resp = c{end}';
if length(resp) < length(correct)
  resp = [resp repmat({''}, 1, length(correct) - length(resp))];
end

fclose(fid);

