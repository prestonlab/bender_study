function b = robust_reg_rows(x, y, include)
%ROBUST_REG_ROWS   Calculate slope between rows of two variables.
%
%  b = robust_reg_rows(x, y, include)

if nargin < 3
    include = true(size(x));
end

n_row = size(x, 1);
b = NaN(n_row, 2);
for i = 1:n_row
    if nnz(include(i,:)) <= 2
        continue
    end
    
    xi = x(i,include(i,:));
    yi = y(i,include(i,:));

    try
        b(i,:) = robustfit(xi, yi);
    catch
        fprintf('Problem fitting subject index %d.\n', i);
    end
end
