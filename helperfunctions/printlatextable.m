function printlatextable(T)
%PRINTLATEXTABLE Summary of this function goes here
%   Detailed explanation goes here
% Convert MATLAB table to cell array
table_data = table2cell(T);
[num_rows, ~] = size(table_data);

% Write LaTeX code to a file
fid = fopen('table.tex', 'w');
for row = 1:num_rows
    formatted_row = cellfun(@(x) sprintf('%5.3g', x), table_data(row,:), 'UniformOutput', false);
    fprintf(fid, '%s \\\\\n', strjoin(formatted_row, ' & '));
end
fclose(fid);
end

