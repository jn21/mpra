function out_table = subset_table(in_table,field,values)
%Subset an input table by the values in a given field

[~, idx] = ismember(in_table{:,field},values);
out_table = in_table(logical(idx),:);

end

