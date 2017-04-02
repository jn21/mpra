function lm_table = convert_data_table_to_lm_table(mpra)
%Takes as input the mpra_data table. Converts it to a table for use with
%fitlm. 

unique_up_prefix = unique(mpra{:,'upstream_prefix'});
unique_dn_prefix = unique(mpra{:,'dnstream_prefix'});

unique_up_prefix_str = cellfun(@(s) sprintf('up_%s',s),unique_up_prefix,'Uni',false);
unique_dn_prefix_str = cellfun(@(s) sprintf('dn_%s',s),unique_dn_prefix,'Uni',false);

num_up_id = numel(unique_up_prefix);
num_dn_id = numel(unique_dn_prefix);

num_rows = numel(mpra{:,'upstream_prefix'});
num_vars = num_up_id + num_dn_id + 5;

lm_table = array2table(zeros(num_rows,num_vars));

lm_table.Properties.VariableNames = [unique_up_prefix_str' unique_dn_prefix_str' 'PAS_change' 'addStrongPAS' 'U1_change' 'E_ratio' 'P_ratio'];

for ii = 1:num_up_id
    [~,locb] = ismember(mpra{:,'upstream_prefix'},unique_up_prefix(ii));
    lm_table{logical(locb),sprintf('up_%s',unique_up_prefix{ii})} = 1;
end

for ii = 1:num_dn_id
    [~,locb] = ismember(mpra{:,'dnstream_prefix'},unique_dn_prefix(ii));
    lm_table{logical(locb),sprintf('dn_%s',unique_dn_prefix{ii})} = 1;
end

% lm_table{:,'upstream_is_reverse'} = mpra{:,'upstream_is_reverse'};
% lm_table{:,'dnstream_is_reverse'} = mpra{:,'dnstream_is_reverse'};
lm_table{:,'U1_change'} = mpra{:,'dnstream_num_addU1'} - mpra{:,'dnstream_num_delU1'};
%lm_table{:,'addPAS'} = mpra{:,'dnstream_addPAS'};
lm_table{:,'PAS_change'} = mpra{:,'dnstream_addPAS'} - mpra{:,'dnstream_delPAS'};
lm_table{:,'addStrongPAS'} = mpra{:,'dnstream_addStrongPAS'};
lm_table{:,'E_ratio'} = mpra{:,'E_ratio_avg_rep'};
lm_table{:,'P_ratio'} = mpra{:,'P_ratio_avg_rep'};

end

