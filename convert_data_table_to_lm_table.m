function [lm_table,mpra] = convert_data_table_to_lm_table(mpra,id_only,remove_reverse)
%Takes as input the mpra_data table. Converts it to a table for use with
%fitlm. 
%
%Input
%   mpra: data table
%   id_only: (Boolean) If true, only uses id's as predictors. Also only
%       includes observations that do not have downstream modifications. 

if id_only
    mpra = subset_table(mpra,'dnstream_is_modified',0);
end

if remove_reverse
    mpra = subset_table(mpra,'dnstream_is_reverse',0);
    mpra = subset_table(mpra,'upstream_is_reverse',0);
end

unique_up_prefix = unique(mpra{:,'upstream_prefix'});
unique_dn_prefix = unique(mpra{:,'dnstream_prefix'});

unique_up_prefix_str = cellfun(@(s) sprintf('up_%s',s),unique_up_prefix,'Uni',false);
unique_dn_prefix_str = cellfun(@(s) sprintf('dn_%s',s),unique_dn_prefix,'Uni',false);

num_up_id = numel(unique_up_prefix);
num_dn_id = numel(unique_dn_prefix);

num_rows = numel(mpra{:,'upstream_prefix'});

if id_only
    num_vars = num_up_id + num_dn_id + 2;
else
    num_vars = num_up_id + num_dn_id + 5;
end

lm_table = array2table(zeros(num_rows,num_vars));

if id_only
    lm_table.Properties.VariableNames = [unique_up_prefix_str' unique_dn_prefix_str' 'promoter_activity' 'enhancer_activity'];
else
    lm_table.Properties.VariableNames = [unique_up_prefix_str' unique_dn_prefix_str' 'PAS_change' 'addStrongPAS' 'U1_change' 'E_ratio' 'P_ratio'];
end

%% ID's
for ii = 1:num_up_id
    [~,locb] = ismember(mpra{:,'upstream_prefix'},unique_up_prefix(ii));
    lm_table{logical(locb),sprintf('up_%s',unique_up_prefix{ii})} = 1;
end

for ii = 1:num_dn_id
    [~,locb] = ismember(mpra{:,'dnstream_prefix'},unique_dn_prefix(ii));
    lm_table{logical(locb),sprintf('dn_%s',unique_dn_prefix{ii})} = 1;
end

%% Dnstream signals
if ~id_only
    lm_table{:,'U1_change'} = mpra{:,'dnstream_addU1'} - mpra{:,'dnstream_num_delU1'};
    lm_table{:,'PAS_change'} = mpra{:,'dnstream_addPAS'} - mpra{:,'dnstream_delPAS'};
    lm_table{:,'addStrongPAS'} = mpra{:,'dnstream_addStrongPAS'};
end

%% Response variables
lm_table{:,'promoter_activity'} = mpra{:,'promoter_activity'};
lm_table{:,'enhancer_activity'} = mpra{:,'enhancer_activity'};

end

