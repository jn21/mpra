function T = fast_parse_construct_name(construct_names)
%Parse the construct name to extract relevant annotations

cell_var_names = {'oligo_id',...
    'construct',...
    'construct_no_oligo_id',...
    'upstream_full_id',...
    'dnstream_full_id',...
    'upstream_prefix',...
    'dnstream_prefix',...
    'upstream_region_id',...
    'dnstream_region_id',...
    'upstream_sequence_relation_to_tss',...
    'dnstream_sequence_relation_to_tss'};
array_var_names = {'upstream_is_reverse',...
    'dnstream_is_reverse',...
    'dnstream_addPAS',...
    'dnstream_delPAS',...
    'dnstream_addStrongPAS',...
    'dnstream_num_delU1',...
    'dnstream_addU1'};

T1 = cell2table(cell(length(construct_names),length(cell_var_names)));
T1.Properties.VariableNames = cell_var_names;

T2 = array2table(zeros(length(construct_names),length(array_var_names)));
T2.Properties.VariableNames = array_var_names;
    
T = [T1 T2];

%% Part 1 - construct names
T{:,'construct'} = construct_names;

%% Part 2 - id's
temp = cellfun(@(s) strsplit(s,{':','-'}),construct_names,'uni',false);
temp = vertcat(temp{:});
%temp = [temp{:}]

T{:,'oligo_id'} = temp(:,1);
T{:,'construct_no_oligo_id'} = strcat(temp(:,2),'-',temp(:,3));
T{:,'upstream_full_id'} = temp(:,2);
T{:,'dnstream_full_id'} = temp(:,3);

%% Part 3 - regions
up_temp = cellfun(@(s) strsplit(s,'_'), temp(:,2), 'uni', false);
dn_temp = cellfun(@(s) strsplit(s,'_'), temp(:,3), 'uni', false);

T{:,'upstream_region_id'} = cellfun(@(s) s{1}, up_temp, 'uni', false);
T{:,'dnstream_region_id'} = cellfun(@(s) s{1}, dn_temp, 'uni', false);

%% Part 4 - remove pas/u1 from id's, but keep reverse
strings_to_remove = {'_addPAS','_addStrongPAS','_delU1','_addU1','_delPAS'};
dnstream_prefix_temp = temp(:,3);

for ii = 1:length(strings_to_remove)
    dnstream_prefix_temp = strrep(dnstream_prefix_temp,...
        strings_to_remove{ii},...
        '');
end

T{:,'upstream_prefix'} = T{:,'upstream_full_id'};
T{:,'dnstream_prefix'} = dnstream_prefix_temp;

%% Part 5 - is the sequence from an enhancer or a promoter
up_gt_10kb_temp = cellfun(@(s) strfind(s,'gt'),T{:,'upstream_prefix'},'uni',false);
up_gt_10kb_idx = logical(cellfun(@(s) length(s),up_gt_10kb_temp));
T{up_gt_10kb_idx,'upstream_sequence_relation_to_tss'} = {'gt_10000nt'};
T{~up_gt_10kb_idx,'upstream_sequence_relation_to_tss'} = {'wt_100nt'};

dn_gt_10kb_temp = cellfun(@(s) strfind(s,'gt'),T{:,'dnstream_prefix'},'uni',false);
dn_gt_10kb_idx = logical(cellfun(@(s) length(s),dn_gt_10kb_temp));
T{dn_gt_10kb_idx,'dnstream_sequence_relation_to_tss'} = {'gt_10000nt'};
T{~dn_gt_10kb_idx,'dnstream_sequence_relation_to_tss'} = {'wt_100nt'};

%% Part 6 - Reverse?
up_is_reverse_temp = cellfun(@(s) strfind(s,'reverse'),T{:,'upstream_full_id'},'uni',false);
up_is_reverse_idx = logical(cellfun(@(s) length(s), up_is_reverse_temp));
T{up_is_reverse_idx,'upstream_is_reverse'} = 1;

dn_is_reverse_temp = cellfun(@(s) strfind(s,'reverse'),T{:,'dnstream_full_id'},'uni',false);
dn_is_reverse_idx = logical(cellfun(@(s) length(s), dn_is_reverse_temp));
T{dn_is_reverse_idx,'dnstream_is_reverse'} = 1;

%% Part 7 - PAS
dn_has_addpas_temp = cellfun(@(s) strfind(s,'addPAS'),T{:,'dnstream_full_id'},'uni',false);
dn_has_addpas_idx = logical(cellfun(@(s) length(s), dn_has_addpas_temp));
T{dn_has_addpas_idx,'dnstream_addPAS'} = 1;

dn_has_delpas_temp = cellfun(@(s) strfind(s,'delPAS'),T{:,'dnstream_full_id'},'uni',false);
dn_has_delpas_idx = logical(cellfun(@(s) length(s), dn_has_delpas_temp));
T{dn_has_delpas_idx,'dnstream_delPAS'} = 1;

dn_has_strong_pas_temp = cellfun(@(s) strfind(s,'addStrongPAS'),T{:,'dnstream_full_id'},'uni',false);
dn_has_strong_pas_idx = logical(cellfun(@(s) length(s), dn_has_strong_pas_temp));
T{dn_has_strong_pas_idx,'dnstream_addStrongPAS'} = 1;

%% Part 8 - U1
dn_num_add_U1_temp = cellfun(@(s) strfind(s,'addU1'),T{:,'dnstream_full_id'},'uni',false);
dn_num_add_U1 = cellfun(@(s) length(s), dn_num_add_U1_temp);
T{:,'dnstream_addU1'} = dn_num_add_U1;

dn_num_del_U1_temp = cellfun(@(s) strfind(s,'delU1'),T{:,'dnstream_full_id'},'uni',false);
dn_num_del_U1 = cellfun(@(s) length(s), dn_num_del_U1_temp);
T{:,'dnstream_num_delU1'} = dn_num_del_U1;

end

