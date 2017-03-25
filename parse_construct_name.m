function T = parse_construct_name(construct_names)
%Parse the construct name to extract relevant annotations

cell_var_names = {'construct',...
    'oligo_id',...
    'upstream_prefix',...
    'dnstream_prefix',...
    'upstream_id',...
    'dnstream_id'};
array_var_names = {'upstream_sequence_gt_10kb_of_tss',...
    'upstream_sequence_wt_100nt_of_tss',...
    'dnstream_sequence_gt_10kb_of_tss',...
    'dnstream_sequence_wt_100nt_of_tss',...
    'upstream_is_reverse',...
    'dnstream_is_reverse',...
    'dnstream_addPAS',...
    'dnstream_addStrongPAS',...
    'dnstream_num_delU1',...
    'dnstream_num_addU1'};
T1 = cell2table(cell(length(construct_names),length(cell_var_names)));
T1.Properties.VariableNames = cell_var_names;

T2 = array2table(zeros(length(construct_names),length(array_var_names)));
T2.Properties.VariableNames = array_var_names;
    
T = [T1 T2];

for ii = 1:length(construct_names)
    ii
    T{ii,'construct'} = construct_names(ii);
    
    temp = strsplit(construct_names{ii},{':','-'});
    T{ii,'oligo_id'} = temp(1);
    T{ii,'upstream_id'} = temp(2);
    T{ii,'dnstream_id'} = temp(3);
    
    %Up and dn stream prefixes
    up_temp = strsplit(temp{2},'_');
    T{ii,'upstream_prefix'} = up_temp(1);
    dn_temp = strsplit(temp{3},'_');
    T{ii,'dnstream_prefix'} = dn_temp(1);
     
    %is the upstream sequence derived from a promoter or enhancer
    if ~isempty(regexp(temp{2},'within_100nt','once'))
        T{ii,'upstream_sequence_gt_10kb_of_tss'} = false;
        T{ii,'upstream_sequence_wt_100nt_of_tss'} = true;
    else
        T{ii,'upstream_sequence_gt_10kb_of_tss'} = true;
        T{ii,'upstream_sequence_wt_100nt_of_tss'} = false;
    end
    
    %is the dnstream sequence derived from a promoter or enhancer
    if ~isempty(regexp(temp{3},'within_100nt','once'))
        T{ii,'dnstream_sequence_gt_10kb_of_tss'} = false;
        T{ii,'dnstream_sequence_wt_100nt_of_tss'} = true;
    else 
        T{ii,'dnstream_sequence_gt_10kb_of_tss'} = true;
        T{ii,'dnstream_sequence_wt_100nt_of_tss'} = false;
    end
    
    %reverse?
    if ~isempty(regexp(temp{2},'reverse','once'))
        T{ii,'upstream_is_reverse'} = true;
    else
        T{ii,'upstream_is_reverse'} = false;
    end
    if ~isempty(regexp(temp{3},'reverse','once'))
        T{ii,'dnstream_is_reverse'} = true;
    else
        T{ii,'dnstream_is_reverse'} = false;
    end
    
    %PAS
    if ~isempty(regexp(temp{3},'addPAS','once'))
        T{ii,'dnstream_addPAS'} = true;
    else
        T{ii,'dnstream_addPAS'} = false;
    end
    
    %Strong PAS
    if ~isempty(regexp(temp{3},'addStrongPAS','once'))
        T{ii,'dnstream_addStrongPAS'} = true;
    else
        T{ii,'dnstream_addStrongPAS'} = false;
    end
    
    %U1
    num_delU1 = numel(strfind(temp{3},'delU1'));
    num_addU1 = numel(strfind(temp{3},'addU1'));
    T{ii,'dnstream_num_delU1'} = num_delU1;
    T{ii,'dnstream_num_addU1'} = num_addU1;
end

end

