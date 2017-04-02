mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');

%% 

delU1_dnstream_idx = logical(mpra_data{:,'dnstream_num_delU1'});
delU1_dnstream_prefix = unique(mpra_data{delU1_dnstream_idx,'dnstream_prefix'});


% Utility
upstream_reverse_idx = logical(mpra_data{:,'upstream_is_reverse'});
dnstream_reverse_idx = logical(mpra_data{:,'dnstream_is_reverse'});
dnstream_add_pas_idx = logical(mpra_data{:,'dnstream_addPAS'});
dnstream_del_pas_idx = logical(mpra_data{:,'dnstream_delPAS'});
dnstream_strong_pas_idx = logical(mpra_data{:,'dnstream_addStrongPAS'});
good_idx = ~(upstream_reverse_idx | ...
             dnstream_reverse_idx | ...
             dnstream_add_pas_idx | ...
             dnstream_del_pas_idx | ...
             dnstream_strong_pas_idx);

for ii = 1:20%length(delU1_dnstream_prefix)
    
    delU1_dnstream_prefix(ii)
    this_id_idx = strcmp(delU1_dnstream_prefix(ii),mpra_data{:,'dnstream_prefix'});
    this_table = mpra_data(this_id_idx & good_idx,:);
    
    num_U1_del = length(unique(this_table{:,'dnstream_num_delU1'}));
    
    %Subset to upstream ids that appear at all del U1 levels
    up_ids = this_table{:,'upstream_full_id'};
    appearances = countmember(unique(up_ids),up_ids);
    up_ids_to_use = unique(up_ids(appearances == num_U1_del));
    
    [~, up_idx] = ismember(this_table{:,'upstream_full_id'},up_ids_to_use);
    this_table = this_table(logical(up_idx),:);
    
    
    if num_U1_del >= 2
        figure
        boxplot(this_table{:,'E_ratio_avg_rep'},...
            this_table{:,'dnstream_num_delU1'},...
            'labels',this_table{:,'dnstream_num_delU1'});
        title_str = sprintf('%s \n number of upstream ids = %d',...
            delU1_dnstream_prefix{ii},length(up_ids_to_use));
        title(title_str,'Interpreter','none')
        xlabel('number U1 deletions')
        ylabel('E Barcode ratio (Promoter Activity)')
        ax = gca;
        gca.XLabelRotation = 45;
    end
    
    
end