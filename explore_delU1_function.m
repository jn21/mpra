function res = explore_delU1_function(activity_type,num_del_u1,remove_up_reverse)

mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot_normalized.txt','Delimiter','\t');

%remove infinities
isfinite_idx = isfinite(mpra_data{:,'enhancer_activity'}) & isfinite(mpra_data{:,'promoter_activity'});
mpra_data = mpra_data(isfinite_idx,:);

if remove_up_reverse
    mpra_data = subset_table(mpra_data,'upstream_is_reverse',0);
end

mpra_data_no_addPAS = subset_table(mpra_data,'dnstream_addPAS',0);

%% U1 deletion
U1_del = subset_table(mpra_data_no_addPAS,'dnstream_num_delU1',num_del_u1);

U1_del_corresponding_constructs = cellfun(@(s) strrep(s,'_delU1',''),U1_del{:,'construct_no_oligo_id'},...
    'Uni',false);   

[~,locb] = ismember(U1_del_corresponding_constructs,mpra_data_no_addPAS{:,'construct_no_oligo_id'});

zero_del_correspond_to_U1_del = subset_table(mpra_data_no_addPAS,'construct_no_oligo_id',U1_del_corresponding_constructs);

U1_del = U1_del(logical(locb),:);

zero_del_correspond_to_U1_del = sortrows(zero_del_correspond_to_U1_del,{'dnstream_full_id','upstream_full_id'});
U1_del = sortrows(U1_del,{'dnstream_full_id','upstream_full_id'});

diff_data = U1_del{:,activity_type} - zero_del_correspond_to_U1_del{:,activity_type};
up_ids = U1_del{:,'upstream_full_id'};
dn_ids = U1_del{:,'dnstream_full_id'};
pval = signrank(diff_data);

res = struct('diff_data',diff_data,...
    'up_ids',up_ids,...
    'dn_ids',dn_ids,...
    'pval',pval);

%     size(zero_del_correspond_to_U1_del)
%     size(U1_del)
%     histogram(zero_del_correspond_to_U1_del{:,'E_ratio_avg_rep'} - U1_del{:,'E_ratio_avg_rep'})
%     median(zero_del_correspond_to_U1_del{:,'E_ratio_avg_rep'} - U1_del{:,'E_ratio_avg_rep'})

% mk_dnstream_signal_plot(U1_del{:,activity_type},...
%     zero_del_correspond_to_U1_del{:,activity_type},...
%     sprintf('Number of U1 deletions: %d',num_u1_del(ii)),...
%     jj,...
%     false)
%     end
% end

end
%% 
% delU1_dnstream_idx = logical(mpra_data{:,'dnstream_num_delU1'});
% delU1_dnstream_prefix = unique(mpra_data{delU1_dnstream_idx,'dnstream_prefix'});
% 
% %% Utility
% upstream_reverse_idx = logical(mpra_data{:,'upstream_is_reverse'});
% dnstream_reverse_idx = logical(mpra_data{:,'dnstream_is_reverse'});
% dnstream_add_pas_idx = logical(mpra_data{:,'dnstream_addPAS'});
% dnstream_del_pas_idx = logical(mpra_data{:,'dnstream_delPAS'});
% dnstream_strong_pas_idx = logical(mpra_data{:,'dnstream_addStrongPAS'});
% good_idx = ~(upstream_reverse_idx | ...
%              dnstream_reverse_idx | ...
%              dnstream_add_pas_idx | ...
%              dnstream_del_pas_idx | ...
%              dnstream_strong_pas_idx);
%          
% %% Examine all id's
% good_table = mpra_data(good_idx,:);
% 
% zero_del_idx = (good_table{:,'dnstream_num_delU1'} == 0);
% one_del_idx = (good_table{:,'dnstream_num_delU1'} == 1);
% two_del_idx = (good_table{:,'dnstream_num_delU1'} == 2);
%     
% dnstream_prefix_to_use = intersect(good_table{one_del_idx,'dnstream_prefix'},good_table{two_del_idx,'dnstream_prefix'});
% [~, dnstream_to_use_idx] = ismember(good_table{:,'dnstream_prefix'},dnstream_prefix_to_use);
% dnstream_to_use_idx = logical(dnstream_to_use_idx);
% 
% up_id_zero = good_table{zero_del_idx,'upstream_full_id'};
% up_id_one = good_table{one_del_idx,'upstream_full_id'};
% up_id_two = good_table{two_del_idx,'upstream_full_id'};
% 
% up_id_to_use = intersect(up_id_zero,intersect(up_id_one,up_id_two));
% [~, up_id_to_use_idx] = ismember(good_table{:,'upstream_full_id'},up_id_to_use);
% up_id_to_use_idx = logical(up_id_to_use_idx);
% 
% idx_to_use = dnstream_to_use_idx & up_id_to_use_idx ;
% 
% to_use_table = good_table(idx_to_use,:);
% boxplot(to_use_table{:,'P_ratio_avg_rep'},to_use_table{:,'dnstream_num_delU1'});
% xlabel('Number of U1 deletions')
% %ylabel('Promoter Activity (E Barcode ratio)')

%% Examine individual dnstream id's
% for ii = 1:length(delU1_dnstream_prefix)
%     this_id_idx = strcmp(delU1_dnstream_prefix(ii),mpra_data{:,'dnstream_prefix'});
%     this_table = mpra_data(this_id_idx & good_idx,:);
%     
%     num_U1_del = length(unique(this_table{:,'dnstream_num_delU1'}));
%     
%     %Subset to upstream ids that appear at all del U1 levels
%     up_ids = this_table{:,'upstream_full_id'};
%     appearances = countmember(unique(up_ids),up_ids);
%     up_ids_to_use = unique(up_ids(appearances == num_U1_del));
%     
%     [~, up_idx] = ismember(this_table{:,'upstream_full_id'},up_ids_to_use);
%     this_table = this_table(logical(up_idx),:);
%     
%     if num_U1_del >= 3
%         figure
%         boxplot(this_table{:,'E_ratio_avg_rep'},...
%             this_table{:,'dnstream_num_delU1'},...
%             'labels',this_table{:,'dnstream_num_delU1'});
%         title_str = sprintf('%s \n number of upstream ids = %d',...
%             delU1_dnstream_prefix{ii},length(up_ids_to_use));
%         title(title_str,'Interpreter','none');
%         xlabel('number U1 deletions')
%         ylabel('E Barcode ratio (Promoter Activity)')
%         ax = gca;
%         gca.XLabelRotation = 45;
%         saveas(gcf,fullfile('~/Documents/mpra/fig/U1_gt_2_del',delU1_dnstream_prefix{ii}),'png');
%     end
%     
%     
% end