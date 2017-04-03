function T = explore_which_upstream_are_significant

mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');

%e_or_p_ratio = 'E_ratio_avg_rep';

res = struct;

%% Build null distribution
up_is_reverse_idx = logical(mpra_data{:,'upstream_is_reverse'});
null_table = mpra_data(up_is_reverse_idx,:);

%% 
unique_upstream_ids = unique(mpra_data{:,'upstream_full_id'});

% pval_table = horzcat(cell2table(cell(length(unique_up_non_reverse_ids),1)),array2table(zeros(length(unique_up_non_reverse_ids),2)));
% pval_table.Properties.VariableNames = {'upstream_full_id','nominal_pval','q_val'};

for ii = 1:length(unique_upstream_ids)
    
    this_table_id_idx = strcmp(mpra_data{:,'upstream_full_id'},unique_upstream_ids{ii});
    
    %get common dnstream id's
    this_dnstream_ids = mpra_data{this_table_id_idx,'dnstream_full_id'};
    common_dnstream_ids = intersect(this_dnstream_ids,null_table{:,'dnstream_full_id'});
    
    %Subset appropriately
    this_table = subset_table(mpra_data,'upstream_full_id',unique_upstream_ids{ii});
    this_table = subset_table(this_table,'dnstream_full_id',common_dnstream_ids);
    this_null_table = subset_table(null_table,'dnstream_full_id',common_dnstream_ids);
    
    %test
    pval_e_ratio = ranksum(this_table{:,'E_ratio_avg_rep'},...
        this_null_table{:,'E_ratio_avg_rep'},...
        'tail','right');
    
    pval_p_ratio = ranksum(this_table{:,'P_ratio_avg_rep'},...
        this_null_table{:,'P_ratio_avg_rep'},...
        'tail','right');
    
    %get relation to tss
    s = strsplit(unique_upstream_ids{ii},'_');
    relation_to_tss = strcat(s(2),'_',s(3));
    
    %store
    res(ii).upstream_full_id = unique_upstream_ids{ii};
    res(ii).relation_to_tss = relation_to_tss;
    res(ii).median_e_ratio = median(this_table{:,'E_ratio_avg_rep'});
    res(ii).median_p_ratio = median(this_table{:,'P_ratio_avg_rep'});
    res(ii).e_ratio_pval = pval_e_ratio;
    res(ii).p_ratio_pval = pval_p_ratio;
    
end

%% Plot
%scatter
T = struct2table(res);
figure
scatter(T{:,'median_e_ratio'},T{:,'median_p_ratio'})
xlabel('median E ratio (Promoter Activity)')
ylabel('median P ratio (Enhancer Activity)')
title('Promoter and Enhancer Activity of Upstream Sequences')


p_ratio_neg_control_median = median(null_table{:,'P_ratio_avg_rep'});
e_ratio_neg_control_median = median(null_table{:,'E_ratio_avg_rep'});

hold on
ax = gca;
plot([e_ratio_neg_control_median e_ratio_neg_control_median], ax.YLim, 'r:')
hold on
plot(ax.XLim, [p_ratio_neg_control_median p_ratio_neg_control_median], 'r:')

%boxplot
figure
subplot(1,2,1)
boxplot(T{:,'median_e_ratio'},T{:,'relation_to_tss'})
title('E Barcode Ratio (Promoter Activity)')

subplot(1,2,2)
boxplot(T{:,'median_p_ratio'},T{:,'relation_to_tss'})
title('P Barcode Ratio (Enhancer Activity)')

%% Full boxplot
