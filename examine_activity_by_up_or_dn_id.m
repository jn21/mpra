function T = examine_activity_by_up_or_dn_id(up_or_dn_id,gfp_or_construct_ratio)

%% Setup variable names
if strcmp(up_or_dn_id,'up')
    reverse_str = 'upstream_is_reverse';
    id_str = 'upstream_full_id';
    other_id_str = 'dnstream_full_id';
    region_str = 'upstream_region_id';
elseif strcmp(up_or_dn_id,'dn')
    reverse_str = 'dnstream_is_reverse';
    id_str = 'dnstream_full_id';
    other_id_str = 'upstream_full_id';
    region_str = 'dnstream_region_id';
else
    error('must be either up or dn')
end

if strcmp(gfp_or_construct_ratio,'gfp')
    barcode_type = 'P_ratio_avg_rep';
    title_str = 'Enhancer Activity (GFP Barcode Ratio)';
elseif strcmp(gfp_or_construct_ratio,'construct')
    barcode_type = 'E_ratio_avg_rep';
    title_str = 'Promoter Activity (Construct Barcode Ratio)';
else
    error('must be either up or dn')
end

res = struct;

%% Load and subset data
mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');
mpra_data = mpra_data(~mpra_data{:,'dnstream_is_modified'},:);
isfinite_idx = isfinite(mpra_data{:,barcode_type});
mpra_data = mpra_data(isfinite_idx,:);

%% Build null distribution
up_is_reverse_idx = logical(mpra_data{:,reverse_str});
null_table = mpra_data(up_is_reverse_idx,:);

%% 
unique_ids = unique(mpra_data{:,id_str});

for ii = 1:length(unique_ids)
    
    this_table_id_idx = strcmp(mpra_data{:,id_str},unique_ids{ii});
    
    %get common dnstream id's
    this_dnstream_ids = mpra_data{this_table_id_idx,other_id_str};
    common_dnstream_ids = intersect(this_dnstream_ids,null_table{:,other_id_str});
    
    %Subset appropriately
    this_table = subset_table(mpra_data,id_str,unique_ids{ii});
    this_table = subset_table(this_table,other_id_str,common_dnstream_ids);
    this_null_table = subset_table(null_table,other_id_str,common_dnstream_ids);
    
    %test
    pval = ranksum(this_table{:,barcode_type},...
        this_null_table{:,barcode_type},...
        'tail','right');
    
    %get relation to tss
    s = strsplit(unique_ids{ii},'_');
    relation_to_tss = strcat(s(2),'_',s(3));
    
    %store
    res(ii).full_id = unique_ids{ii};
    res(ii).relation_to_tss = relation_to_tss{1};
    res(ii).is_reverse = this_table{1,reverse_str};
    res(ii).region_id = this_table{1,region_str};
    res(ii).median_ratio = median(this_table{:,barcode_type});
    res(ii).neg_control_median_ratio = median(null_table{:,barcode_type});
    res(ii).pval = pval;
    
end

T = struct2table(res);

%% Anova
[p_anova,anova_tab,anova_stats] = anova1(mpra_data{:,barcode_type},...
    mpra_data{:,id_str},...
    'off');
anova_table = cell2table(anova_tab);
anova_var_explained = cell2mat(anova_table{2,2})/cell2mat(anova_table{4,2});

%% Plot
%boxplot
figure
subplot(1,2,1)
boxplot(T{:,'median_ratio'},T{:,'relation_to_tss'})
title(title_str)

%% Full boxplot
figure;

%get sort variable to plot reverse sequences first
temp = cellfun(@(s) strfind(s,'reverse'), unique_ids, 'uni',false);
temp2 = cellfun(@(s) length(s), temp);
[~,sort_idx] = sort(temp2,'descend');

boxplot(mpra_data{:,barcode_type},...
    mpra_data{:,id_str},...
    'grouporder',unique_ids(sort_idx),...
    'outliersize',1,...
    'symbol','')
ax = gca;
ax.Box = 'off';
ax.XTick = [];
ax.XTickLabel = [];

title_str2 = sprintf('%s \n Variance Explained (ANOVA) = %d %%',...
    title_str,round(100*anova_var_explained));
title(title_str2)
ylabel(title_str)
xlabel(id_str,'Interpreter','none')

end