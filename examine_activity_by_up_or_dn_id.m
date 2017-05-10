function T = examine_activity_by_up_or_dn_id(up_or_dn_id,activity_type,save_plot)

%% Setup variable names
if strcmp(up_or_dn_id,'up')
    reverse_str = 'upstream_is_reverse';
    id_str = 'upstream_full_id';
    other_id_str = 'dnstream_full_id';
    region_str = 'upstream_region_id';
    derived_str = 'upstream_sequence_derivation';
    id_str_for_title = 'Upstream Sequence Medians';
    id_str_for_label = 'Upstream Sequence';
elseif strcmp(up_or_dn_id,'dn')
    reverse_str = 'dnstream_is_reverse';
    id_str = 'dnstream_full_id';
    other_id_str = 'upstream_full_id';
    region_str = 'dnstream_region_id';
    derived_str = 'dnstream_sequence_derivation';
    id_str_for_title = 'Downstream Sequence Medians';
    id_str_for_label = 'Downstream Sequence';
else
    error('must be either up or dn')
end

if strcmp(activity_type,'enhancer_activity')
    %activity_type = 'P_ratio_avg_rep';
    activity_str = 'Enhancer Activity';
elseif strcmp(activity_type,'promoter_activity')
    %activity_type = 'E_ratio_avg_rep';
    activity_str = 'Promoter Activity';
else
    error('must be either up or dn')
end

res = struct;

%% Load and subset data
mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');
mpra_data = mpra_data(~mpra_data{:,'dnstream_is_modified'},:);
isfinite_idx = isfinite(mpra_data{:,activity_type});
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
    pval = ranksum(this_table{:,activity_type},...
        this_null_table{:,activity_type},...
        'tail','right');
    
    %Examine distributions
%     if mod(ii,9) == 3 || mod(ii,9) == 4 || mod(ii,9) == 5
%         figure
%         histogram(this_table{:,activity_type},'BinWidth',.05)
%     end
    
    %get relation to tss
%     s = strsplit(unique_ids{ii},'_');
%     relation_to_tss = strcat(s(2),'_',s(3));
    %derivation = this_table{1,derived_str};
    
    %store
    res(ii).full_id = unique_ids{ii};
    res(ii).derivation = this_table{1,derived_str};
    res(ii).is_reverse = this_table{1,reverse_str};
    res(ii).region_id = this_table{1,region_str};
    res(ii).median_ratio = median(this_table{:,activity_type});
    res(ii).neg_control_median_ratio = median(null_table{:,activity_type});
    res(ii).pval = pval;
    
end

[~,qval] = mafdr([res.pval]);

T = struct2table(res);

T{:,'qval'} = qval';

%% Anova
[p_anova,anova_tab,anova_stats] = anova1(mpra_data{:,activity_type},...
    mpra_data{:,id_str},...
    'off');
anova_table = cell2table(anova_tab);
anova_var_explained = cell2mat(anova_table{2,2})/cell2mat(anova_table{4,2});

%% Boxplot - medians by promoter/enhancer derived
%boxplot
bar_fig = figure;
boxplot(T{:,'median_ratio'},T{:,'derivation'})
T_no_reverse = subset_table(T,'is_reverse',0);
promoter_derived_data = subset_table(T_no_reverse,'derivation','promoter_derived');
enhancer_derived_data = subset_table(T_no_reverse,'derivation','enhancer_derived');

pval = ranksum(promoter_derived_data{:,'median_ratio'},enhancer_derived_data{:,'median_ratio'});

title_str = sprintf('%s \n %s \n p = %.2g',activity_str,id_str_for_title,pval);

title(title_str,'Interpreter','none')


bar_fig.PaperPositionMode = 'auto';
bar_fig.Units = 'Normalized';
bar_fig.Position = [0 0 .25 .6];
%fig.OuterPosition = [0 0 1 .65];
%set(gca, 'LooseInset', [.03 .03 .03 .08]);

if save_plot
    save_str = sprintf('~/Documents/mpra/fig/up_dn_analysis/bar_%s_%s',...
        up_or_dn_id,activity_type);
    saveas(bar_fig,save_str,'png')
end


%% Full boxplot - ANOVA
anova_fig = figure;

%get sort variable to plot reverse sequences first
reverse_temp = cellfun(@(s) strfind(s,'reverse'), unique_ids, 'uni',false);
prom_derived_temp = cellfun(@(s) strfind(s,'gt'), unique_ids, 'uni',false);
reverse_temp2 = cellfun(@(s) length(s), reverse_temp);
prom_derived_temp2 = cellfun(@(s) length(s), prom_derived_temp);

%[reverse_temp2, prom_derived_temp2]

[~,sort_idx] = sortrows([reverse_temp2, prom_derived_temp2],[-1 2]);

boxplot(mpra_data{:,activity_type},...
    mpra_data{:,id_str},...
    'grouporder',unique_ids(sort_idx),...
    'outliersize',1,...
    'symbol','r.')
f = gca;
f.Box = 'off';
%ax.XTick = [];
f.TickLength = [.01 .01];
f.XTickLabelRotation = 90;
f.XTickLabel = ' ';

title_str2 = sprintf('%s \n Group by: %s \n Variance Explained (ANOVA) = %d %%',...
    activity_str,id_str_for_label,round(100*anova_var_explained));
title(title_str2,'Interpreter','None','FontSize',15)
ylabel(activity_str,'FontSize',15)
xlabel(id_str_for_label,'Interpreter','none','FontSize',15)

anova_fig.PaperPositionMode = 'auto';
anova_fig.Units = 'Normalized';
anova_fig.Position = [0 0 1 .65];
%fig.OuterPosition = [0 0 1 .65];
set(gca, 'LooseInset', [.03 .03 .03 .08]);
    
if save_plot
    save_str = sprintf('~/Documents/mpra/fig/up_dn_analysis/anova_%s_%s',...
        up_or_dn_id,activity_type);
    saveas(anova_fig,save_str,'png')
end

end