% Explore the depletion of 'E' barcodes under an downstream construct having
% a strong PAS signal

%
mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');

%Get indices of constructs with add strong PAS annotation
strong_pas_idx = logical(mpra_data{:,'dnstream_addStrongPAS'});
%strong_pas_idx = logical(mpra_data{:,'dnstream_num_delU1'});
%strong_pas_idx = logical(mpra_data{:,'dnstream_num_addU1'});
%strong_pas_idx = logical(mpra_data{:,'dnstream_addPAS'});
%strong_pas_idx = logical(mpra_data{:,'dnstream_delPAS'});

add_strong_pas_dnstream_ids = unique(mpra_data{strong_pas_idx,'dnstream_full_id'});

%funny hack to get the corresponding dnstream id without addStrongPAS. Just
%delete '_addStrongPAS' from the dnstream_id
 corresponding_dnstream_ids = cellfun(@(s) strrep(s,'_addStrongPAS',''),add_strong_pas_dnstream_ids,...
     'Uni',false);
% corresponding_dnstream_ids = cellfun(@(s) strrep(s,'_delU1',''),add_strong_pas_dnstream_ids,...
%     'Uni',false);
% corresponding_dnstream_ids = cellfun(@(s) strrep(s,'_addU1',''),add_strong_pas_dnstream_ids,...
%     'Uni',false);
% corresponding_dnstream_ids = cellfun(@(s) strrep(s,'_addPAS',''),add_strong_pas_dnstream_ids,...
%     'Uni',false);
% corresponding_dnstream_ids = cellfun(@(s) strrep(s,'_delPAS',''),add_strong_pas_dnstream_ids,...
%     'Uni',false);

plot_val = 'E_ratio_avg_rep';

full_pas_table = table;
full_corr_table = table;

%% Make a diagnostic plot for each dnstream id
for ii = 1:length(add_strong_pas_dnstream_ids)
    
    %Get barcode data for constructs that have an addStrongPAS and
    %their counterparts
    temp_pas_idx = strcmp(mpra_data{:,'dnstream_full_id'},add_strong_pas_dnstream_ids{ii});
    temp_corr_idx = strcmp(mpra_data{:,'dnstream_full_id'},corresponding_dnstream_ids{ii});
    temp_pas_table = mpra_data(temp_pas_idx,{'dnstream_full_id','upstream_full_id',plot_val});
    temp_corr_table = mpra_data(temp_corr_idx,{'dnstream_full_id','upstream_full_id',plot_val});
    
    % remove infinite rows
    pas_is_finite_idx = isfinite(temp_pas_table{:,plot_val});
    temp_pas_table = temp_pas_table(pas_is_finite_idx,:);
    corr_is_finite_idx = isfinite(temp_corr_table{:,plot_val});
    temp_corr_table = temp_corr_table(corr_is_finite_idx,:);
    
    %subset to common upstream ids
    common_upstream_ids = intersect(temp_pas_table{:,'upstream_full_id'},temp_corr_table{:,'upstream_full_id'});
    temp_idx1 = ismember(temp_pas_table{:,'upstream_full_id'},common_upstream_ids);
    temp_idx2 = ismember(temp_corr_table{:,'upstream_full_id'},common_upstream_ids);
    
    temp_pas_table = temp_pas_table(temp_idx1,:);
    temp_corr_table = temp_corr_table(temp_idx2,:);
    
    temp_pas_table = sortrows(temp_pas_table,2);
    temp_corr_table = sortrows(temp_corr_table,2);
    
    full_pas_table = vertcat(full_pas_table,temp_pas_table);
    full_corr_table = vertcat(full_corr_table,temp_corr_table);
    
    assert(isequal(temp_pas_table{:,'upstream_full_id'},temp_corr_table{:,'upstream_full_id'}),...
        'error - ids not equal');
    
    
    
    %Plot
%     figure;
%     subplot(1,2,1)
%     scatter(temp_pas_table{:,plot_val},temp_corr_table{:,plot_val})
%     xlabel(add_strong_pas_dnstream_ids{ii},'Interpreter','None')
%     ylabel(corresponding_dnstream_ids{ii},'Interpreter','None')
%     title(plot_val,'Interpreter','None')
%     ax = gca;
%     max_lim = max([ax.XLim ax.YLim]);
%     min_lim = min([ax.XLim ax.YLim]);
%     
%     ax.XLim = [min_lim max_lim];
%     ax.YLim = [min_lim max_lim];
%     hold on
%     plot([min_lim max_lim], [min_lim max_lim], 'r:')
%     
%     
%     subplot(1,2,2)
%     boxplot([temp_pas_table{:,plot_val} temp_corr_table{:,plot_val}],...
%         'labels',{'addStrongPAS','no_PAS'})
%     
end

%% 
figure;
subplot(1,2,1)
scatter(full_pas_table{:,plot_val},full_corr_table{:,plot_val})
xlabel('addPas','Interpreter','None')
ylabel('no_addition','Interpreter','None')
title(plot_val,'Interpreter','None')
ax = gca;
max_lim = max([ax.XLim ax.YLim]);
min_lim = min([ax.XLim ax.YLim]);

ax.XLim = [min_lim max_lim];
ax.YLim = [min_lim max_lim];
hold on
plot([min_lim max_lim], [min_lim max_lim], 'r:')


subplot(1,2,2)
boxplot([full_pas_table{:,plot_val} full_corr_table{:,plot_val}],...
    'labels',{'addPAS','no_PAS'})

[p,h,stats] = signrank(full_pas_table{:,plot_val},full_corr_table{:,plot_val})

figure
histogram(full_pas_table{:,plot_val} - full_corr_table{:,plot_val},...
    'BinWidth',.1)
