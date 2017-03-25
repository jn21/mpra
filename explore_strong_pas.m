% Explore the depletion of 'E' barcodes under an downstream construct having
% a strong PAS signal

%
mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');

%Get indices of constructs with add strong PAS annotation
strong_pas_idx = logical(mpra_data{:,'dnstream_addStrongPAS'});
%strong_pas_idx = logical(mpra_data{:,'dnstream_num_addU1'});

add_strong_pas_dnstream_ids = unique(mpra_data{strong_pas_idx,'dnstream_id'});

%funny hack to get the corresponding dnstream id without addStrongPAS. Just
%delete '_addStrongPAS' from the dnstream_id
 corresponding_dnstream_ids = cellfun(@(s) strrep(s,'_addStrongPAS',''),add_strong_pas_dnstream_ids,...
     'Uni',false);
% corresponding_dnstream_ids = cellfun(@(s) strrep(s,'_addU1',''),add_strong_pas_dnstream_ids,...
%     'Uni',false);

%% Make a diagnostic plot for each dnstream id
for ii = 1:length(add_strong_pas_dnstream_ids)
    
    %Get 'E' barcode data for constructs that have an addStrongPAS and
    %their counterparts
    temp_pas_idx = strcmp(mpra_data{:,'dnstream_id'},add_strong_pas_dnstream_ids{ii});
    temp_corr_idx = strcmp(mpra_data{:,'dnstream_id'},corresponding_dnstream_ids{ii});
    temp_pas_table = mpra_data(temp_pas_idx,{'dnstream_id','upstream_id','E_ratio_avg_rep'});
    temp_corr_table = mpra_data(temp_corr_idx,{'dnstream_id','upstream_id','E_ratio_avg_rep'});
    
    %subset to common upstream ids
    common_upstream_ids = intersect(temp_pas_table{:,'upstream_id'},temp_corr_table{:,'upstream_id'});
    temp_idx1 = ismember(temp_pas_table{:,'upstream_id'},common_upstream_ids);
    temp_idx2 = ismember(temp_corr_table{:,'upstream_id'},common_upstream_ids);
    
    temp_pas_table = temp_pas_table(temp_idx1,:);
    temp_corr_table = temp_corr_table(temp_idx2,:);
    
    temp_pas_table = sortrows(temp_pas_table,2);
    temp_corr_table = sortrows(temp_corr_table,2);
    
    assert(isequal(temp_pas_table{:,'upstream_id'},temp_corr_table{:,'upstream_id'}),...
        'error - ids not equal');
    
    %Plot
    figure;
    subplot(1,2,1)
    scatter(temp_pas_table{:,'E_ratio_avg_rep'},temp_corr_table{:,'E_ratio_avg_rep'})
    xlabel(add_strong_pas_dnstream_ids{ii},'Interpreter','None')
    ylabel(corresponding_dnstream_ids{ii},'Interpreter','None')
    title('E_ratio_avg_rep','Interpreter','None')
    ax = gca;
    max_lim = max([ax.XLim ax.YLim]);
    min_lim = min([ax.XLim ax.YLim]);
    
    ax.XLim = [min_lim max_lim];
    ax.YLim = [min_lim max_lim];
    hold on
    plot([min_lim max_lim], [min_lim max_lim], 'r:')
    
    
    subplot(1,2,2)
    boxplot([temp_pas_table{:,'E_ratio_avg_rep'} temp_corr_table{:,'E_ratio_avg_rep'}],...
        'labels',{'addStrongPAS','no_PAS'})
    

end