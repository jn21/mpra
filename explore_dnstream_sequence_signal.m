function explore_dnstream_sequence_signal(signal,e_or_p_ratio,only_good_idx)
%
% Input
%       signal: (string) field from mpra_data_table: Eg dnstream_addPAS
%       e_or_p_ratio: (string) which ratio to explore ('E_ratio_avg_rep' or
%           'P_ratio_avg_rep')
%       only_good_idx = remove reverse constructs?

mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');

signal_idx = logical(mpra_data{:,['dnstream_' signal]});
signal_dnstream_ids = unique(mpra_data{signal_idx,'dnstream_full_id'});
corresponding_dnstream_ids = cellfun(@(s) strrep(s,['_' signal],''),signal_dnstream_ids,...
    'Uni',false);

full_signal_table = table;
full_corr_table = table;

for ii = 1:length(signal_dnstream_ids)
    
    %Get barcode data for constructs that have the signal and
    %their counterparts
    temp_signal_idx = strcmp(mpra_data{:,'dnstream_full_id'},signal_dnstream_ids{ii});
    temp_corr_idx = strcmp(mpra_data{:,'dnstream_full_id'},corresponding_dnstream_ids{ii});
    temp_signal_table = mpra_data(temp_signal_idx,{'dnstream_full_id','upstream_full_id',e_or_p_ratio});
    temp_corr_table = mpra_data(temp_corr_idx,{'dnstream_full_id','upstream_full_id',e_or_p_ratio});
    
    % remove infinite rows
    signal_is_finite_idx = isfinite(temp_signal_table{:,e_or_p_ratio});
    temp_signal_table = temp_signal_table(signal_is_finite_idx,:);
    corr_is_finite_idx = isfinite(temp_corr_table{:,e_or_p_ratio});
    temp_corr_table = temp_corr_table(corr_is_finite_idx,:);
    
    %subset to common upstream ids
    common_upstream_ids = intersect(temp_signal_table{:,'upstream_full_id'},temp_corr_table{:,'upstream_full_id'});
    temp_idx1 = ismember(temp_signal_table{:,'upstream_full_id'},common_upstream_ids);
    temp_idx2 = ismember(temp_corr_table{:,'upstream_full_id'},common_upstream_ids);
    
    temp_signal_table = temp_signal_table(temp_idx1,:);
    temp_corr_table = temp_corr_table(temp_idx2,:);
    
    temp_signal_table = sortrows(temp_signal_table,2);
    temp_corr_table = sortrows(temp_corr_table,2);
    
    full_signal_table = vertcat(full_signal_table,temp_signal_table);
    full_corr_table = vertcat(full_corr_table,temp_corr_table);
    
    assert(isequal(temp_signal_table{:,'upstream_full_id'},temp_corr_table{:,'upstream_full_id'}),...
        'error - ids not equal');
end

figure;
subplot(2,2,1)
scatter(full_corr_table{:,e_or_p_ratio},full_signal_table{:,e_or_p_ratio})
xlabel('no_signal','Interpreter','None')
ylabel(signal,'Interpreter','None')
title(e_or_p_ratio,'Interpreter','None')
ax = gca;
max_lim = max([ax.XLim ax.YLim]);
min_lim = min([ax.XLim ax.YLim]);
grid on

ax.XLim = [min_lim max_lim];
ax.YLim = [min_lim max_lim];
hold on
plot([min_lim max_lim], [min_lim max_lim], 'r:')


subplot(2,2,2)
boxplot([full_signal_table{:,e_or_p_ratio} full_corr_table{:,e_or_p_ratio}],...
    'labels',{signal,'no_signal'})

[p,h,stats] = signrank(full_signal_table{:,e_or_p_ratio},full_corr_table{:,e_or_p_ratio})

subplot(2,2,3)
histogram(full_signal_table{:,e_or_p_ratio} - full_corr_table{:,e_or_p_ratio},...
    'BinWidth',.1)
xlabel('signal - no_signal','Interpreter','None')
ylabel('histogram')


end

