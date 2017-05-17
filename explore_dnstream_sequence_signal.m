function res = explore_dnstream_sequence_signal(activity_type,signal,remove_up_reverse,remove_dn_reverse)
%
% Input
%       signal: (string) field from mpra_data_table: Eg dnstream_addPAS
%       activity_type: (string) which ratio to explore ('E_ratio_avg_rep' or
%           'P_ratio_avg_rep')
%       remove_reverse = (boolean) remove reverse constructs

mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot_normalized.txt','Delimiter','\t');

%remove infinities
isfinite_idx = isfinite(mpra_data{:,'enhancer_activity'}) & isfinite(mpra_data{:,'promoter_activity'});
mpra_data = mpra_data(isfinite_idx,:);

if remove_up_reverse
    mpra_data = subset_table(mpra_data,'upstream_is_reverse',0);
end

if remove_dn_reverse
    mpra_data = subset_table(mpra_data,'dnstream_is_reverse',0);
end

signal_idx = logical(mpra_data{:,['dnstream_' signal]});
signal_dnstream_ids = unique(mpra_data{signal_idx,'dnstream_full_id'});
corresponding_dnstream_ids = cellfun(@(s) strrep(s,['_' signal],''),signal_dnstream_ids,...
    'Uni',false);

full_signal_table = table;
full_corr_table = table;

for ii = 1:length(signal_dnstream_ids)
    
    %Get barcode data for constructs that have the signal and
    %their counterparts
    temp_signal_table = subset_table(mpra_data,'dnstream_full_id',signal_dnstream_ids{ii});
    temp_corr_table = subset_table(mpra_data,'dnstream_full_id',corresponding_dnstream_ids{ii});
    
    % remove infinite rows
    signal_is_finite_idx = isfinite(temp_signal_table{:,activity_type});
    temp_signal_table = temp_signal_table(signal_is_finite_idx,:);
    corr_is_finite_idx = isfinite(temp_corr_table{:,activity_type});
    temp_corr_table = temp_corr_table(corr_is_finite_idx,:);
    
    %subset to common upstream ids
    common_upstream_ids = intersect(temp_signal_table{:,'upstream_full_id'},temp_corr_table{:,'upstream_full_id'});

    temp_signal_table = subset_table(temp_signal_table,'upstream_full_id',common_upstream_ids);
    temp_corr_table = subset_table(temp_corr_table,'upstream_full_id',common_upstream_ids);
    
    temp_signal_table = sortrows(temp_signal_table, 'upstream_full_id');
    temp_corr_table = sortrows(temp_corr_table, 'upstream_full_id');
    
    full_signal_table = vertcat(full_signal_table,temp_signal_table);
    full_corr_table = vertcat(full_corr_table,temp_corr_table);
    
    assert(isequal(temp_signal_table{:,'upstream_full_id'},temp_corr_table{:,'upstream_full_id'}),...
        'error - ids not equal');
end

diff_data = full_signal_table{:,activity_type} - full_corr_table{:,activity_type};
up_ids = full_signal_table{:,'upstream_full_id'};
dn_ids = full_signal_table{:,'dnstream_full_id'};
pval = signrank(diff_data);

res = struct('diff_data',diff_data,...
    'signal_data',full_signal_table{:,activity_type},...
    'nonsignal_data',full_corr_table{:,activity_type},...
    'up_ids',up_ids,...
    'dn_ids',dn_ids,...
    'pval',pval);

% switch e_or_p_ratio
%     case 'E_ratio_avg_rep'
%         e_or_p_ratio_str = 'Promoter Activity (E Barcode Ratio)';
%     case 'P_ratio_avg_rep'
%         e_or_p_ratio_str = 'Enhancer Activity (P Barcode Ratio)';
% end

% f = figure;
% set(f,'Units','normalized')
% set(f,'Position',[0 0 1 1])
% subplot(2,2,1)
% scatter(full_corr_table{:,e_or_p_ratio},full_signal_table{:,e_or_p_ratio})
% xlabel('no_signal','Interpreter','None')
% ylabel(signal,'Interpreter','None')
% %title(e_or_p_ratio,'Interpreter','None')
% ax = gca;
% max_lim = max([ax.XLim ax.YLim]);
% min_lim = min([ax.XLim ax.YLim]);
% grid on
% 
% ax.XLim = [min_lim max_lim];
% ax.YLim = [min_lim max_lim];
% hold on
% plot([min_lim max_lim], [min_lim max_lim], 'r:')
% 
% subplot(2,2,2)
% boxplot([full_signal_table{:,e_or_p_ratio} full_corr_table{:,e_or_p_ratio}],...
%     'labels',{signal,'no_signal'})
% ylabel(e_or_p_ratio_str)
% grid on
% 
% [pval,~,~] = signrank(full_signal_table{:,e_or_p_ratio},full_corr_table{:,e_or_p_ratio});
% 
% subplot(2,2,3)
% data = full_signal_table{:,e_or_p_ratio} - full_corr_table{:,e_or_p_ratio};
% med_data = median(data);
% histogram(data,...
%     'BinWidth',.1)
% hold on
% ax = gca;
% plot([med_data med_data], [0 ax.YLim(2)])
% xlabel([signal ' - no_signal'],'Interpreter','None')
% ylabel('histogram')
% 
% title_str = sprintf(['%s \n',...
%     '%s \n',...
%     '%d constructs \n',...
%     'median absolute difference = %.2g \n', ...
%     'sign rank pvalue = %.2g'],...
%     signal,e_or_p_ratio_str,length(data),med_data,pval);
% suptitle(title_str)
% 
% if save_fig
%     saveas(f,fullfile('~/Documents/mpra/fig/dnstream_signal/',[signal '_' e_or_p_ratio]),'png');
%     %print(f,'-r0','-dpng',fullfile('~/Documents/mpra/fig/dnstream_signal/',[signal '_' e_or_p_ratio]))
% end

end

