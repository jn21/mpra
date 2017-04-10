
activity_type = {'E_ratio_avg_rep','P_ratio_avg_rep'};

%activity_table = table;

signals = {'addPAS','addStrongPAS','delPAS','addU1'};

for kk = 1:length(activity_type)
    activity_table = table;
    
    for ii = 1:length(signals)
        if strcmp(signals{ii},'addStrongPAS')
            remove_reverse = false;
        else
            remove_reverse = true;
        end

        [this_diff,this_pval] = explore_dnstream_sequence_signal(activity_type{kk},signals{ii},remove_reverse);

        %this_struct = struct('diff_data',this_diff,'signal',repmat(signals(ii),1,length(this_diff)));

        this_table = table(this_diff,repmat(signals(ii),length(this_diff),1),...
            'VariableNames',{'diff','signal'});

        activity_table = vertcat(activity_table,this_table);
        
        this_result_str = sprintf([...
            'n = %d \n',...
            'median = %.2g \n',...
            'pval = %.2g'],...
            length(this_diff),median(this_diff),this_pval);
        
        result_struct(ii).signal = signals{ii};
        result_struct(ii).string = this_result_str;
    end

    for jj = 1:3
        [this_diff,this_pval] = explore_delU1_function(activity_type{kk},jj);

        temp_str = ['delU1' strjoin(repmat({'_delU1'},1,jj-1),'')];
        this_signal = repmat({temp_str},length(this_diff),1);
        
        this_table = table(this_diff,this_signal,...
            'VariableNames',{'diff','signal'});

        activity_table = vertcat(activity_table,this_table);
        
        this_result_str = sprintf([...
            'n = %d \n',...
            'median = %.2g \n',...
            'pval = %.2g'],...
            length(this_diff),median(this_diff),this_pval);
        
        result_struct(ii+jj).signal = temp_str;
        result_struct(ii+jj).string = this_result_str;
    end
    
    switch activity_type{kk}
        case 'E_ratio_avg_rep'
            title_str = 'Promoter Activity (Construct Barcode Ratio)';
            offset = .1;
        case 'P_ratio_avg_rep'
            title_str = 'Enhancer Activity (GFP Barcode Ratio)';
            offset = .025;
    end

    f = figure;
    
%     label_order = {'addPAS','addStrongPAS','delPAS','delU1','delU1_delU1','delU1_delU1_delU1','addU1'}
%     {result_struct.signal}
    label_order = fliplr(1:7);
    %[~,pos_idx] = ismember({result_struct.signal},label_order);
    boxplot(activity_table{:,'diff'},activity_table{:,'signal'},...
        'orientation','horizontal',...
        'positions',label_order)
    title(title_str,'Interpreter','None')
    xlabel('Difference')
    ax = gca;
    hold on
    plot([0 0],ax.YLim, 'k:')
    %ax.Position = [0 0 .8 1];
    
    %ax2 = axes('Position',get(ax,'Position'),'YAxisLocation', 'right', 'Color','none','XTick',[], 'YTickLabel','');
    ax2 = axes('YAxisLocation', 'right', 'Color','none','XTick',[], 'YTickLabel','');

    linkaxes([ax ax2],'xy')
    
    axes(ax2)
    %offset = .1;
    %idx = ismember({result_struct.signal},label_order)
    text(zeros(ii+jj,1)+ax.XLim(2)+offset,ax.YTick,{result_struct(label_order).string})
    
    
%     [hx, hy] = format_ticks(ax2,'',{result_struct.string})
%     set(hy,'YAxis','Location','right')
end




% function mk_dnstream_signal_plot(signal_data,no_signal_data,signal_str,activity_str,save_fig)
% 
% switch activity_str{1}
%     case 'E_ratio_avg_rep'
%         activity_str = 'Promoter Activity (E Barcode Ratio)';
%     case 'P_ratio_avg_rep'
%         activity_str = 'Enhancer Activity (P Barcode Ratio)';
% end
% 
% f = figure;
% set(f,'Units','normalized')
% set(f,'Position',[0 0 1 1])
% subplot(2,2,1)
% scatter(no_signal_data,signal_data)
% xlabel('no_signal','Interpreter','None')
% ylabel(signal_str,'Interpreter','None')
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
% boxplot([signal_data no_signal_data],...
%     'labels',{signal_str,'no_signal'})
% ylabel(activity_str)
% grid on
% 
% [pval,~,~] = signrank(signal_data,no_signal_data);
% 
% subplot(2,2,3)
% diff_data = signal_data - no_signal_data;
% med_diff_data = median(diff_data);
% histogram(diff_data,...
%     'BinWidth',.1)
% hold on
% ax = gca;
% plot([med_diff_data med_diff_data], [0 ax.YLim(2)])
% xlabel([signal_str ' - no_signal'],'Interpreter','None')
% ylabel('histogram')
% 
% title_str = sprintf(['%s \n',...
%     '%s \n',...
%     '%d constructs \n',...
%     'median absolute difference = %.2g \n', ...
%     'sign rank pvalue = %.2g'],...
%     signal_str,activity_str,length(diff_data),med_diff_data,pval);
% suptitle(title_str)
% 
% if save_fig
%     saveas(f,fullfile('~/Documents/mpra/fig/dnstream_signal/',[signal_str '_' activity_str]),'png');
%     %print(f,'-r0','-dpng',fullfile('~/Documents/mpra/fig/dnstream_signal/',[signal '_' e_or_p_ratio]))
% end
% 
% 
% end

