
activity_type = {'promoter_activity','enhancer_activity'};

figdir = '~/Documents/mpra/fig/downstream_signal';

%signals = {'addPAS','addStrongPAS','delPAS','addU1'};
signals = {'addPAS','addStrongPAS','addU1'};

run_anova = false;

for kk = 1:length(activity_type)
    activity_table = table;
    
    switch activity_type{kk}
        case 'promoter_activity'
            activity_str = 'Promoter Activity';
            offset = .1;
        case 'enhancer_activity'
            activity_str = 'Enhancer Activity';
            offset = .025;
    end
    
    for ii = 1:length(signals)
        if strcmp(signals{ii},'addStrongPAS')
            remove_dn_reverse = false;
            remove_up_reverse = true;
        else
            remove_dn_reverse = false;
            remove_up_reverse = true;
        end

        %[this_diff,this_pval] = explore_dnstream_sequence_signal(activity_type{kk},signals{ii},remove_reverse);
        res = explore_dnstream_sequence_signal(activity_type{kk},signals{ii},remove_up_reverse,remove_dn_reverse);

        %this_struct = struct('diff_data',this_diff,'signal',repmat(signals(ii),1,length(this_diff)));

        this_diff = res(1).diff_data;
        this_table = table(this_diff,repmat(signals(ii),length(this_diff),1),...
            'VariableNames',{'diff','signal'});

        %this signal's anova
        if run_anova
            figure
            [anova_p, anova_tab] = anova1(this_diff,{res.dn_ids},'off');
            anova_table = cell2table(anova_tab);
            anova_var_explained = cell2mat(anova_table{2,2})/cell2mat(anova_table{4,2});

            boxplot(this_diff,{res.dn_ids},...
                'OutlierSize',1,...
                'symbol','r.')
            title(sprintf('%s \n Signal: %s \n Variance Explained (ANOVA) = %.2g %%',activity_str,signals{ii},anova_var_explained*100),...
                'Interpreter','none')
            ax = gca;
            ax.XTickLabel = '';
            ylabel('Signal Effect')
            xlabel('Downstream Sequence')

            fig_str = [activity_type{kk} '_' signals{ii} '_anova'];
            saveas(gcf,fullfile(figdir,fig_str),'png')
        end
        
%         figure;
%         scatter(res(1).nonsignal_data,res(1).signal_data)
%         xlabel('Native Sequence')
%         ylabel('Modified with Signal')
%         title(sprintf('%s \n %s',activity_str,signals{ii}),'Interpreter','none')
%         grid on
%         ax = gca;
%         ax.XLim = [-8 0];
%         ax.YLim = [-8 0];
%         hold on
%         plot([-8 0], [-8 0], 'r:')
        
        activity_table = vertcat(activity_table,this_table);
        
        this_result_str = sprintf([...
            'n = %d \n',...
            'median = %.2g \n',...
            'p = %.2g'],...
            length(this_diff),median(this_diff),res(1).pval);
        
        result_struct(ii).signal = signals{ii};
        result_struct(ii).string = this_result_str;
    end

    %U1 Deletions
    for jj = 1%:3
        %[this_diff,this_pval] = explore_delU1_function(activity_type{kk},jj);
        res = explore_delU1_function(activity_type{kk},jj,true);
        this_diff = res(1).diff_data;

        temp_str = ['delU1' strjoin(repmat({'_delU1'},1,jj-1),'')];
        this_signal = repmat({temp_str},length(this_diff),1);
        
        %this signal's anova
        if run_anova
            figure
            [anova_p, anova_tab] = anova1(this_diff,{res.dn_ids},'off');
            anova_table = cell2table(anova_tab);
            anova_var_explained = cell2mat(anova_table{2,2})/cell2mat(anova_table{4,2});

            boxplot(this_diff,{res.dn_ids},...
                'OutlierSize',1,...
                'symbol','r.')
            title(sprintf('%s \n Signal: %s \n Variance Explained (ANOVA) = %.2g %%',activity_str,temp_str,anova_var_explained*100),...
                'Interpreter','none')
            ax = gca;
            ax.XTickLabels = '';
            ylabel('Signal Effect')
            xlabel('Downstream Sequence')

            fig_str = [activity_type{kk} '_' temp_str '_anova'];
            saveas(gcf,fullfile(figdir,fig_str),'png')
        end
        
        this_table = table(this_diff,this_signal,...
            'VariableNames',{'diff','signal'});

        activity_table = vertcat(activity_table,this_table);
        
        this_result_str = sprintf([...
            'n = %d \n',...
            'median = %.2g \n',...
            'p = %.2g'],...
            length(this_diff),median(this_diff),res(1).pval);
        
        result_struct(ii+jj).signal = temp_str;
        result_struct(ii+jj).string = this_result_str;
    end

    %activity_table = subset_table(activity_table,'signal',{'addPAS','addStrongPAS','addU1','delU1'});
    f = figure;
    
    f.PaperPositionMode = 'auto';
    f.Units = 'Normalized';
    f.OuterPosition = [0 0 .85 1];
    f.Position = [0 0 .6 1];
    %set(gca, 'LooseInset', [.13 .2 .095 .075]);
    
%     label_order = {'addPAS','addStrongPAS','delPAS','delU1','delU1_delU1','delU1_delU1_delU1','addU1'}
%     {result_struct.signal}
    
    %label_order = fliplr(1:7);
    label_order = fliplr(1:4);
    %[~,pos_idx] = ismember({result_struct.signal},label_order);
    boxplot(activity_table{:,'diff'},activity_table{:,'signal'},...
        'orientation','horizontal',...
        'positions',label_order)
    xlabel('Signal Effect')
    ax = gca;
    ax.FontSize = 18;
    hold on
    plot([0 0],[ax.YLim(1) - .5 ax.YLim(2)], 'k:')
    title(activity_str,'Interpreter','None','FontSize',24)
    %ax.Position = [0 0 .8 1];
    
    %ax2 = axes('Position',get(ax,'Position'),'YAxisLocation', 'right', 'Color','none','XTick',[], 'YTickLabel','');
    ax2 = axes('YAxisLocation', 'right', 'Color','none','XTick',[], 'YTickLabel','');
    
    linkaxes([ax ax2],'xy')
    
    axes(ax2)
    %offset = .1;
    %idx = ismember({result_struct.signal},label_order)
    text(zeros(ii+jj,1)+ax.XLim(2)+offset,ax.YTick,{result_struct(label_order).string},...
        'FontSize',14)

    fig_str = [activity_type{kk} '_downstream_signal_summary'];
    %saveas(f,fullfile(figdir,fig_str),'png')
    
    
%     [hx, hy] = format_ticks(ax2,'',{result_struct.string})
%     set(hy,'YAxis','Location','right')
end


