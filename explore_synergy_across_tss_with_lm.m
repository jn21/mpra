mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');
isfinite_idx = isfinite(mpra_data{:,'P_ratio_avg_rep'}) & isfinite(mpra_data{:,'E_ratio_avg_rep'});
mpra_data = mpra_data(isfinite_idx,:);

%% Run Regressions
[lm_table, mpra_table_for_lm] = convert_data_table_to_lm_table(mpra_data,true);
enhancer_activity = lm_table;
enhancer_activity.E_ratio = [];
promoter_activity = lm_table;
promoter_activity.P_ratio = [];

promoter_activity_lm = fitlm(promoter_activity,...
            'ResponseVar','E_ratio');
        
promoter_activity_lm_no_intercept = fitlm(promoter_activity,...
            'ResponseVar','E_ratio',...
            'Intercept',false);
        
enhancer_activity_lm = fitlm(enhancer_activity,...
    'ResponseVar','P_ratio');

%% Make Plots Comparing Residuals
same_region_idx = logical(mpra_table_for_lm{:,'is_intact_sequence'});
promoter_derived_idx = strcmp(mpra_table_for_lm{:,'upstream_sequence_relation_to_tss'},'wt_100nt');

%group_idx
% 1 : different region
% 2 : same region - enhancer derived
% 3 : same region - promoter derived
group_idx = 1 + same_region_idx .* (1 + promoter_derived_idx);

pval_promoter_12 = ranksum(promoter_activity_lm.Residuals{group_idx == 1,'Raw'},...
        promoter_activity_lm.Residuals{group_idx == 2,'Raw'});
pval_promoter_13 = ranksum(promoter_activity_lm.Residuals{group_idx == 1,'Raw'},...
        promoter_activity_lm.Residuals{group_idx == 3,'Raw'});
pval_promoter_23 = ranksum(promoter_activity_lm.Residuals{group_idx == 2,'Raw'},...
        promoter_activity_lm.Residuals{group_idx == 3,'Raw'});

pval_enhancer_12 = ranksum(enhancer_activity_lm.Residuals{group_idx == 1,'Raw'},...
        enhancer_activity_lm.Residuals{group_idx == 2,'Raw'});
pval_enhancer_13 = ranksum(enhancer_activity_lm.Residuals{group_idx == 1,'Raw'},...
        enhancer_activity_lm.Residuals{group_idx == 3,'Raw'});
pval_enhancer_23 = ranksum(enhancer_activity_lm.Residuals{group_idx == 2,'Raw'},...
        enhancer_activity_lm.Residuals{group_idx == 3,'Raw'});
    
% pval_enhancer_activity = ranksum(enhancer_activity_lm.Residuals{same_region_idx,'Raw'},...
%         enhancer_activity_lm.Residuals{~same_region_idx,'Raw'});

% promoter_activity_effect = median(promoter_activity_lm.Residuals{same_region_idx,'Raw'}) - ...
%     median(promoter_activity_lm.Residuals{~same_region_idx,'Raw'});
% 
% enhancer_activity_effect = median(enhancer_activity_lm.Residuals{same_region_idx,'Raw'}) - ...
%     median(enhancer_activity_lm.Residuals{~same_region_idx,'Raw'});


labels = {'Different Region',...
    sprintf('   Same Region \nEnhancer Derived'),...
    sprintf('   Same Region \nPromoter Derived')};

figure
subplot(1,2,1)
boxplot(promoter_activity_lm.Residuals{:,'Raw'},...
    group_idx)
title('Promoter Activity Residuals')
ylabel('Residuals')
grid on
ax = gca;
ax.XTickLabel = [];
text(ax.XTick - .35,repmat(ax.YLim(1),3,1) - .2,labels,'FontSize',12)
H = sigstar({[1,3], [2,3]},[ pval_promoter_13, pval_promoter_23]);

subplot(1,2,2)
boxplot(enhancer_activity_lm.Residuals{:,'Raw'},...
    group_idx)
title('Enhancer Activity Residuals')
grid on
ylabel('Residuals')
ax = gca;
ax.XTickLabel = [];
text(ax.XTick - .35,repmat(ax.YLim(1),3,1) - .05,labels,'FontSize',12)

% figure
% subplot(1,2,1)
% boxplot(promoter_activity_lm.Residuals{:,'Raw'},...
%     mpra_table_for_lm{:,'is_intact_sequence'},...
%     'Labels',{'Different Region','Same Region'})
% title(sprintf('Promoter Activity Residuals \n p = %.2g',pval_promoter_activity))
% grid on
% 
% subplot(1,2,2)
% boxplot(enhancer_activity_lm.Residuals{:,'Raw'},...
%     mpra_table_for_lm{:,'is_intact_sequence'},...
%     'Labels',{'Different Region','Same Region'})
% title(sprintf('Enhancer Activity Residuals \n p = %.2g',pval_enhancer_activity))
% grid on

%% Distributions
promoter_activity_lm.Residuals{group_idx == 3,'Raw'}