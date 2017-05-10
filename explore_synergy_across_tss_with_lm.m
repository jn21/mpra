mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');
isfinite_idx = isfinite(mpra_data{:,'enhancer_activity'}) & isfinite(mpra_data{:,'promoter_activity'});
mpra_data = mpra_data(isfinite_idx,:);

%% Run Regressions
[lm_table, mpra_table_for_lm] = convert_data_table_to_lm_table(mpra_data,true,true);
enhancer_activity_lm = lm_table;
enhancer_activity_lm.promoter_activity = [];
promoter_activity_lm = lm_table;
promoter_activity_lm.enhancer_activity = [];

promoter_activity_lm = fitlm(promoter_activity_lm,...
            'ResponseVar','promoter_activity');
        
% promoter_activity_lm_no_intercept = fitlm(promoter_activity_lm,...
%             'ResponseVar','E_ratio',...
%             'Intercept',false);

enhancer_activity_lm = fitlm(enhancer_activity_lm,...
    'ResponseVar','enhancer_activity');

%% Make Plots Comparing Residuals
same_region_idx = logical(mpra_table_for_lm{:,'is_intact_sequence'});
promoter_derived_idx = strcmp(mpra_table_for_lm{:,'upstream_sequence_derivation'},'promoter_derived');

%group_idx
% 1 : different region
% 2 : same region - enhancer derived
% 3 : same region - promoter derived
group_idx = 1 + same_region_idx .* (1 + promoter_derived_idx);

p_promoter_diff = signtest(promoter_activity_lm.Residuals{group_idx == 1,'Raw'})
p_promoter_enhancer_derived = signtest(promoter_activity_lm.Residuals{group_idx == 2,'Raw'})
p_promoter_promoter_derived = signtest(promoter_activity_lm.Residuals{group_idx == 3,'Raw'})

p_enhancer_diff = signtest(enhancer_activity_lm.Residuals{group_idx == 1,'Raw'})
p_enhancer_enhancer_derived = signtest(enhancer_activity_lm.Residuals{group_idx == 2,'Raw'})
p_enhancer_promoter_derived = signtest(enhancer_activity_lm.Residuals{group_idx == 3,'Raw'})

pval_promoter_12 = ranksum(promoter_activity_lm.Residuals{group_idx == 1,'Raw'},...
        promoter_activity_lm.Residuals{group_idx == 2,'Raw'})
pval_promoter_13 = ranksum(promoter_activity_lm.Residuals{group_idx == 1,'Raw'},...
        promoter_activity_lm.Residuals{group_idx == 3,'Raw'})
pval_promoter_23 = ranksum(promoter_activity_lm.Residuals{group_idx == 2,'Raw'},...
        promoter_activity_lm.Residuals{group_idx == 3,'Raw'})
% 
% pval_enhancer_12 = ranksum(enhancer_activity_lm.Residuals{group_idx == 1,'Raw'},...
%         enhancer_activity_lm.Residuals{group_idx == 2,'Raw'});
% pval_enhancer_13 = ranksum(enhancer_activity_lm.Residuals{group_idx == 1,'Raw'},...
%         enhancer_activity_lm.Residuals{group_idx == 3,'Raw'});
% pval_enhancer_23 = ranksum(enhancer_activity_lm.Residuals{group_idx == 2,'Raw'},...
%         enhancer_activity_lm.Residuals{group_idx == 3,'Raw'});
    
% pval_enhancer_activity = ranksum(enhancer_activity_lm.Residuals{same_region_idx,'Raw'},...
%         enhancer_activity_lm.Residuals{~same_region_idx,'Raw'});

% promoter_activity_effect = median(promoter_activity_lm.Residuals{same_region_idx,'Raw'}) - ...
%     median(promoter_activity_lm.Residuals{~same_region_idx,'Raw'});
% 
% enhancer_activity_effect = median(enhancer_activity_lm.Residuals{same_region_idx,'Raw'}) - ...
%     median(enhancer_activity_lm.Residuals{~same_region_idx,'Raw'});


labels = {...
    sprintf('Different Region \nn = %d', nnz(group_idx == 1)),...
    sprintf('Same Region \nEnhancer Derived\nn = %d', nnz(group_idx == 2)),...
    sprintf('Same Region \nPromoter Derived\nn = %d', nnz(group_idx == 3))};

figure
subplot(1,2,1)
boxplot(promoter_activity_lm.Residuals{:,'Raw'},...
    group_idx)
title_str = sprintf('Promoter Activity \n R^2 = %.2f',...
    promoter_activity_lm.Rsquared.Ordinary);
title(title_str)
ylabel('Residuals')
grid on
ax = gca;
ax.XTickLabel = [];
text(ax.XTick - .35,repmat(ax.YLim(1),3,1) - .3,labels,'FontSize',12)
%H = sigstar({[1,3], [2,3]},[ pval_promoter_13, pval_promoter_23]);

subplot(1,2,2)
boxplot(enhancer_activity_lm.Residuals{:,'Raw'},...
    group_idx)
title_str = sprintf('Enhancer Activity \n R^2 = %.2f',...
    enhancer_activity_lm.Rsquared.Ordinary);
title(title_str)
grid on
ylabel('Residuals')
ax = gca;
ax.XTickLabel = [];
text(ax.XTick - .35,repmat(ax.YLim(1),3,1) - .1,labels,'FontSize',12)

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
figure
scatter(promoter_activity_lm.Fitted(group_idx == 1),...
    promoter_activity_lm.Residuals{group_idx == 1,'Raw'},...
    'DisplayName','Different Region')
hold on
scatter(promoter_activity_lm.Fitted(group_idx == 2),...
    promoter_activity_lm.Residuals{group_idx == 2,'Raw'},...
    'filled',...
    'DisplayName','Same Region - Enhancer Derived')
hold on
scatter(promoter_activity_lm.Fitted(group_idx == 3),...
    promoter_activity_lm.Residuals{group_idx == 3,'Raw'},...
    'filled',...
    'DisplayName','Same Region - Promoter Derived')
legend show
xlabel('Fitted')
ylabel('Residual')
title('Promoter Activity')
saveas(gcf,'~/Documents/mpra/fig/tss_synergy/fitted_vs_residual_promoter_activity','png')

figure
scatter(enhancer_activity_lm.Fitted(group_idx == 1),...
    enhancer_activity_lm.Residuals{group_idx == 1,'Raw'},...
    'DisplayName','Different Region')
hold on
scatter(enhancer_activity_lm.Fitted(group_idx == 2),...
    enhancer_activity_lm.Residuals{group_idx == 2,'Raw'},...
    'filled',...
    'DisplayName','Same Region - Enhancer Derived')
hold on
scatter(enhancer_activity_lm.Fitted(group_idx == 3),...
    enhancer_activity_lm.Residuals{group_idx == 3,'Raw'},...
    'filled',...
    'DisplayName','Same Region - Promoter Derived')
legend show
xlabel('Fitted')
ylabel('Residual')
title('Enhancer Activity')
saveas(gcf,'~/Documents/mpra/fig/tss_synergy/fitted_vs_residual_enhancer_activity','png')
