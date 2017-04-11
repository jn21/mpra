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
pval_promoter_activity = ranksum(promoter_activity_lm.Residuals{same_region_idx,'Raw'},...
        promoter_activity_lm.Residuals{~same_region_idx,'Raw'});
    
promoter_activity_effect = median(promoter_activity_lm.Residuals{same_region_idx,'Raw'}) - ...
    median(promoter_activity_lm.Residuals{~same_region_idx,'Raw'});

enhancer_activity_effect = median(enhancer_activity_lm.Residuals{same_region_idx,'Raw'}) - ...
    median(enhancer_activity_lm.Residuals{~same_region_idx,'Raw'});

pval_enhancer_activity = ranksum(enhancer_activity_lm.Residuals{same_region_idx,'Raw'},...
        enhancer_activity_lm.Residuals{~same_region_idx,'Raw'});

figure
subplot(1,2,1)
boxplot(promoter_activity_lm.Residuals{:,'Raw'},...
    mpra_table_for_lm{:,'is_intact_sequence'},...
    'Labels',{'Different Region','Same Region'})
title(sprintf('Promoter Activity \n p = %.2g',pval_promoter_activity))
grid on

subplot(1,2,2)
boxplot(enhancer_activity_lm.Residuals{:,'Raw'},...
    mpra_table_for_lm{:,'is_intact_sequence'},...
    'Labels',{'Different Region','Same Region'})
title(sprintf('Enhancer Activity \n p = %.2g',pval_enhancer_activity))
grid on