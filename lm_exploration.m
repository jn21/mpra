mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');
prefix_annot = readtable('~/Documents/mpra/data/prefix_is_from_promoter_annot.txt','Delimiter','\t');

%% Remove some low 'E' count entries
% CUTOFF = 20;
% rep1_idx = (mpra_data{:,'Rep1_ETotal'} < CUTOFF);
% rep2_idx = (mpra_data{:,'Rep2_ETotal'} < CUTOFF);
% mpra_data = mpra_data(~(rep1_idx | rep2_idx),:);

%% Remove non-finite entries
idx = isfinite(mpra_data{:,'P_ratio_avg_rep'}) & isfinite(mpra_data{:,'E_ratio_avg_rep'});
mpra_data = mpra_data(idx,:);

%% Linear regression
lm_table = convert_data_table_to_lm_table(mpra_data);

pratio_table = lm_table; 
pratio_table.E_ratio = [];

eratio_table = lm_table; 
eratio_table.P_ratio = [];

% pratio_lm = fitlm(pratio_table,'ResponseVar','P_ratio');
% eratio_lm = fitlm(eratio_table,'ResponseVar','E_ratio');

%% ANOVA
[p_anova,anova_tab,anova_stats] = anova1(mpra_data{:,'E_ratio_avg_rep'},...
    mpra_data{:,'dnstream_full_id'},...
    'off');

%% Run linear regression with only upstream/dnstream ids
% tag = 'up';
% tag_idx = strcmp(tag,...
%         cellfun(@(s) s(1:2),pratio_table.Properties.VariableNames,'uni',false));
% pratio_lm_tag = fitlm(pratio_table,...
%     'ResponseVar','P_ratio',...
%     'PredictorVars',pratio_table.Properties.VariableNames(tag_idx),...
%     'Intercept',false);

[lm_table, mpra_table_for_lm] = convert_data_table_to_lm_table(mpra_data,true);
enhancer_activity = lm_table;
enhancer_activity.E_ratio = [];
promoter_activity = lm_table;
promoter_activity.P_ratio = [];

promoter_activity_lm = fitlm(promoter_activity,...
            'ResponseVar','E_ratio');
        
enhancer_activity_lm = fitlm(enhancer_activity,...
    'ResponseVar','P_ratio');

%%
figure
subplot(1,2,1)
boxplot(promoter_activity_lm.Residuals{:,'Raw'},mpra_table_for_lm{:,'is_intact_sequence'})
title('Promoter Activity')
grid on

subplot(1,2,2)
boxplot(enhancer_activity_lm.Residuals{:,'Raw'},mpra_table_for_lm{:,'is_intact_sequence'})
title('Enhancer Activity')
grid on



%% Make coefficients table
assert(isequal(eratio_lm.Coefficients.Properties.RowNames,pratio_lm.Coefficients.Properties.RowNames),...
    'row ordering not the same');
coeff_table = table(pratio_lm.CoefficientNames',...
    pratio_lm.Coefficients{:,'Estimate'},...
    eratio_lm.Coefficients{:,'Estimate'},...
    'VariableNames',{'name','pratio_coeff','eratio_coeff'});

%% add some metadata - eg if each sequence is from a genomic 'promoter' or 'enhancer'
coeff_name = coeff_table{:,'name'};
for ii = 1:length(coeff_name)
    ii
    
    this_name = coeff_table{ii,'name'};
    splits = strsplit(this_name{1},'_');
    
    if length(splits) > 1 && ismember(splits{1},{'up','dn'})
        coeff_table{ii,'tag'} = splits(1);
        
        %prefix annot
        idx = find(prefix_annot{:,'prefix'} == str2double(splits{2}));
        coeff_table{ii,'is_within_100nt_of_tss'} = prefix_annot{idx,'is_within_100nt_of_tss'};
    else
        coeff_table{ii,'tag'} = this_name;
        coeff_table{ii,'is_within_100nt_of_tss'} = NaN;
    end
    
%     if strcmp(splits{end},'reverse')
%         coeff_table{ii,'is_reverse'} = 1
        
    coeff_table{ii,'is_reverse'} = logical(strcmp(splits{end},'reverse'));
     
end

%% Write table
writetable(coeff_table,'~/Documents/mpra/data/temp_coeff_table.txt','Delimiter','\t');

%% Plotting
scatter(eratio_lm.Coefficients{:,'Estimate'},pratio_lm.Coefficients{:,'Estimate'})
xlabel('construct_ratio','Interpreter','none')
ylabel('GFP_ratio','Interpreter','none')

%% Get correlations of coeffs
up_idx = strcmp(coeff_table{:,'tag'},'up');
[up_corr, up_pval] = corr(coeff_table{up_idx,'pratio_coeff'},coeff_table{up_idx,'eratio_coeff'})


dn_idx = strcmp(coeff_table{:,'tag'},'dn');
[dn_corr, dn_pval] = corr(coeff_table{dn_idx,'pratio_coeff'},coeff_table{dn_idx,'eratio_coeff'})

%% Examine coefficients for wt_100nt/gt_10kb enrichment

p_or_e_ratio = 'eratio_coeff';
up_or_dn = 'up';

up_idx = strcmp(coeff_table{:,'tag'},'up');
dn_idx = strcmp(coeff_table{:,'tag'},'dn');
is_within_idx = coeff_table{:,'is_within_100nt_of_tss'};
is_reverse_idx = coeff_table{:,'is_reverse'};

switch up_or_dn
    case 'up'
        up_or_dn_idx = up_idx;
    case 'dn'
        up_or_dn_idx = dn_idx;
end

rs = ranksum(coeff_table{up_or_dn_idx & (is_within_idx == 1) & ~is_reverse_idx,p_or_e_ratio},coeff_table{up_or_dn_idx & (is_within_idx == 0) & ~is_reverse_idx,p_or_e_ratio});
boxplot(coeff_table{up_or_dn_idx,p_or_e_ratio},coeff_table{up_or_dn_idx,'is_within_100nt_of_tss'},...
    'Labels',{'gt_10kb','wt_100nt'})
title_str = sprintf('which_barcode = %s, up_or_dn_coeff = %s \n ranksum pval = %.3g',p_or_e_ratio,up_or_dn,rs);
title(title_str,'Interpreter','none')

