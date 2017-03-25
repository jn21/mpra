mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');
prefix_annot = readtable('~/Documents/mpra/data/prefix_is_from_promoter_annot.txt','Delimiter','\t');

%% Remove some low 'E' count entries
CUTOFF = 20;
rep1_idx = (mpra_data{:,'Rep1_ETotal'} < CUTOFF);
rep2_idx = (mpra_data{:,'Rep2_ETotal'} < CUTOFF);
mpra_data = mpra_data(~(rep1_idx | rep2_idx),:);

%% Linear regression
lm_table = convert_data_table_to_lm_table(mpra_data);

pratio_table = lm_table; 
pratio_table.E_ratio = [];

eratio_table = lm_table; 
eratio_table.P_ratio = [];

pratio_lm = fitlm(pratio_table,'ResponseVar','P_ratio');
eratio_lm = fitlm(eratio_table,'ResponseVar','E_ratio');

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
    
    this_name = coeff_table{ii,'name'};
    splits = strsplit(this_name{1},'_');
    
    if length(splits) > 1 && strcmp(splits{2},'prefix')
        coeff_table{ii,'tag'} = splits(1);
        
        %prefix annot
        idx = find(prefix_annot{:,'prefix'} == str2double(splits{3}));
        coeff_table{ii,'is_within_100nt_of_tss'} = prefix_annot{idx,'is_within_100nt_of_tss'};
    else
        coeff_table{ii,'tag'} = this_name;
        coeff_table{ii,'is_within_100nt_of_tss'} = NaN;
    end
     
end

%% Write table
writetable(coeff_table,'~/Documents/mpra/data/temp_coeff_table.txt','Delimiter','\t');

%% Plotting
scatter(eratio_lm.Coefficients{:,'Estimate'},pratio_lm.Coefficients{:,'Estimate'})
xlabel('construct_ratio','Interpreter','none')
ylabel('GFP_ratio','Interpreter','none')

%% Residuals