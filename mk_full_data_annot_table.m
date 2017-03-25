mpra_data = readtable('~/Documents/mpra/data/insert_totals.txt','Delimiter','\t');

construct_names = mpra_data{:,'construct'};

construct_names_parsed = parse_construct_name(construct_names);

T = join(mpra_data,construct_names_parsed,'Keys','construct');

writetable(T,'~/Documents/mpra/data/insert_totals_with_annot.txt','Delimiter','\t');

%% Do some basic calculations
raw_data = readtable('~/Documents/mpra/data/insert_totals_with_annot.txt',...
    'Delimiter','\t');

raw_data{:,'rep1_E_ratio'} = log2(raw_data{:,'Rep1_ETotal'} ./ raw_data{:,'DNAInput_ETotal'});
raw_data{:,'rep2_E_ratio'} = log2(raw_data{:,'Rep2_ETotal'} ./ raw_data{:,'DNAInput_ETotal'});

raw_data{:,'E_ratio_avg_rep'} = (raw_data{:,'rep1_E_ratio'} + raw_data{:,'rep2_E_ratio'})/2;


raw_data{:,'rep1_P_ratio'} = log2(raw_data{:,'Rep1_PTotal'} ./ raw_data{:,'DNAInput_PTotal'});
raw_data{:,'rep2_P_ratio'} = log2(raw_data{:,'Rep2_PTotal'} ./ raw_data{:,'DNAInput_PTotal'});

raw_data{:,'P_ratio_avg_rep'} = (raw_data{:,'rep1_P_ratio'} + raw_data{:,'rep2_P_ratio'})/2;

writetable(raw_data,'~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');
