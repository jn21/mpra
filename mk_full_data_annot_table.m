mpra_data = readtable('~/Documents/mpra/data/insert_totals.txt','Delimiter','\t');

construct_names = mpra_data{:,'construct'};

% construct_names_parsed = parse_construct_name(construct_names);
construct_names_parsed = fast_parse_construct_name(construct_names);

T = join(mpra_data,construct_names_parsed,'Keys','construct');

writetable(T,'~/Documents/mpra/data/insert_totals_with_annot.txt','Delimiter','\t');

raw_data = readtable('~/Documents/mpra/data/insert_totals_with_annot.txt',...
    'Delimiter','\t');

%% Add some combined annotations
raw_data{:,'dnstream_is_modified'} = raw_data{:,'dnstream_addPAS'} | ...
                                     raw_data{:,'dnstream_addStrongPAS'} | ...
                                     raw_data{:,'dnstream_num_delU1'} | ...
                                     raw_data{:,'dnstream_addU1'} | ...
                                     raw_data{:,'dnstream_delPAS'};

raw_data{:,'is_intact_sequence'} = ~logical(raw_data{:,'dnstream_is_modified'}) &...
                                   ~logical(raw_data{:,'dnstream_is_reverse'}) &...
                                   ~logical(raw_data{:,'upstream_is_reverse'}) &...
                                   raw_data{:,'upstream_region_id'} == raw_data{:,'dnstream_region_id'};
                               
%% Do some basic calculations                               
raw_data{:,'rep1_A_ratio'} = log2(raw_data{:,'Rep1_ATotal'} ./ raw_data{:,'DNAInput_ATotal'});
raw_data{:,'rep2_A_ratio'} = log2(raw_data{:,'Rep2_ATotal'} ./ raw_data{:,'DNAInput_ATotal'});

raw_data{:,'promoter_activity'} = (raw_data{:,'rep1_A_ratio'} + raw_data{:,'rep2_A_ratio'})/2;

raw_data{:,'rep1_B_ratio'} = log2(raw_data{:,'Rep1_BTotal'} ./ raw_data{:,'DNAInput_BTotal'});
raw_data{:,'rep2_B_ratio'} = log2(raw_data{:,'Rep2_BTotal'} ./ raw_data{:,'DNAInput_BTotal'});

raw_data{:,'enhancer_activity'} = (raw_data{:,'rep1_B_ratio'} + raw_data{:,'rep2_B_ratio'})/2;


                               
%% Write result
writetable(raw_data,'~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');
